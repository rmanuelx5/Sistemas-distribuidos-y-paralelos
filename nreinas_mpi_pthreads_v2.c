#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <sys/time.h>
#include <mpi.h>

#define MAXSIZE 24
#define MINSIZE 2

double dwalltime();

int SIZE, SIZEE;
int MASK, TOPBIT;

int T;          /* pthreads por proceso MPI */
int MPI_RANK;
int MPI_SIZE;

/*
 * Trabajo distribuido.
 * Cada entrada empaqueta dos columnas en un int:
 *   bits 16..31 -> colA  (col fila 0 o fila 1 segun la fase)
 *   bits  0..15 -> colB  (col de la fila siguiente)
 * -1 indica celda de padding.
 *
 * Todas las tareas de ambas fases se concentran en un unico arreglo
 * y se distribuyen con un solo MPI_Scatter.
 */
int *WORK_BOUNDS       = NULL;
int  WORK_COUNT        = 0;   /* entradas que recibe cada proceso */
int  WORK_INDEX        = 0;   /* proximo indice a tomar (protegido por mutex) */
int  WORK_GLOBAL_OFFSET = 0;  /* indice global del primer elemento de este proceso */
int  FASE1_COUNT       = 0;   /* cuantas tareas globales pertenecen a la fase 1 */

pthread_mutex_t mutex;

long int THREAD_COUNT8[MAXSIZE];
long int THREAD_COUNT4[MAXSIZE];
long int THREAD_COUNT2[MAXSIZE];
double   tiempo_hilos[MAXSIZE];

/* ---------- helpers de empaquetado ---------- */
static inline int pack(int a, int b)  { return (a << 16) | (b & 0xFFFF); }
static inline int unpack0(int v)      { return (v >> 16) & 0xFFFF; }
static inline int unpack1(int v)      { return  v        & 0xFFFF; }

/* ---------- estructura por hilo ---------- */
typedef struct {
    int BOARD[MAXSIZE];
    int *BOARDE;
    int *BOARD1;
    int *BOARD2;
    int SIDEMASK;
    int LASTMASK;
    int ENDBIT;
    int BOUND1;
    int BOUND2;
    long int COUNT8;
    long int COUNT4;
    long int COUNT2;
} structHilo;

/* ---------- Check y Backtrack sin modificar ---------- */

void Check(structHilo *s)
{
    int *own, *you, bit, ptn;

    if (*(s->BOARD2) == 1) {
        for (ptn=2, own=s->BOARD+1; own<=s->BOARDE; own++, ptn<<=1) {
            bit = 1;
            for (you=s->BOARDE; *you!=ptn && *own>=bit; you--)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
        if (own > s->BOARDE) { s->COUNT2++; return; }
    }

    if (*(s->BOARDE) == s->ENDBIT) {
        for (you=s->BOARDE-1, own=s->BOARD+1; own<=s->BOARDE; own++, you--) {
            bit = 1;
            for (ptn=TOPBIT; ptn!=*you && *own>=bit; ptn>>=1)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
        if (own > s->BOARDE) { s->COUNT4++; return; }
    }

    if (*(s->BOARD1) == TOPBIT) {
        for (ptn=TOPBIT>>1, own=s->BOARD+1; own<=s->BOARDE; own++, ptn>>=1) {
            bit = 1;
            for (you=s->BOARD; *you!=ptn && *own>=bit; you++)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
    }
    s->COUNT8++;
}

void Backtrack2(int y, int left, int down, int right, structHilo *s)
{
    int bitmap, bit;
    bitmap = MASK & ~(left | down | right);
    if (y == SIZEE) {
        if (bitmap) {
            if (!(bitmap & s->LASTMASK)) {
                s->BOARD[y] = bitmap;
                Check(s);
            }
        }
    } else {
        if (y < s->BOUND1) {
            bitmap |= s->SIDEMASK;
            bitmap ^= s->SIDEMASK;
        } else if (y == s->BOUND2) {
            if (!(down & s->SIDEMASK)) return;
            if ((down & s->SIDEMASK) != s->SIDEMASK) bitmap &= s->SIDEMASK;
        }
        while (bitmap) {
            bitmap ^= s->BOARD[y] = bit = -bitmap & bitmap;
            Backtrack2(y+1, (left|bit)<<1, down|bit, (right|bit)>>1, s);
        }
    }
}

void Backtrack1(int y, int left, int down, int right, structHilo *s)
{
    int bitmap, bit;
    bitmap = MASK & ~(left | down | right);
    if (y == SIZEE) {
        if (bitmap) {
            s->BOARD[y] = bitmap;
            s->COUNT8++;
        }
    } else {
        if (y < s->BOUND1) {
            bitmap |= 2;
            bitmap ^= 2;
        }
        while (bitmap) {
            bitmap ^= s->BOARD[y] = bit = -bitmap & bitmap;
            Backtrack1(y+1, (left|bit)<<1, down|bit, (right|bit)>>1, s);
        }
    }
}

/* ---------- Worker ---------- */

void *worker(void *arg){
    int id = *(int *)arg;
    tiempo_hilos[id] = dwalltime();

    pthread_mutex_lock(&mutex);
    while (WORK_INDEX < WORK_COUNT) {
        int local_idx = WORK_INDEX++;
        pthread_mutex_unlock(&mutex);

        int packed = WORK_BOUNDS[local_idx];
        if (packed == -1) {
            /* celda de padding, ignorar */
            pthread_mutex_lock(&mutex);
            continue;
        }

        /* La fase se determina por la posicion global de la tarea */
        int global_idx = WORK_GLOBAL_OFFSET + local_idx;
        int phase = (global_idx < FASE1_COUNT) ? 1 : 2;

        structHilo s = {0};
        s.BOARDE = &s.BOARD[SIZEE];

        if (phase == 1) {
            /*
             * Fase 1: primera reina en la esquina (fila 0, columna 0).
             * packed = (col1, col2): columnas para filas 1 y 2.
             * BOUND1 = col1 (restriccion de simetria original).
             *
             * Fijamos tres reinas y llamamos a Backtrack1 desde fila 3.
             */
            int col1 = unpack0(packed);
            int col2 = unpack1(packed);

            int bit0 = 1;           /* fila 0 */
            int bit1 = 1 << col1;  /* fila 1 */
            int bit2 = 1 << col2;  /* fila 2 */

            s.BOUND1   = col1;
            s.BOARD[0] = bit0;
            s.BOARD[1] = bit1;
            s.BOARD[2] = bit2;

            /*
             * Propagacion de diagonales:
             *   Al entrar a fila y+1: L = (L_prev | bit_y) << 1
             *                          D = D_prev | bit_y
             *                          R = (R_prev | bit_y) >> 1
             */
            int L1 = (bit0 << 1) & MASK;
            int D1 =  bit0;
            int R1 = (bit0 >> 1) & MASK;

            int L2 = ((L1 | bit1) << 1) & MASK;
            int D2 =   D1 | bit1;
            int R2 = ((R1 | bit1) >> 1) & MASK;

            int L3 = ((L2 | bit2) << 1) & MASK;
            int D3 =   D2 | bit2;
            int R3 = ((R2 | bit2) >> 1) & MASK;

            Backtrack1(3, L3, D3, R3, &s);

        } else {
            /*
             * Fase 2: primera reina en el interior.
             * packed = (col0, col1): columnas para filas 0 y 1.
             * BOUND1 = col0, BOUND2 = SIZE-1-col0 (igual que antes).
             *
             * Fijamos dos reinas y llamamos a Backtrack2 desde fila 2.
             */
            int col0 = unpack0(packed);
            int col1 = unpack1(packed);

            int bit0 = 1 << col0;
            int bit1 = 1 << col1;

            s.BOUND1   = col0;
            s.BOUND2   = SIZE - 1 - col0;
            s.BOARD1   = &s.BOARD[s.BOUND1];
            s.BOARD2   = &s.BOARD[s.BOUND2];
            s.SIDEMASK = TOPBIT | 1;
            s.LASTMASK = TOPBIT | 1;
            s.ENDBIT   = TOPBIT >> 1;

            for (int i = 1; i < s.BOUND1; i++) {
                s.LASTMASK |= (s.LASTMASK >> 1) | (s.LASTMASK << 1);
                s.ENDBIT >>= 1;
            }

            s.BOARD[0] = bit0;
            s.BOARD[1] = bit1;

            int L1 = (bit0 << 1) & MASK;
            int D1 =  bit0;
            int R1 = (bit0 >> 1) & MASK;

            int L2 = ((L1 | bit1) << 1) & MASK;
            int D2 =   D1 | bit1;
            int R2 = ((R1 | bit1) >> 1) & MASK;

            Backtrack2(2, L2, D2, R2, &s);
        }

        THREAD_COUNT8[id] += s.COUNT8;
        THREAD_COUNT4[id] += s.COUNT4;
        THREAD_COUNT2[id] += s.COUNT2;

        pthread_mutex_lock(&mutex);
    }
    pthread_mutex_unlock(&mutex);

    tiempo_hilos[id] = dwalltime() - tiempo_hilos[id];
    return NULL;
}

/* ---------- Logica de hilos ---------- */

void logica_hilos(long int *lc8, long int *lc4, long int *lc2)
{
    pthread_t threads[T];
    for (int i = 0; i < T; i++)
        THREAD_COUNT8[i] = THREAD_COUNT4[i] = THREAD_COUNT2[i] = 0;

    pthread_mutex_init(&mutex, NULL);

    int ids[T];
    for (int i = 0; i < T; i++) {
        ids[i] = i;
        pthread_create(&threads[i], NULL, worker, &ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
        *lc8 += THREAD_COUNT8[i];
        *lc4 += THREAD_COUNT4[i];
        *lc2 += THREAD_COUNT2[i];
    }
    pthread_mutex_destroy(&mutex);
}

/* ---------- Generacion de tareas ---------- */

/*
 * Solo rank 0 llama a esta funcion.
 * Genera todas las tareas de fase 1 y fase 2 en un unico arreglo.
 * La validez de cada posicion se verifica propagando las diagonales,
 * igual que lo hace el worker, para no generar combinaciones imposibles.
 *
 * Fase 1: bit0=1 (esquina). Para cada col1 en [2, SIZEE-1] y cada
 *         col2 en [0, SIZE-1] que no colisione con bit0 y bit1.
 *
 * Fase 2: col0 en [1, (SIZE-2)/2]. Para cada col1 en [0, SIZE-1]
 *         que no colisione con bit0.
 */
static int *generar_tareas(int *total, int *f1_count)
{
    int cap = SIZE * SIZE * 2 + 16;
    int *buf = malloc(sizeof(int) * cap);
    int n = 0;

    /* --- Fase 1 --- */
    for (int col1 = 2; col1 <= SIZEE - 1; col1++) {
        int bit0 = 1;
        int bit1 = 1 << col1;

        int L1 = (bit0 << 1) & MASK;
        int D1 =  bit0;
        int R1 = (bit0 >> 1) & MASK;

        int L2 = ((L1 | bit1) << 1) & MASK;
        int D2 =   D1 | bit1;
        int R2 = ((R1 | bit1) >> 1) & MASK;
        int forbidden2 = L2 | D2 | R2;

        for (int col2 = 0; col2 < SIZE; col2++) {
            int bit2 = 1 << col2;
            if (bit2 & forbidden2) continue;
            buf[n++] = pack(col1, col2);
        }
    }
    *f1_count = n;

    /* --- Fase 2 --- */
    for (int col0 = 1; col0 <= (SIZE - 2) / 2; col0++) {
        int bit0 = 1 << col0;

        int L1 = (bit0 << 1) & MASK;
        int D1 =  bit0;
        int R1 = (bit0 >> 1) & MASK;
        int forbidden1 = L1 | D1 | R1;

        for (int col1 = 0; col1 < SIZE; col1++) {
            int bit1 = 1 << col1;
            if (bit1 & forbidden1) continue;
            buf[n++] = pack(col0, col1);
        }
    }
    *total = n;

    return buf;
}

/* ---------- Preparar y distribuir trabajo (unica comunicacion) ---------- */

void preparar_y_distribuir_trabajo(void)
{
    int total_tareas = 0;
    int fase1_count  = 0;
    int *todas       = NULL;

    if (MPI_RANK == 0)
        todas = generar_tareas(&total_tareas, &fase1_count);

    /*
     * Difundir los dos contadores. Todos los nodos pueden calcular
     * WORK_COUNT a partir de ellos sin comunicacion extra.
     */
    int info[2] = {total_tareas, fase1_count};
    MPI_Bcast(info, 2, MPI_INT, 0, MPI_COMM_WORLD);
    total_tareas = info[0];
    fase1_count  = info[1];
    FASE1_COUNT  = fase1_count;

    WORK_COUNT = (total_tareas + MPI_SIZE - 1) / MPI_SIZE;
    int total_celdas = WORK_COUNT * MPI_SIZE;

    int *buf_scatter = NULL;
    if (MPI_RANK == 0) {
        buf_scatter = malloc(sizeof(int) * total_celdas);
        for (int i = 0; i < total_celdas; i++)
            buf_scatter[i] = (i < total_tareas) ? todas[i] : -1;
        free(todas);
    }

    if (WORK_BOUNDS != NULL) { free(WORK_BOUNDS); WORK_BOUNDS = NULL; }
    WORK_BOUNDS = malloc(sizeof(int) * WORK_COUNT);

    MPI_Scatter(buf_scatter, WORK_COUNT, MPI_INT,
                WORK_BOUNDS,  WORK_COUNT, MPI_INT,
                0, MPI_COMM_WORLD);

    if (MPI_RANK == 0) free(buf_scatter);

    WORK_GLOBAL_OFFSET = MPI_RANK * WORK_COUNT;
    WORK_INDEX = 0;
}

/* ---------- NQueens ---------- */

double NQueens(void)
{
    long int LOCAL_COUNT8  = 0, LOCAL_COUNT4  = 0, LOCAL_COUNT2  = 0;
    long int GLOBAL_COUNT8 = 0, GLOBAL_COUNT4 = 0, GLOBAL_COUNT2 = 0;
    double sumaLocal = 0.0, sumaGlobal = 0.0;
    double tIni, tFin;
    int totalThreads = MPI_SIZE * T;

    SIZEE  = SIZE - 1;
    TOPBIT = 1 << SIZEE;
    MASK   = (1 << SIZE) - 1;

    tIni = dwalltime();

    preparar_y_distribuir_trabajo();
    logica_hilos(&LOCAL_COUNT8, &LOCAL_COUNT4, &LOCAL_COUNT2);

    /* Un unico reduce con arreglo de 3 elementos */
    long int local_counts[3]  = {LOCAL_COUNT8, LOCAL_COUNT4, LOCAL_COUNT2};
    long int global_counts[3] = {0, 0, 0};
    MPI_Reduce(local_counts, global_counts, 3, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    tFin = dwalltime();

    if (MPI_RANK == 0) {
        GLOBAL_COUNT8 = global_counts[0];
        GLOBAL_COUNT4 = global_counts[1];
        GLOBAL_COUNT2 = global_counts[2];
        long int UNIQUE = GLOBAL_COUNT8 + GLOBAL_COUNT4 + GLOBAL_COUNT2;
        long int TOTAL  = GLOBAL_COUNT8 * 8 + GLOBAL_COUNT4 * 4 + GLOBAL_COUNT2 * 2;
        printf("Numero de resultados: %lu - Unicas: %lu\n", TOTAL, UNIQUE);
    }

    /* Tiempos por hilo (fuera de la toma de tiempo principal) */
    for (int r = 0; r < MPI_SIZE; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (MPI_RANK == r) {
            for (int i = 0; i < T; i++) {
                sumaLocal += tiempo_hilos[i];
                printf("Rank %d, Thread %d: Tiempo = %f segundos\n",
                       MPI_RANK, i, tiempo_hilos[i]);
            }
        }
    }
    MPI_Reduce(&sumaLocal, &sumaGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (MPI_RANK == 0)
        printf("Tiempo promedio total por hilo = %f segundos\n",
               sumaGlobal / totalThreads);

    return tFin - tIni;
}

/* ---------- Main ---------- */

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

    if (argc < 2) {
        if (MPI_RANK == 0)
            fprintf(stderr, "Uso: %s N [pthreads_por_proceso]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    SIZE = atoi(argv[1]);
    T    = (argc > 2) ? atoi(argv[2]) : 1;

    if (SIZE < MINSIZE || SIZE > MAXSIZE) {
        if (MPI_RANK == 0)
            fprintf(stderr, "N debe estar entre %d y %d\n", MINSIZE, MAXSIZE);
        MPI_Finalize();
        return 1;
    }

    if (T < 1)       T = 1;
    if (T > MAXSIZE) T = MAXSIZE;

    MPI_Barrier(MPI_COMM_WORLD);

    double tiempo = NQueens();

    MPI_Barrier(MPI_COMM_WORLD);

    if (WORK_BOUNDS != NULL) { free(WORK_BOUNDS); WORK_BOUNDS = NULL; }

    if (MPI_RANK == 0) {
        printf("Tamano del tablero: %d\n", SIZE);
        printf("Procesos MPI: %d\n", MPI_SIZE);
        printf("pthreads por proceso MPI: %d\n", T);
        printf("Total de workers: %d\n", MPI_SIZE * T);
        printf("Tiempo Total: %f segundos\n", tiempo);
    }

    MPI_Finalize();
    return 0;
}

double dwalltime()
{
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}
