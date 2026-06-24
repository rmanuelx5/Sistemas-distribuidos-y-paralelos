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

// pthreads por cada proceso MPI
int T;

// datos MPI
int MPI_RANK;
int MPI_SIZE;

// trabajo recibido por cada proceso
int *WORK_BOUNDS = NULL;
int WORK_COUNT   = 0;
int WORK_INDEX   = 0;
int PHASE;

// exclusion mutua para tomar trabajo
pthread_mutex_t mutex;

// contadores locales por hilo
long int THREAD_COUNT8[MAXSIZE];
long int THREAD_COUNT4[MAXSIZE];
long int THREAD_COUNT2[MAXSIZE];

// tiempo de ejecucion de cada hilo
double tiempo_hilos[MAXSIZE];

// variables privadas de cada hilo
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

int tomar_trabajo(void)
{
    int bound;
    pthread_mutex_lock(&mutex);
    if (WORK_INDEX >= WORK_COUNT) {
        bound = -1;
    } else {
        bound = WORK_BOUNDS[WORK_INDEX];
        WORK_INDEX++;
    }
    pthread_mutex_unlock(&mutex);
    return bound;
}

void Display(structHilo *s)
{
    int y, bit;
    printf("N= %d\n", SIZE);
    for (y = 0; y < SIZE; y++) {
        for (bit = TOPBIT; bit; bit >>= 1)
            printf("%s ", (s->BOARD[y] & bit) ? "Q" : "-");
        printf("\n");
    }
    printf("\n");
}

void Check(structHilo *s)
{
    int *own, *you, bit, ptn;

    /* 90-degree rotation */
    if (*(s->BOARD2) == 1) {
        for (ptn = 2, own = s->BOARD + 1; own <= s->BOARDE; own++, ptn <<= 1) {
            bit = 1;
            for (you = s->BOARDE; *you != ptn && *own >= bit; you--)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
        if (own > s->BOARDE) { s->COUNT2++; return; }
    }

    /* 180-degree rotation */
    if (*(s->BOARDE) == s->ENDBIT) {
        for (you = s->BOARDE - 1, own = s->BOARD + 1; own <= s->BOARDE; own++, you--) {
            bit = 1;
            for (ptn = TOPBIT; ptn != *you && *own >= bit; ptn >>= 1)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
        if (own > s->BOARDE) { s->COUNT4++; return; }
    }

    /* 270-degree rotation */
    if (*(s->BOARD1) == TOPBIT) {
        for (ptn = TOPBIT >> 1, own = s->BOARD + 1; own <= s->BOARDE; own++, ptn >>= 1) {
            bit = 1;
            for (you = s->BOARD; *you != ptn && *own >= bit; you++)
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
            Backtrack2(y + 1, (left | bit) << 1, down | bit, (right | bit) >> 1, s);
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
            Backtrack1(y + 1, (left | bit) << 1, down | bit, (right | bit) >> 1, s);
        }
    }
}

void *worker(void *arg)
{
    int id = *(int *)arg;
    tiempo_hilos[id] = dwalltime();

    int bound = tomar_trabajo();

    while (bound != -1) {
        structHilo s = {0};
        int bit;

        s.BOARDE = &s.BOARD[SIZEE];
        s.BOUND1 = bound;

        if (PHASE == 1) {
            s.BOARD[0] = 1;
            s.BOARD[1] = bit = 1 << s.BOUND1;
            Backtrack1(2, (2 | bit) << 1, 1 | bit, bit >> 1, &s);
        } else {
            s.BOUND2     = SIZE - 1 - s.BOUND1;
            s.BOARD1     = &s.BOARD[s.BOUND1];
            s.BOARD2     = &s.BOARD[s.BOUND2];
            s.SIDEMASK   = TOPBIT | 1;
            s.LASTMASK   = TOPBIT | 1;
            s.ENDBIT     = TOPBIT >> 1;

            for (int i = 1; i < s.BOUND1; i++) {
                s.LASTMASK |= (s.LASTMASK >> 1) | (s.LASTMASK << 1);
                s.ENDBIT   >>= 1;
            }

            s.BOARD[0] = bit = 1 << s.BOUND1;
            Backtrack2(1, bit << 1, bit, bit >> 1, &s);
        }

        THREAD_COUNT8[id] += s.COUNT8;
        THREAD_COUNT4[id] += s.COUNT4;
        THREAD_COUNT2[id] += s.COUNT2;

        bound = tomar_trabajo();
    }

    tiempo_hilos[id] = dwalltime() - tiempo_hilos[id];
    return NULL;
}

void logica_hilos(long int *local_count8, long int *local_count4, long int *local_count2)
{
    pthread_t threads[T];

    for (int i = 0; i < T; i++) {
        THREAD_COUNT8[i] = 0;
        THREAD_COUNT4[i] = 0;
        THREAD_COUNT2[i] = 0;
    }

    pthread_mutex_init(&mutex, NULL);

    int threads_ids[T];
    for (int i = 0; i < T; i++) {
        threads_ids[i] = i;
        pthread_create(&threads[i], NULL, &worker, (void *)&threads_ids[i]);
    }

    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
        *local_count8 += THREAD_COUNT8[i];
        *local_count4 += THREAD_COUNT4[i];
        *local_count2 += THREAD_COUNT2[i];
    }

    pthread_mutex_destroy(&mutex);
}

/*
 * preparar_trabajo: distribuye los valores [inicio..fin] entre los procesos MPI
 * usando MPI_Scatterv para que cada proceso reciba exactamente su porcion,
 * sin padding ni trabajos fantasma marcados con -1.
 */
void preparar_trabajo(int phase, int inicio, int fin)
{
    PHASE      = phase;
    WORK_INDEX = 0;

    if (WORK_BOUNDS != NULL) {
        free(WORK_BOUNDS);
        WORK_BOUNDS = NULL;
    }

    int total = (fin >= inicio) ? (fin - inicio + 1) : 0;

    int *sendcounts = calloc(MPI_SIZE, sizeof(int));
    int *displs     = calloc(MPI_SIZE, sizeof(int));
    int *all        = NULL;

    /* Solo root construye el arreglo completo y los sendcounts reales */
    if (MPI_RANK == 0) {
        int base  = total / MPI_SIZE;
        int resto = total % MPI_SIZE;
        int offset = 0;

        all = malloc(sizeof(int) * (total > 0 ? total : 1));

        for (int p = 0; p < MPI_SIZE; p++) {
            sendcounts[p] = base + (p < resto ? 1 : 0);
            displs[p]     = offset;
            for (int i = 0; i < sendcounts[p]; i++)
                all[offset + i] = inicio + offset + i;
            offset += sendcounts[p];
        }
    }

    /*
     * Broadcast de sendcounts y displs para que cada proceso sepa
     * cuantos elementos va a recibir (necesario antes de Scatterv).
     */
    MPI_Bcast(sendcounts, MPI_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs,     MPI_SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    WORK_COUNT = sendcounts[MPI_RANK];

    if (WORK_COUNT > 0)
        WORK_BOUNDS = malloc(sizeof(int) * WORK_COUNT);

    /* Scatterv reparte porciones de tamano variable sin padding */
    MPI_Scatterv(
        all, sendcounts, displs, MPI_INT,
        WORK_BOUNDS, WORK_COUNT, MPI_INT,
        0, MPI_COMM_WORLD
    );

    if (MPI_RANK == 0) free(all);
    free(sendcounts);
    free(displs);
}

void division_trabajo(long int *local_count8, long int *local_count4, long int *local_count2)
{
    /* Fase 1: reina en esquina, BOUND1 = 2 .. SIZE-2 */
    preparar_trabajo(1, 2, SIZEE - 1);
    logica_hilos(local_count8, local_count4, local_count2);

    /* Fase 2: reina interior, BOUND1 = 1 .. (SIZE-2)/2 */
    preparar_trabajo(2, 1, (SIZE - 2) / 2);
    logica_hilos(local_count8, local_count4, local_count2);
}

double NQueens(void)
{
    long int LOCAL_COUNT8  = 0, LOCAL_COUNT4  = 0, LOCAL_COUNT2  = 0;
    long int GLOBAL_COUNT8 = 0, GLOBAL_COUNT4 = 0, GLOBAL_COUNT2 = 0;
    long int TOTAL, UNIQUE;

    double sumaLocal  = 0.0;
    double sumaGlobal = 0.0;
    int totalThreads  = MPI_SIZE * T;
    double tIni, tFin;

    SIZEE  = SIZE - 1;
    TOPBIT = 1 << SIZEE;
    MASK   = (1 << SIZE) - 1;

    tIni = dwalltime();

    division_trabajo(&LOCAL_COUNT8, &LOCAL_COUNT4, &LOCAL_COUNT2);

    /* Reduce las tres cuentas en una sola llamada usando arreglos */
    long int local_counts[3]  = { LOCAL_COUNT8,  LOCAL_COUNT4,  LOCAL_COUNT2  };
    long int global_counts[3] = { 0, 0, 0 };
    MPI_Reduce(local_counts, global_counts, 3, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    tFin = dwalltime();

    if (MPI_RANK == 0) {
        GLOBAL_COUNT8 = global_counts[0];
        GLOBAL_COUNT4 = global_counts[1];
        GLOBAL_COUNT2 = global_counts[2];

        UNIQUE = GLOBAL_COUNT8 + GLOBAL_COUNT4 + GLOBAL_COUNT2;
        TOTAL  = GLOBAL_COUNT8 * 8 + GLOBAL_COUNT4 * 4 + GLOBAL_COUNT2 * 2;

        printf("Numero de resultados: %ld - Unicas: %ld\n", TOTAL, UNIQUE);
    }

    /* Imprime el tiempo de cada hilo por proceso (en orden de rank) */
    for (int r = 0; r < MPI_SIZE; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (MPI_RANK == r) {
            for (int i = 0; i < T; i++) {
                sumaLocal += tiempo_hilos[i];
                printf("Rank %d, Thread %d: Tiempo = %f segundos\n", MPI_RANK, i, tiempo_hilos[i]);
            }
        }
    }

    MPI_Reduce(&sumaLocal, &sumaGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (MPI_RANK == 0)
        printf("Tiempo promedio total por hilo = %f segundos\n", sumaGlobal / totalThreads);

    return (tFin - tIni);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

    /*
     * CORRECCIÓN CLAVE: solo root valida y lee los parametros,
     * luego los difunde a todos los procesos con MPI_Bcast.
     * Sin esto, SIZE y T quedan en 0 en los ranks != 0,
     * provocando resultados incorrectos (basura o valores enormes).
     */
    int error = 1;

    if (MPI_RANK == 0) {
        if (argc < 2) {
            fprintf(stderr, "Uso: %s N [pthreads_por_proceso]\n", argv[0]);
            error = 0;
        } else {
            SIZE = atoi(argv[1]);
            T    = (argc > 2) ? atoi(argv[2]) : 1;

            if (SIZE < MINSIZE || SIZE > MAXSIZE) {
                fprintf(stderr, "N debe estar entre %d y %d\n", MINSIZE, MAXSIZE);
                error = 0;
            }
            if (T < 1)       T = 1;
            if (T > MAXSIZE) T = MAXSIZE;
        }
    }

    MPI_Bcast(&error, 1, MPI_INT,  0, MPI_COMM_WORLD);
    if (!error) { MPI_Finalize(); return 1; }

    /* Todos los procesos necesitan SIZE y T para funcionar correctamente */
    MPI_Bcast(&SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T,    1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    double tiempo = NQueens();

    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_RANK == 0) {
        printf("Tamaño del tablero:        %d\n", SIZE);
        printf("Procesos MPI:              %d\n", MPI_SIZE);
        printf("pthreads por proceso MPI:  %d\n", T);
        printf("Total de workers:          %d\n", MPI_SIZE * T);
        printf("Tiempo Total:              %f segundos\n", tiempo);
    }

    if (WORK_BOUNDS != NULL) {
        free(WORK_BOUNDS);
        WORK_BOUNDS = NULL;
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