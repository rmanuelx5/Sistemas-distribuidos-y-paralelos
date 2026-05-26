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
int WORK_COUNT = 0;
int WORK_INDEX = 0;
int PHASE;
//exclusion mutua para tomar trabajo 
pthread_mutex_t mutex;

// contadores locales
long int THREAD_COUNT8[MAXSIZE];
long int THREAD_COUNT4[MAXSIZE];
long int THREAD_COUNT2[MAXSIZE];

//para que cada hilo guarde su tiempo de ejecucion
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

//mostrar tablero en pantalla
void Display(structHilo *s)
{
    int y, bit;

    printf("N= %d\n", SIZE);
    for (y=0; y<SIZE; y++) {
        for (bit=TOPBIT; bit; bit>>=1)
            printf("%s ", (s->BOARD[y] & bit)? "Q": "-");
        printf("\n");
    }
    printf("\n");
}

//verifica simetrias
void Check(structHilo *s)
{
    int *own, *you, bit, ptn;

    /* 90-degree rotation */
    if (*(s->BOARD2) == 1) {
        for (ptn=2, own=s->BOARD+1; own<=s->BOARDE; own++, ptn<<=1) {
            bit = 1;
            for (you=s->BOARDE; *you!=ptn && *own>=bit; you--)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
        if (own > s->BOARDE) {
            s->COUNT2++;
            return;
        }
    }

    /* 180-degree rotation */
    if (*(s->BOARDE) == s->ENDBIT) {
        for (you=s->BOARDE-1, own=s->BOARD+1; own<=s->BOARDE; own++, you--) {
            bit = 1;
            for (ptn=TOPBIT; ptn!=*you && *own>=bit; ptn>>=1)
                bit <<= 1;
            if (*own > bit) return;
            if (*own < bit) break;
        }
        if (own > s->BOARDE) {
            s->COUNT4++;
            return;
        }
    }

    /* 270-degree rotation */
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

/**********************************************/
/* First queen is inside                      */
/**********************************************/
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
            Backtrack2(y+1, (left | bit)<<1, down | bit, (right | bit)>>1, s);
        }
    }
}

/**********************************************/
/* First queen is in the corner               */
/**********************************************/
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
            Backtrack1(y+1, (left | bit)<<1, down | bit, (right | bit)>>1, s);
        }
    }
}

//funcion de cada pthread
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
            Backtrack1(2, (2 | bit)<<1, 1 | bit, bit>>1, &s);
        } else {
            s.BOUND2 = SIZE - 1 - s.BOUND1;
            s.BOARD1 = &s.BOARD[s.BOUND1];
            s.BOARD2 = &s.BOARD[s.BOUND2];
            s.SIDEMASK = TOPBIT | 1;
            s.LASTMASK = TOPBIT | 1;
            s.ENDBIT = TOPBIT >> 1;

            for (int i = 1; i < s.BOUND1; i++) {
                s.LASTMASK |= (s.LASTMASK >> 1) | (s.LASTMASK << 1);
                s.ENDBIT >>= 1;
            }

            s.BOARD[0] = bit = 1 << s.BOUND1;
            Backtrack2(1, bit<<1, bit, bit>>1, &s);
        }

        THREAD_COUNT8[id] += s.COUNT8;
        THREAD_COUNT4[id] += s.COUNT4;
        THREAD_COUNT2[id] += s.COUNT2;

        bound = tomar_trabajo();
    }
    
    tiempo_hilos[id] = dwalltime() - tiempo_hilos[id];

    return NULL;
}

void division_trabajo(int phase, int inicio, int fin)
{
    int total_trabajos = 0;
    int *todos_los_bounds = NULL;

    PHASE = phase;
    WORK_INDEX = 0;

    if (WORK_BOUNDS != NULL) {
        free(WORK_BOUNDS);
        WORK_BOUNDS = NULL;
    }

    if (MPI_RANK == 0) {
        if (fin >= inicio) {
            total_trabajos = fin - inicio + 1;
        }
        //variable para mandar a cada proceso la cantidad de trabajos que le toca hacer
        WORK_COUNT = (total_trabajos + MPI_SIZE - 1) / MPI_SIZE; //la division de trabajos entre procesos, redondeando hacia arriba  
    }

    //se manda a cada proceso la cantidad de trabajo 
    MPI_Bcast(&WORK_COUNT, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //se reserva espacio de manera dinamica segun la cantidad recibida 
    if (WORK_COUNT > 0) {
        WORK_BOUNDS = malloc(sizeof(int) * WORK_COUNT);

        //mensaje de error por si no se pudo reservar memoria
        if (WORK_BOUNDS == NULL) {
            fprintf(stderr, "Rank %d: error reservando WORK_BOUNDS\n", MPI_RANK);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    //solo root prepara el arreglo con los bounds y luego se reparte con scatter
    if (MPI_RANK == 0 && WORK_COUNT > 0) {
        int total_celdas = WORK_COUNT * MPI_SIZE; //total celdas = trabajo + padding si hubiese (marcado con -1)


        todos_los_bounds = malloc(sizeof(int) * total_celdas);
        if (todos_los_bounds == NULL) {
            fprintf(stderr, "Error reservando memoria para MPI_Scatter\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //padding con -1 para los procesos que no reciban trabajo 
        for (int i = 0; i < total_celdas; i++) {
            todos_los_bounds[i] = -1;
        }

        //llenar el arreglo con los bounds reales a repartir
        for (int i = 0; i < total_trabajos; i++) {
            todos_los_bounds[i] = inicio + i;
        }
    }

    //scatter envia a cada proceso su parte
    //MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
    if (WORK_COUNT > 0) {
        MPI_Scatter(todos_los_bounds, WORK_COUNT, MPI_INT, WORK_BOUNDS, WORK_COUNT, MPI_INT, 0, MPI_COMM_WORLD);
    }

    //root libera el arreglo temporal 
    if (MPI_RANK == 0) {
        free(todos_los_bounds);
    }
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
    for(int i=0;i<T;i++){
        threads_ids[i]=i;
        pthread_create(&threads[i],NULL, &worker,(void*)&threads_ids[i]);
    }

    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
        *local_count8 += THREAD_COUNT8[i];
        *local_count4 += THREAD_COUNT4[i];
        *local_count2 += THREAD_COUNT2[i];
    }

    pthread_mutex_destroy(&mutex);
}


/**********************************************/
/* Search of N-Queens                         */
/**********************************************/
void NQueens(void)
{
    long int LOCAL_COUNT8 = 0;
    long int LOCAL_COUNT4 = 0;
    long int LOCAL_COUNT2 = 0;

    long int GLOBAL_COUNT8 = 0;
    long int GLOBAL_COUNT4 = 0;
    long int GLOBAL_COUNT2 = 0;

    long int TOTAL, UNIQUE;

    SIZEE  = SIZE - 1;
    TOPBIT = 1 << SIZEE;
    MASK   = (1 << SIZE) - 1;

    //
    division_trabajo(1, 2, SIZEE - 1);
    logica_hilos(&LOCAL_COUNT8, &LOCAL_COUNT4, &LOCAL_COUNT2);

    /* Fase 2 del algoritmo original: BOUND1 = 1 .. (SIZE-2)/2 */
    division_trabajo(2, 1, (SIZE - 2) / 2);
    logica_hilos(&LOCAL_COUNT8, &LOCAL_COUNT4, &LOCAL_COUNT2);

    /* MPI_Reduce: junta los resultados parciales de todos los procesos. */
    MPI_Reduce(&LOCAL_COUNT8, &GLOBAL_COUNT8, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LOCAL_COUNT4, &GLOBAL_COUNT4, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LOCAL_COUNT2, &GLOBAL_COUNT2, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    double sumaLocal = 0.0;
    double sumaGlobal = 0.0;
    int totalThreads = MPI_SIZE * T;

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

    if (MPI_RANK == 0) {
        printf("Tiempo promedio total por hilo = %f segundos\n", sumaGlobal / totalThreads);
    }
    

    if (MPI_RANK == 0) {
        UNIQUE = GLOBAL_COUNT8 + GLOBAL_COUNT4 + GLOBAL_COUNT2;
        TOTAL  = GLOBAL_COUNT8 * 8 + GLOBAL_COUNT4 * 4 + GLOBAL_COUNT2 * 2;

        printf("Numero de resultados: %lu - Unicas: %lu\n", TOTAL, UNIQUE);
    }
}


int main(int argc, char *argv[])
{
    double tIni, tFin;
    int error = 1; //flag para indicar si los parametros son validos

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

    // El proceso root valida parametros y los distribuye con broadcast
    if (MPI_RANK == 0) {
        if (argc < 2) {
            fprintf(stderr, "Uso: %s N [pthreads_por_proceso]\n", argv[0]);
            error = 0;
        } else {
            SIZE = atoi(argv[1]);
            T = (argc > 2) ? atoi(argv[2]) : 1;

            if (SIZE < MINSIZE || SIZE > MAXSIZE) {
                fprintf(stderr, "N debe estar entre %d y %d\n", MINSIZE, MAXSIZE);
                error = 0;
            }

            if (T < 1) T = 1;
            if (T > MAXSIZE) T = MAXSIZE;
        }
    }
    
    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!error) {
        MPI_Finalize();
        return 1;
    }

    // Todos los procesos reciben los parametros ya validados 
    MPI_Bcast(&SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // Se mide solo el tiempo de computo
    tIni = dwalltime();

    NQueens();

    MPI_Barrier(MPI_COMM_WORLD);
    
    tFin = dwalltime();

    // Imprime configuracion, resultados y tiempos
    if (MPI_RANK == 0) {
        printf("Tamaño del tablero: %d\n", SIZE);
        printf("Procesos MPI: %d\n", MPI_SIZE);
        printf("pthreads por proceso MPI: %d\n", T);
        printf("Total de workers: %d\n", MPI_SIZE * T);
        printf("Tiempo Total: %f segundos\n", tFin - tIni);
    }

    // limpiar recursos
    // if (WORK_BOUNDS != NULL) {
    //     free(WORK_BOUNDS);
    //     WORK_BOUNDS = NULL;
    // }

    MPI_Finalize();
    return 0;
}

/* Queda por compatibilidad con el archivo original. MPI_Wtime se usa en main. */
double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}
