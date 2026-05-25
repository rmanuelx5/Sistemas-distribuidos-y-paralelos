/**************************************************************************/
/* N-Queens Solutions ver3.1 - MPI + pthreads                               */
/* Base: Takaken July/2003                                                  */
/* MPI reparte bounds entre procesos; pthreads reparte dentro de cada nodo.  */
/* Regla conservada: a cada pthread solo se le pasa su id como argumento.    */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <sys/time.h>
#include <mpi.h>

#define MAXSIZE 24
#define MINSIZE 2


//holis

double dwalltime();

int SIZE, SIZEE;
int MASK, TOPBIT;

//cant hilos por proceso MPI
int T;

// MPI
int WORLD_RANK;
int WORLD_SIZE;


//compartidas
int NEXT_BOUND;
int END_BOUND;
int PHASE;
pthread_mutex_t mutex;

// array de contadores 
long int THREAD_COUNT8[MAXSIZE];
long int THREAD_COUNT4[MAXSIZE];
long int THREAD_COUNT2[MAXSIZE];

// variables de cada hilo, no compartido
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

// Tomar siguiente trabajo                    
int tomar_trabajo(void)
{
    int bound;

    //bloquea mutex
    pthread_mutex_lock(&mutex);
    if (NEXT_BOUND > END_BOUND) {
        bound = -1;
    } else {
        bound = NEXT_BOUND;

        /* MPI: cada proceso toma solamente sus bounds:
           rank, rank+WORLD_SIZE, rank+2*WORLD_SIZE, ...
           Ejemplo fase con inicio 2:
           rank 0 -> 2, 2+P, ...
           rank 1 -> 3, 3+P, ...
        */
        NEXT_BOUND += WORLD_SIZE;
    }
    //larga mutex
    pthread_mutex_unlock(&mutex);

    return bound;
}

/**********************************************/
/* Display the Board Image                    */
/**********************************************/
void Display(structHilo *s)
{
    int  y, bit;

    printf("N= %d\n", SIZE);
    for (y=0; y<SIZE; y++) {
        for (bit=TOPBIT; bit; bit>>=1)
            printf("%s ", (s->BOARD[y] & bit)? "Q": "-");
        printf("\n");
    }
    printf("\n");
}

/**********************************************/
/* Check Unique Solutions                     */
/**********************************************/
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

/**********************************************/
/* Worker: recibe solamente su id             */
/**********************************************/
void *worker(void *arg)
{
    int id = (int)(intptr_t)arg;
    int bound = tomar_trabajo();
    
    while (bound != -1) {
        structHilo s = {0};
        int bit;

        s.BOARDE = &s.BOARD[SIZEE];
        s.BOUND1 = bound;

        if (PHASE == 1) {
            /* Equivale al primer for del NQueens original. */
            s.BOARD[0] = 1;
            s.BOARD[1] = bit = 1 << s.BOUND1;
            Backtrack1(2, (2 | bit)<<1, 1 | bit, bit>>1, &s);
        } else {
            /* Equivale al segundo for del NQueens original. */
            s.BOUND2 = SIZE - 1 - s.BOUND1;
            s.BOARD1 = &s.BOARD[s.BOUND1];
            s.BOARD2 = &s.BOARD[s.BOUND2];
            s.SIDEMASK = TOPBIT | 1;
            s.LASTMASK = TOPBIT | 1;
            s.ENDBIT = TOPBIT >> 1;

            /* Reproduce el estado que tendrian LASTMASK y ENDBIT
               al llegar a este BOUND1 en el for secuencial original. */
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

    return NULL;
}

/**********************************************/
/* Ejecutar una fase                          */
/**********************************************/
void ejecutar_fase(long int *total_count8, long int *total_count4, long int *total_count2)
{
    pthread_t threads[MAXSIZE];

    for (int i = 0; i < T; i++) {
        THREAD_COUNT8[i] = 0;
        THREAD_COUNT4[i] = 0;
        THREAD_COUNT2[i] = 0;
    }

    pthread_mutex_init(&mutex, NULL);

    for (int i = 0; i < T; i++) {
        pthread_create(&threads[i], NULL, worker, (void *)(intptr_t)i);
    }

    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
        *total_count8 += THREAD_COUNT8[i];
        *total_count4 += THREAD_COUNT4[i];
        *total_count2 += THREAD_COUNT2[i];
    }

    pthread_mutex_destroy(&mutex);
}

/**********************************************/
/* Search of N-Queens                         */
/**********************************************/
void NQueens(void)
{
    long int COUNT8 = 0;
    long int COUNT4 = 0;
    long int COUNT2 = 0;
    long int TOTAL, UNIQUE;

    SIZEE  = SIZE - 1;
    TOPBIT = 1 << SIZEE;
    MASK   = (1 << SIZE) - 1;

    /* 0:000000001 */
    /* 1:011111100 */

    
    PHASE = 1;
    NEXT_BOUND = 2;
    END_BOUND = SIZEE - 1;

    ejecutar_fase(&COUNT8, &COUNT4, &COUNT2);


    PHASE = 2;
    NEXT_BOUND = 1;
    END_BOUND = (SIZE - 2) / 2;

    /* 0:000001110 */
    ejecutar_fase(&COUNT8, &COUNT4, &COUNT2);

    UNIQUE = COUNT8     + COUNT4     + COUNT2;
    TOTAL  = COUNT8 * 8 + COUNT4 * 4 + COUNT2 * 2;

    printf("Numero de resultados: %lu - Unicas: %lu\n", TOTAL, UNIQUE);
}

/**********************************************/
/* N-Queens Solutions MAIN                    */
/**********************************************/
int main(int argc, char *argv[])
{
    double tIni, tFin;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &WORLD_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &WORLD_SIZE);

    if (argc < 2) {
        if (WORLD_RANK == 0) {
            fprintf(stderr, "Uso: %s N [threads_por_proceso]\n", argv[0]);
            fprintf(stderr, "Ejemplo: mpirun -np 4 %s 14 2\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    SIZE = atoi(argv[1]);
    if (SIZE < MINSIZE || SIZE > MAXSIZE) {
        if (WORLD_RANK == 0) {
            fprintf(stderr, "N debe estar entre %d y %d\n", MINSIZE, MAXSIZE);
        }
        MPI_Finalize();
        return 1;
    }

    T = (argc > 2) ? atoi(argv[2]) : 1;
    if (T > MAXSIZE) T = MAXSIZE;

    MPI_Barrier(MPI_COMM_WORLD);
    tIni = dwalltime();

    NQueens();

    MPI_Barrier(MPI_COMM_WORLD);
    tFin = dwalltime();

    if (WORLD_RANK == 0) {
        printf("Procesos MPI: %d - pthreads por proceso: %d\n", WORLD_SIZE, T);
        printf("Tiempo Total: %f segundos\n", tFin - tIni);
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