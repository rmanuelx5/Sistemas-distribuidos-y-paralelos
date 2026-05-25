/**************************************************************************/
/* N-Queens Solutions ver3.1 - pthreads                                    */
/* Base: Takaken July/2003                                                  */
/* Regla pedida: a cada pthread solo se le pasa su id como argumento.        */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <sys/time.h>

#define MAXSIZE 24
#define MINSIZE 2

/* Time in seconds from some point in the past */
double dwalltime();

int SIZE, SIZEE;
int MASK, TOPBIT;

/* Cantidad de hilos */
int NUM_THREADS;

/* Trabajo global compartido: los hilos no reciben esto por argumento. */
int NEXT_BOUND;
int END_BOUND;
int PHASE;
pthread_mutex_t WORK_MUTEX;

/* Contadores globales por hilo. Cada hilo usa su posicion segun su id. */
long int THREAD_COUNT8[MAXSIZE];
long int THREAD_COUNT4[MAXSIZE];
long int THREAD_COUNT2[MAXSIZE];

/* Estado local de busqueda de cada hilo. */
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
} ThreadState;

/**********************************************/
/* Tomar siguiente trabajo                    */
/**********************************************/
int tomar_trabajo(void)
{
    int bound;

    pthread_mutex_lock(&WORK_MUTEX);
    if (NEXT_BOUND > END_BOUND) {
        bound = -1;
    } else {
        bound = NEXT_BOUND;
        NEXT_BOUND++;
    }
    pthread_mutex_unlock(&WORK_MUTEX);

    return bound;
}

/**********************************************/
/* Check Unique Solutions                     */
/**********************************************/
void Check(ThreadState *s)
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
void Backtrack2(int y, int left, int down, int right, ThreadState *s)
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
void Backtrack1(int y, int left, int down, int right, ThreadState *s)
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
        ThreadState s = {0};
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
void ejecutar_fase(int phase, int start_bound, int end_bound,
                   long int *total_count8, long int *total_count4, long int *total_count2)
{
    pthread_t threads[MAXSIZE];

    /* CAMBIADO: en pthread_create se pasa solamente el id del hilo, no args. */
    PHASE = phase;
    NEXT_BOUND = start_bound;
    END_BOUND = end_bound;

    for (int i = 0; i < NUM_THREADS; i++) {
        THREAD_COUNT8[i] = 0;
        THREAD_COUNT4[i] = 0;
        THREAD_COUNT2[i] = 0;
    }

    pthread_mutex_init(&WORK_MUTEX, NULL);

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_create(&threads[i], NULL, worker, (void *)(intptr_t)i);
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        *total_count8 += THREAD_COUNT8[i];
        *total_count4 += THREAD_COUNT4[i];
        *total_count2 += THREAD_COUNT2[i];
    }

    pthread_mutex_destroy(&WORK_MUTEX);
}

/**********************************************/
/* Search of N-Queens                         */
/**********************************************/
void NQueens(void)
{
    long int TOTAL, UNIQUE;
    
    long int COUNT8 = 0;
    long int COUNT4 = 0;
    long int COUNT2 = 0;
    SIZEE  = SIZE - 1;
    
    TOPBIT = 1 << SIZEE;
    MASK   = (1 << SIZE) - 1;

    /* 0:000000001 */
    /* 1:011111100 */
    ejecutar_fase(1, 2, SIZEE - 1, &COUNT8, &COUNT4, &COUNT2);

    /* 0:000001110 */
    ejecutar_fase(2, 1, (SIZE - 2) / 2, &COUNT8, &COUNT4, &COUNT2);

    UNIQUE = COUNT8     + COUNT4     + COUNT2;
    TOTAL  = COUNT8 * 8 + COUNT4 * 4 + COUNT2 * 2;
/**************************************************************************/
/* N-Queens Solutions ver3.1 - pthreads                                    */
/* Base: Takaken July/2003                                                  */
/* Regla pedida: a cada pthread solo se le pasa su id como argumento.        */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <sys/time.h>

#define MAXSIZE 24
#define MINSIZE 2


double dwalltime();

int SIZE, SIZEE;
int MASK, TOPBIT;

//cant hilos
int T;


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
        NEXT_BOUND++;
    }
    //larga mutex
    pthread_mutex_unlock(&mutex);

    return bound;
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
void ejecutar_fase(int phase, int start_bound, int end_bound,
                   long int *total_count8, long int *total_count4, long int *total_count2)
{
    pthread_t threads[MAXSIZE];

    PHASE = phase;
    NEXT_BOUND = start_bound;
    END_BOUND = end_bound;

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
    ejecutar_fase(1, 2, SIZEE - 1, &COUNT8, &COUNT4, &COUNT2);

    /* 0:000001110 */
    ejecutar_fase(2, 1, (SIZE - 2) / 2, &COUNT8, &COUNT4, &COUNT2);

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

    if (argc < 2) {
        fprintf(stderr, "Uso: %s N [threads]\n", argv[0]);
        return 1;
    }

    SIZE = atoi(argv[1]);
    if (SIZE < MINSIZE || SIZE > MAXSIZE) {
        fprintf(stderr, "N debe estar entre %d y %d\n", MINSIZE, MAXSIZE);
        return 1;
    }

    T = (argc > 2) ? atoi(argv[2]) : 16;
    if (T < 1) T = 1;
    if (T > MAXSIZE) T = MAXSIZE;

    tIni = dwalltime();
    NQueens();
    tFin = dwalltime();

    printf("Tiempo Total: %f segundos\n", tFin - tIni);
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
    printf("Numero de resultados: %lu - Unicas: %lu\n", TOTAL, UNIQUE);
}

/**********************************************/
/* N-Queens Solutions MAIN                    */
/**********************************************/
int main(int argc, char *argv[])
{
    double tIni, tFin;

    if (argc < 2) {
        fprintf(stderr, "Uso: %s N [threads]\n", argv[0]);
        return 1;
    }

    SIZE = atoi(argv[1]);
    
    if (SIZE < MINSIZE || SIZE > MAXSIZE) {
        fprintf(stderr, "N debe estar entre %d y %d\n", MINSIZE, MAXSIZE);
        return 1;
    }

    NUM_THREADS = (argc > 2) ? atoi(argv[2]) : 16;

    if (NUM_THREADS < 1) NUM_THREADS = 1;
    if (NUM_THREADS > MAXSIZE) NUM_THREADS = MAXSIZE;

    tIni = dwalltime();
    NQueens();
    tFin = dwalltime();

    printf("Tiempo Total: %f segundos\n", tFin - tIni);
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