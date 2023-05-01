#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include <omp.h>


#define _OPENMP
#define CHUNK 1
#define SCHEDULE guided


int main(int argc, char* argv[]) {
    int i, N;
    struct timeval T1, T2;
    long delta_ms;
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    #if defined(_OPENMP)
        omp_set_dynamic(0);
        const int THREADS = atoi(argv[2]);
        omp_set_num_threads(THREADS);
    #endif
    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
    int N2 = N/2; /* N2 равен N/2*/
    double* restrict M1 = (double*) malloc(N * sizeof(double)); 
    double* restrict M2 = (double*) malloc(N2 * sizeof(double)); /* Массивы M1 разм N и M2 разм N2*/
    double* restrict M2_old = (double*) malloc(N2 * sizeof(double));
    double A = 490.0; /* Ф*И*О */
    double min = 1; double max = A; double max_2 = max * 10;
    double key, X;
    int j, k;
    //int z, min_s;
    unsigned int seed;
    unsigned int* restrict seedp = &seed;
    unsigned int* restrict seedp1 = &seed;
    for (i=0; i < 100; ++i) { /* 100 экспериментов */
        /* инициализировать начальное значение ГСЧ */
        seed = i;
        X = 0.0;
        /* Заполнить массив исходных данных размером N */
        // GENERATE
        for(j=0; j < N; ++j) {
            M1[j]=((double)rand_r(seedp) / (RAND_MAX)) * (max - min) + min;
        }
        for (k=0; k < N2; ++k) {
            M2[k]=((double)rand_r(seedp1) / (RAND_MAX)) * (max_2 - max) + max;
        }
        /*----------------------------------------------------------------------*/
        #pragma omp parallel default(none) shared(N, N2, M1, M2, M2_old, key, X) private(k)
        {
            // MAP
            #if defined(CHUNK) && defined(SCHEDULE)
            #pragma omp for schedule(SCHEDULE, CHUNK) nowait
            #else
            #pragma omp for nowait
            #endif
            for (k = 0; k < N; ++k) {
                M1[k] = exp(sqrt(M1[k]));
            }
            #if defined(CHUNK) && defined(SCHEDULE)
            #pragma omp for schedule(SCHEDULE, CHUNK) nowait
            #else
            #pragma omp for nowait
            #endif
            for (k = 0; k < N2; ++k) {
                M2_old[k] = M2[k];
            }
            #if defined(CHUNK) && defined(SCHEDULE)
            #pragma omp for schedule(SCHEDULE, CHUNK)
            #else
            #pragma omp for
            #endif
            for(k = 1; k < N2; ++k) {
                M2[k] = M2[k] + M2_old[k-1];
            }
            #if defined(CHUNK) && defined(SCHEDULE)
            #pragma omp for schedule(SCHEDULE, CHUNK)
            #else
            #pragma omp for
            #endif
            for(k = 0; k < N2; ++k) {
                M2[k] = log(fabs(tan(M2[k])));
            }
            // MERGE
            #if defined(CHUNK) && defined(SCHEDULE)
            #pragma omp for schedule(SCHEDULE, CHUNK)
            #else
            #pragma omp for
            #endif
            for(k=0; k < N2; ++k) {
                M2[k]= M1[k] * M2[k];
            }
            /*----------------------------------------------------------------------*/
            // SORT
            /*
            for (k = 0; k < (N2-1); ++k) {
                min_s = k;
                for (j = k + 1; j < N2; ++j) {
                    if (M2[min_s] > M2[j]) {
                        min_s = j;
                    }
                }
                key = M2[min_s];
                for(z = min_s; z > k; --z) {
                    M2[z] = M2[z-1];
                }
                M2[k] = key;
            }
            */
        /*----------------------------------------------------------------------*/
            // REDUCE            
            #pragma omp single
            key = M2[0];
            #pragma omp single
            for (k = 1; k < N2; ++k) {
                if(M2[k] != 0) {
                    if (key == 0 || M2[k] < key) {
                        key = M2[k];
                    }
                }
            }
            #if defined(CHUNK) && defined(SCHEDULE)
            #pragma omp for reduction(+ : X) schedule(SCHEDULE, CHUNK)
            #else
            #pragma omp for reduction(+ : X)
            #endif
            for(k = 0; k < N2; k++) {
                if (((int)(M2[k] / key) % 2) == 0) {
                    X += sin(M2[k]);
                }
            }
            /*----------------------------------------------------------------------*/
        }
        //printf("X = %f", X);
        //printf("\n\n");
        /* Решить поставленную задачу, заполнить массив с результатами*/
        /* Отсортировать массив с результатами указанным методом */
    }
    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
    //printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */
    printf("%ld\n",delta_ms);
    free(M1);
    free(M2);
    free(M2_old);
    return 0;
}
