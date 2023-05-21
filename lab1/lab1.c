#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void selectionSort(double arr[], int start, int end)
{
    int i, j, min_idx;
    // One by one move boundary of unsorted subarray
    for (i = start; i < end-1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i+1; j < end; j++)
          if (arr[j] < arr[min_idx])
            min_idx = j;
 
        // Swap the found minimum element with the first element
           if(min_idx != i)
            swap(&arr[min_idx], &arr[i]);
    }
}

void generate_array(double *restrict m, int size, unsigned int min, unsigned int max, int seed)
{
    //int tmp_seed;
    // #pragma omp for private(tmp_seed)
    for (int i = 0; i < size; ++i)
    {
        unsigned int tmp_seed = sqrt(i + seed);
        m[i] = ((double)rand_r(&tmp_seed) / (RAND_MAX)) * (max - min) + min;
    }
}

int main(int argc, char* argv[]) {
    int i, N;
    struct timeval T1, T2;
    long delta_ms;
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
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
    for (i=0; i < 100; ++i) { /* 100 экспериментов */
        /* инициализировать начальное значение ГСЧ */
        seed = i;
        /* Заполнить массив исходных данных размером N */
        // GENERATE
        generate_array(M1, N, min, max, seed);
        generate_array(M2, N2, max, max_2, seed+2);
        // MAP
        for (j = 0; j < N; ++j) {
            M1[j] = exp(sqrt(M1[j]));
        }
        for (k = 0; k < N2; ++k) {
            M2_old[k] = M2[k];
        }
        for(int k = 1; k < N2; ++k) {
            M2[k] = M2[k] + M2_old[k-1];
        }
        for(int k = 0; k < N2; ++k) {
            M2[k] = log(fabs(tan(M2[k])));
        }
        // MERGE
        for(k=0; k < N2; ++k) {
            M2[k]= M1[k] * M2[k];
        }
        // SORT
        selectionSort(M2, 0, N2);
        // REDUCE
        X = 0.0;
        key = M2[0];
        for (k = 1; k < N2; ++k) {
            if(M2[k] != 0) {
                if (key == 0 || M2[k] < key) {
                    key = M2[k];
                }
            }
        }
        for(int k = 0; k < N2; ++k) {
            if (((int)(M2[k] / key) % 2) == 0) {
                X += sin(M2[k]);
            }
        }
        //printf("X = %f", X);
        //printf("\n\n");
        /* Решить поставленную задачу, заполнить массив с результатами*/
        /* Отсортировать массив с результатами указанным методом */
    }
    printf("X= %f\n", X);
    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
    //printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */
    printf("%f\n", X);
    printf("%ld\n",delta_ms);
    return 0;
}
