//
// Created by Fernando Cores Prado on 4/12/23.
//

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <memory.h>
#include "Standard_MultMat.h"
#include "Matrix.h"
#include "Errors.h"
#include <pthread.h>

double elapsed_std;

struct MultiplyThread{
    int start;
    int end;
    float **matrixA;
    float **matrixB;
    float **result;
    int n;
};
typedef struct MultiplyThread MThread;

void* standardThread(void* arg) {
    MThread* m_thread = (MThread*)arg;
    for (int i = m_thread->start; i < m_thread->end; i++) {
        for (int j = 0; j < m_thread->n; j++) {
            float sum = 0;
            for (int k = 0; k < m_thread->n; k++) {
                sum += m_thread->matrixA[i][k] * m_thread->matrixB[k][j];
            }
            m_thread->result[i][j] = sum;
        }
    }
    pthread_exit(NULL);
}

float ** concurrentStandardMultiplication(float **matrixA, float **matrixB, int n, int num_threads) {
    pthread_t threads[num_threads];
    MThread mThread[num_threads];

    int rows_per_thread = n / num_threads;
    int extra_rows = n % num_threads;
    int start= 0;

    float** result = allocateMatrix(n); 

    struct timespec start_time, finish_time;

    clock_gettime(CLOCK_MONOTONIC, &start_time);

    for (int i = 0; i < num_threads; i++) {
        mThread[i].start= start;
        mThread[i].end= start + rows_per_thread + (i < extra_rows? 1 : 0);
        start= mThread[i].end;

        mThread[i].matrixA = matrixA;
        mThread[i].matrixB = matrixB;
        mThread[i].result = result;
        mThread[i].n = n;

        if (pthread_create(&threads[i], NULL, standardThread, (void*) &mThread[i]) != 0) {
            perror("Failed to create thread");
            for (int j = 0; j < i; j++) {
                pthread_cancel(threads[j]);
            }
            freeMatrix(result, n);
            exit(EXIT_FAILURE);
        }

        if (extra_rows > 0) extra_rows--;
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_time);
    elapsed_std = (finish_time.tv_sec - start_time.tv_sec);
    elapsed_std += (finish_time.tv_nsec - start_time.tv_nsec) / 1000000000.0;

    return result;

}



/*
* Standard Matrix multiplication with O(n^3) time complexity.
*/
float ** standardMultiplication(float ** matrixA,float ** matrixB,int n)
{
    return standardMultiplication_ijk(matrixA,matrixB,n);
    //return standardMultiplication_ikj(matrixA,matrixB,n);
}

/*
* Standard ijk Matrix multiplication with O(n^3) time complexity.
*/
float ** standardMultiplication_ijk(float ** matrixA,float ** matrixB,int n)
{
    struct timespec start, finish;
    float ** result;
    int i,j,k;

    clock_gettime(CLOCK_MONOTONIC, &start);

    result = (float**)malloc(n*sizeof(float *));
    for(i=0;i<n;i++){  
        result[i]=(float*)malloc(n*sizeof(float));
        memset(result[i],0,n*sizeof(float));
        for(j=0;j<n;j++){
            for(k=0;k<n;k++) {
                result[i][j]=result[i][j]+(matrixA[i][k]*matrixB[k][j]);
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed_std = (finish.tv_sec - start.tv_sec);
    elapsed_std += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    return result;
}

/*
* Standard ikj Matrix multiplication with O(n^3) time complexity.
*/
float ** standardMultiplication_ikj(float ** matrixA,float ** matrixB,int n)
{
    struct timespec start, finish;
    float ** result;
    int i,j,k;

    clock_gettime(CLOCK_MONOTONIC, &start);

    result = (float**)malloc(n*sizeof(float *));
    for(i=0;i<n;i++){
        result[i]=(float*)malloc(n*sizeof(float));
        memset(result[i],0,n*sizeof(float));
        for(k=0;k<n;k++) {
            for(j=0;j<n;j++){
                result[i][j]=result[i][j]+(matrixA[i][k]*matrixB[k][j]);
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed_std = (finish.tv_sec - start.tv_sec);
    elapsed_std += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    return result;
}
