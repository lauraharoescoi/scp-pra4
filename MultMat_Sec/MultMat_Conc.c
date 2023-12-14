//
// Created by Fernando Cores Prado on 4/12/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Strassens_MultMat.h"
#include "Standard_MultMat.h"
#include "Matrix.h"
#include "Errors.h"

// Constants
#define DEBUG 0
char *usage_msg = "Usage: MultMat_Conc <MatrixA_File> <MatrixB_File> <ResultMatrix_File> <#Threads>\n";
char *input_folder = "Input";
char *results_folder = "Results";

// Global variables
char debug_msg [256];
struct timespec start, finish;
double elapsed;

/*
* Main function where the execution starts.
*/
int main(int argc, char ** argv)
{
    int n=0;
    char random_matrixA_name[256], random_matrixB_name[256], std_result_matrix_name[256], str_result_matrix_name[256];
    char * matrixA_name, * matrixB_name;
    float **matrixA, **matrixB;
    float **standardRes, **strassenRes;
    int num_threads;

    if (DEBUG) {
        sprintf(debug_msg,"[Main] Start Program with %d parameters.\n",argc);
        printMessage(debug_msg, COLOR_MAGENTA);
    }
    if(argc==5)
    {
        int n1, n2;
        double int_part;
        
        openMatrix(argv[1], &matrixA,&n1);
        openMatrix(argv[2], &matrixB,&n2);
        if (n1!=n2)
            Error("[Main]: Error input matrices have differente dimensions!");
        n=n1;

        int pow = 1;
        while (pow < n) {
            pow = pow * 2;
        }
        if (n!=pow)
            Error("[Main]: Error input matrices are not power of k!");

        matrixA_name = argv[1];
        matrixA_name = argv[2];
        sprintf(std_result_matrix_name, "%s.std", argv[3]);
        sprintf(str_result_matrix_name, "%s.str", argv[3]);
        num_threads = atoi(argv[4]);

    }
    else {
        printMessage(usage_msg,COLOR_RED);
        Error("[Main]: Invalid number of arguments!\n\n");
    }

    if (n<10) {
        printMessage("Matrix A:\n", COLOR_GREEN_B);
        printMatrixC(matrixA, n, COLOR_GREEN_B);

        printMessage("Matrix B:\n", COLOR_GREEN_B);
        printMatrixC(matrixB, n, COLOR_GREEN_B);
    }

    float ** stdRes = concurrentStandardMultiplication(matrixA,matrixB,n, num_threads);
    if (n<10) {
        print("Standard Multiplication Result:\n");
        printMatrix(stdRes, n);
    }
    saveMatrix(std_result_matrix_name, stdRes, n);
    sprintf(debug_msg,"[Standard Mult] Multiplication time: %05.6f.\n",elapsed_std);
    printMessage(debug_msg,COLOR_CYAN_B);


    float ** strassensRes = concurrentStrassensMultiplication(matrixA,matrixB,n, num_threads);
    if (n<10) {
        print("Strassen's Multiplication Result:\n");
        printMatrix(strassensRes, n);
    }
    saveMatrix(str_result_matrix_name, strassensRes, n);
    sprintf(debug_msg,"[Strassen Mult] Multiplication time: %05.6f.\n",elapsed_str);
    printMessage(debug_msg,COLOR_CYAN_B);

    return 0;
}

