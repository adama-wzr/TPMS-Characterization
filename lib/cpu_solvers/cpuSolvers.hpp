#ifndef _CPU_SOLVE
#define _CPU_SOLVE

#include <iostream>
#include <math.h>

#include <data_structures.hpp>

/*

    Gauss-Seidel Implementations:

*/

void pGS3D_inner(meshInfo *mesh, float *Coeff, float *RHS, float *x_vec)
{
    /*
        Function pGS3D_inner:
        Inputs:
            - pointer to mesh struct
            - pointer to Coeff matrix
            - pointer RHS
            - pointer to x_vec
        Outputs:
            - none
        
        Function will complete one iteration of the inner loop of
        the Guass-Seidel iteration.
    */

    // set variables, start parallel CPU code
    long int i = 0;
    float sigma = 0.0;

    #pragma omp parallel private (i, sigma)

    #pragma omp parallel for schedule(auto)
    for(i = 0; i < mesh->nElements; i++)
    {
        sigma = 0.0;

        if (Coeff[i*7 + 0] == 0)
            continue;
        
        for(int j = 1; j < 7; j++)
        {
            if (j == 1)
            {
                // West
                sigma += Coeff[i*7 + j] * x_vec[i - 1];
            }
            else if(j == 2)
            {
                // East
                sigma += Coeff[i*7 + j] * x_vec[i + 1];
            }
            else if(j == 3)
            {
                // South
                sigma += Coeff[i*7 + j] * x_vec[i + mesh->numCellsX];
            }
            else if(j == 4)
            {
                // North
                sigma += Coeff[i*7 + j] * x_vec[i - mesh->numCellsX];
            }
            else if(j == 5)
            {
                // Back
                sigma += Coeff[i*7 + j] * x_vec[i + mesh->numCellsX * mesh->numCellsY]; 
            }
            else if(j == 6)
            {
                // Front
                sigma += Coeff[i*7 + j] * x_vec[i - mesh->numCellsX * mesh->numCellsY];
            }
        }
        // update x
        x_vec[i] = 1/Coeff[i*7 + 0] * (RHS[i] - sigma);
    } // end for

    return;
}


int pGS3D_handle(options *opts, meshInfo *mesh, saveInfo *save, float *Coeff, float *RHS, float *x_vec)
{
    /*
        Function pGS3D_handle:
        Inputs:
            - pointer to options struct
            - pointer to mesh struct
            - pointer to save struct
            - pointer to Coeff matrix
            - pointer to RHS matrix
            - pointer to x_vec
        Outputs:
            - none
        
        Function is a parallel implementation of the Guass-Seidel iterative solver for a 3D
        problem where the coefficient matrix has 7 main diagonals. Solves Ax = b for x, where
        A is Coeff and b is RHS.
    */

    // declare useful variables
    long int nIter = 0;
    int iterToCheck = 1000;
    float conv = 1.0;
    long int i = 0;     // index is private
    float sigma = 0.0;  // update sigma is private

    // get array for comparison
    float *oldX = (float *)malloc(sizeof(float) * mesh->nElements);

    memcpy(oldX, x_vec, sizeof(float) * mesh->nElements);
    
    // main loop
    while(nIter < opts->MAX_ITER && conv > opts->ConvergeCriteria)
    {
        // call solver
        pGS3D_inner(mesh, Coeff, RHS, x_vec);
    
        // update 

        nIter++;

        if(nIter % iterToCheck == 0 && nIter != 0)
        {
            float sum = 0.0;
            long int count = 0;

            for(long int index = 0; index < mesh->nElements; index++)
            {
                if (x_vec[index] == 0)
                    continue;
                
                count++;
                sum += fabs((oldX[index] - x_vec[index])/oldX[index]);
            }
            
            conv = sum / count;
            
            memcpy(oldX, x_vec, sizeof(float) * mesh->nElements);
        }

        if(nIter % 10000 == 0 && opts->verbose == 1)
        {
            printf("Iteration = %ld, pct Change = %1.3e\n", nIter, conv);
        }
    }

    // Memory Management

    free(oldX);

    return 0;
}

void testHello(void)
{
    printf("\nHello!\n");
}

#endif