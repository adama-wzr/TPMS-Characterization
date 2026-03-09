#ifndef _CPU_SOLVE
#define _CPU_SOLVE

#include <iostream>
#include <math.h>
#include <omp.h>
#include <string.h>

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

    #pragma omp parallel for schedule(auto)
    for(int i = 0; i < mesh->nElements; i++)
    {
        float sigma = 0.0;

        for(int j = 1; j < 7; j++)
        {
            if (Coeff[i*7 + j] == 0)
                continue;
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

void pGS3D_innerPB(meshInfo *mesh, float *Coeff, float *RHS, float *x_vec)
{
    /*
        Function pGS3D_innerPB:
        Inputs:
            - pointer to mesh struct
            - pointer to Coeff matrix
            - pointer RHS
            - pointer to x_vec
        Outputs:
            - none
        
        Function will complete one iteration of the inner loop of
        the Guass-Seidel iteration, with periodic boundaries.
        (periodic in y and z only)
    */

    // set variables, start parallel CPU code

    int nCols = mesh->numCellsX;
    int nRows = mesh->numCellsY;
    int nSlices = mesh->numCellsZ;

    #pragma omp parallel for schedule(auto)
    for(int i = 0; i < mesh->nElements; i++)
    {
        float sigma = 0.0;
        
        int mySlice = i/(nRows * nCols);
        int myRow = (i - mySlice * nRows * nCols)/nCols;
        int myCol = i - mySlice * nRows * nCols - myRow * nCols;
        
        for(int j = 1; j < 7; j++)
        {
            if (Coeff[i*7 + j] == 0)
                continue;
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
                if(myRow == nRows - 1)
                {
                    // periodic south
                    sigma += Coeff[i*7 + j] * x_vec[mySlice * nRows * nCols + myCol];
                }
                else
                {
                    // South
                    sigma += Coeff[i*7 + j] * x_vec[i + nCols];
                }
            }
            else if(j == 4)
            {
                if(myRow == 0)
                {
                    // Periodic North
                    sigma += Coeff[i*7 + j] * x_vec[mySlice * nRows * nCols + (nRows - 1) * nCols + myCol];    
                }
                else
                {
                    // North
                    sigma += Coeff[i*7 + j] * x_vec[i - nCols];
                }
            }
            else if(j == 5)
            {
                if(mySlice == nSlices - 1)
                {
                    // periodic Back
                    sigma += Coeff[i*7 + j] * x_vec[myRow * nCols + myCol];
                }
                else
                {
                    // Back
                    sigma += Coeff[i*7 + j] * x_vec[i + nCols * nRows];
                } 
            }
            else if(j == 6)
            {
                if(mySlice == 0)
                {
                    // Periodic Front
                    sigma += Coeff[i*7 + j] * x_vec[(nSlices - 1) * nCols * nRows + myRow * nCols + myCol];    
                }
                else
                {
                    // Front
                    sigma += Coeff[i*7 + j] * x_vec[i - nCols * nRows];
                }
            }
        }
        // update x
        x_vec[i] = 1.0/Coeff[i*7 + 0] * (RHS[i] - sigma);
    } // end for

    return;
}

void pGS3D_allPB(meshInfo *mesh, float *Coeff, float *RHS, float *x_vec)
{
    /*
        Function pGS3D_allPB:
        Inputs:
            - pointer to mesh struct
            - pointer to Coeff matrix
            - pointer RHS
            - pointer to x_vec
        Outputs:
            - none
        
        Function will complete one iteration of the inner loop of
        the Guass-Seidel iteration, with periodic boundaries.
        (periodic in x, y, and z)
    */

    // set variables, start parallel CPU code

    int nCols = mesh->numCellsX;
    int nRows = mesh->numCellsY;
    int nSlices = mesh->numCellsZ;

    #pragma omp parallel for schedule(auto)
    for(int i = 0; i < mesh->nElements; i++)
    {
        float sigma = 0.0;
        
        int mySlice = i/(nRows * nCols);
        int myRow = (i - mySlice * nRows * nCols)/nCols;
        int myCol = i - mySlice * nRows * nCols - myRow * nCols;
        
        for(int j = 1; j < 7; j++)
        {
            if (Coeff[i*7 + j] == 0)
                continue;
            if (j == 1)
            {
                if(myCol == 0)
                {
                    // Periodic West
                    sigma += Coeff[i*7 + j] * x_vec[i + nCols - 1];
                }else
                {
                    // West
                    sigma += Coeff[i*7 + j] * x_vec[i - 1]; 
                }
            }
            else if(j == 2)
            {
                if(myCol == nCols - 1)
                {
                    // Periodic East
                    sigma += Coeff[i*7 + j] * x_vec[i - (nCols + 1)];
                }
                else
                {
                    // East
                    sigma += Coeff[i*7 + j] * x_vec[i + 1];
                }
            }
            else if(j == 3)
            {
                if(myRow == nRows - 1)
                {
                    // periodic south
                    sigma += Coeff[i*7 + j] * x_vec[mySlice * nRows * nCols + myCol];
                }
                else
                {
                    // South
                    sigma += Coeff[i*7 + j] * x_vec[i + nCols];
                }
            }
            else if(j == 4)
            {
                if(myRow == 0)
                {
                    // Periodic North
                    sigma += Coeff[i*7 + j] * x_vec[mySlice * nRows * nCols + (nRows - 1) * nCols + myCol];    
                }
                else
                {
                    // North
                    sigma += Coeff[i*7 + j] * x_vec[i - nCols];
                }
            }
            else if(j == 5)
            {
                if(mySlice == nSlices - 1)
                {
                    // periodic Back
                    sigma += Coeff[i*7 + j] * x_vec[myRow * nCols + myCol];
                }
                else
                {
                    // Back
                    sigma += Coeff[i*7 + j] * x_vec[i + nCols * nRows];
                } 
            }
            else if(j == 6)
            {
                if(mySlice == 0)
                {
                    // Periodic Front
                    sigma += Coeff[i*7 + j] * x_vec[(nSlices - 1) * nCols * nRows + myRow * nCols + myCol];    
                }
                else
                {
                    // Front
                    sigma += Coeff[i*7 + j] * x_vec[i - nCols * nRows];
                }
            }
        }
        // update x
        x_vec[i] = 1.0/Coeff[i*7 + 0] * (RHS[i] - sigma);
    } // end for

    return;
}

/*

    Handles for Solvers:

*/

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

    // set num threads
    omp_set_num_threads(opts->nThreads);

    // declare useful variables
    long int nIter = 0;
    int iterToCheck = 1000;
    float conv = 1.0;

    // get array for comparison
    float *oldX = (float *)malloc(mesh->nElements * sizeof(float));

    memcpy(oldX, x_vec, sizeof(float) * mesh->nElements);
    
    // main loop
    while(nIter < opts->MAX_ITER && conv > opts->ConvergeCriteria)
    {
        // call solver
        if(opts->PB)
        {
            // periodic boundary
            pGS3D_innerPB(mesh, Coeff, RHS, x_vec);
        }
        else
        {
            // no flux boundaries
            pGS3D_inner(mesh, Coeff, RHS, x_vec);
        }
    
        // update iter count
        nIter++;
        // check convergence (if applicable)

        if(nIter % iterToCheck == 0 && nIter != 0)
        {
            float sum = 0.0;
            long int count = 0;

            for(long int index = 0; index < mesh->nElements; index++)
            {
                if (x_vec[index] == 0)
                    continue;
                
                count++;
                sum += fabs((x_vec[index] - oldX[index])/x_vec[index]);
            }
            
            conv = sum / count;
            
            memcpy(oldX, x_vec, sizeof(float) * mesh->nElements);
        }

        // update user (if applicable)
        if(nIter % 1000 == 0 && opts->verbose == 1)
        {
            printf("Iteration = %ld, pct Change = %1.3e\n", nIter, conv);
        }
    }

    // Memory Management

    free(oldX);

    return 0;
}

int pGS3D_SF_handle(options *opts, meshInfo *mesh, saveInfo *save, float *Coeff, float *RHS, float *x_vec)
{
    /*
        Function pGS3D_SF_handle:
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
        A is Coeff and b is RHS. This one is specifically for the Shape Factor calculation
    */

    // set num threads
    omp_set_num_threads(opts->nThreads);

    // declare useful variables
    long int nIter = 0;
    int iterToCheck = 100;
    float conv = 1.0;

    // get array for comparison
    float *oldX = (float *)malloc(mesh->nElements * sizeof(float));

    memcpy(oldX, x_vec, sizeof(float) * mesh->nElements);
    
    // main loop
    while(nIter < opts->MAX_ITER && conv > opts->ConvergeCriteria)
    {
        pGS3D_allPB(mesh, Coeff, RHS, x_vec);
    
        // update iter count
        nIter++;
        // check convergence (if applicable)

        if(nIter % iterToCheck == 0 && nIter != 0)
        {
            float sum = 0.0;
            long int count = 0;

            for(long int index = 0; index < mesh->nElements; index++)
            {
                if (x_vec[index] == 0)
                    continue;
                
                count++;
                sum += fabs((x_vec[index] - oldX[index])/x_vec[index]);
            }
            
            conv = sum / count;
            
            memcpy(oldX, x_vec, sizeof(float) * mesh->nElements);
        }

        // update user (if applicable)
        if(nIter % iterToCheck == 0 && opts->verbose == 1)
        {
            printf("Iteration = %ld, pct Change = %1.3e\n", nIter, conv);
        }
    }

    // Memory Management

    free(oldX);

    return 0;
}

#endif