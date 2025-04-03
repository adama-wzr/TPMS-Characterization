
/*

Some of the GPU related functions.

Last modified: 03/27/2025

Andre Adam.
*/

#ifndef _GPU
#define _GPU

#include "cuda_runtime.h"
#include "cuda.h"
#include "device_launch_parameters.h"
#include <data_structures.hpp>

// CUDA CHECK ERROR

#define CHECK_CUDA(func)                                               \
    {                                                                  \
        cudaError_t status = (func);                                   \
        if (status != cudaSuccess)                                     \
        {                                                              \
            printf("CUDA API failed at line %d with error: %s (%d)\n", \
                   __LINE__, cudaGetErrorString(status), status);      \
            return EXIT_FAILURE;                                       \
        }                                                              \
    }

// Solvers

// 3D GPU SOR

__global__ void JI_SOR3D_kernel(
    float       *A,
    float       *x,
    float       *b,
    float       *xNew,
    long int    nElements,
    int         nCols,
    int         nRows)
{
    /*

        JI_SOR3D_kernel:
        Inputs:
            - pointer to coefficient matrix (device memory)
            - pointer to x-vector (device memory)
            - pointer to b/RHS (device memory)
            - pointer to new x-vector (device memory)
            - long int nElements, passed by value
            - int nCols, passed by value
            - int nRows, passed by value
        Outputs:
            - none.
        
        The function calculates the standard over-relaxed Jacobi Iteration
        for a coefficient matrix with 7 main diagonals. For indexing reference,
        check the discretization function.
    */
    unsigned int myIdx = blockIdx.x * blockDim.x + threadIdx.x;
    float w = 2.0 / 3.0;

    if (myIdx < nElements)
    {
        float sigma = 0;
        for (int j = 1; j < 7; j++)
        {
            if (A[myIdx * 7 + j] > 1e-15)
            {
                if (j == 1)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx - 1];
                }
                else if (j == 2)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx + 1];
                }
                else if (j == 3)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx + nCols];
                }
                else if (j == 4)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx - nCols];
                }
                else if (j == 5)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx + nCols * nRows];
                }
                else if (j == 6)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx - nCols * nRows];
                }
            }
        }
        xNew[myIdx] = (1.0 - w) * x[myIdx] + w / A[myIdx * 7 + 0] * (b[myIdx] - sigma);
    }
}

__global__ void JI_SOR3D_kernelPB(
    float *A,
    float *x,
    float *b,
    float *xNew,
    long int nElements,
    int nCols,
    int nRows,
    int nSlices)
{
    /*

        JI_SOR3D_kernelPB:
        Inputs:
            - pointer to coefficient matrix (device memory)
            - pointer to x-vector (device memory)
            - pointer to b/RHS (device memory)
            - pointer to new x-vector (device memory)
            - long int nElements, passed by value
            - int nCols, passed by value
            - int nRows, passed by value
            - int nSlices, passed by value
        Outputs:
            - none.
        
        The function calculates the standard over-relaxed Jacobi Iteration
        for a coefficient matrix with 7 main diagonals. For indexing reference,
        check the discretization function. This function assumes a periodic
        boundary conditions for the directions parallel to the mass flux.
    */
    unsigned int myIdx = blockIdx.x * blockDim.x + threadIdx.x;
    int mySlice = myIdx / (nCols * nRows);
    int myRow = (myIdx - mySlice * nRows * nCols)/nCols;
    int myCol = myIdx - mySlice * nRows * nCols - myRow * nCols;
    float w = 2.0 / 3.0;

    if (myIdx < nElements)
    {
        float sigma = 0;
        for (int j = 1; j < 7; j++)
        {
            if (A[myIdx * 7 + j] > 1e-15)
            {
                if (j == 1)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx - 1];
                }
                else if (j == 2)
                {
                    sigma += A[myIdx * 7 + j] * x[myIdx + 1];
                }
                else if (j == 3)
                {
                    if (myRow == nRows - 1)
                    {
                        // Periodic South
                        sigma += A[myIdx * 7 + j] * x[mySlice * nRows * nCols + myCol];
                    } else
                    {
                        sigma += A[myIdx * 7 + j] * x[myIdx + nCols];
                    }
                }
                else if (j == 4)
                {
                    if (myRow == 0)
                    {
                        // Periodic North
                        sigma += A[myIdx * 7 + j] * x[mySlice * nRows * nCols + (nRows - 1) * nCols + myCol];
                    } else
                    {
                        sigma += A[myIdx * 7 + j] * x[myIdx - nCols];
                    }   
                }
                else if (j == 5)
                {
                    if (mySlice == nSlices - 1)
                    {
                        // Periodic Back
                        sigma += A[myIdx * 7 + j] * x[myRow * nCols + myCol];
                    } else
                    {
                        sigma += A[myIdx * 7 + j] * x[myIdx + nCols * nRows];
                    }
                    
                }
                else if (j == 6)
                {
                    if (mySlice == 0)
                    {
                        // Periodic Back
                        sigma += A[myIdx * 7 + j] * x[(nSlices - 1) * nRows * nCols + myRow * nCols + myCol];
                    } else
                    {
                        sigma += A[myIdx * 7 + j] * x[myIdx - nCols * nRows];
                    }
                }
            }
        }
        xNew[myIdx] = (1.0 - w) * x[myIdx] + w / A[myIdx * 7 + 0] * (b[myIdx] - sigma);
    }
}

#endif