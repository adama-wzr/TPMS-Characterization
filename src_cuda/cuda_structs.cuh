
#ifndef _CUDA_STRUC
#define _CUDA_STRUC

#include "cuda.h"
#include "cuda_runtime.h"



// GPU Context (multiGPU)

typedef struct
{
    long int dataN;     // number of data points
    // device buffers
    float *d_Coeff, *d_Temp, *d_RHS, *d_X;
    // Stream for asynchronous execution
    cudaStream_t stream;
    // offset in data storage
    long int xOffset;
    // num blocks
    int nBlocks;
} MGPU_Solve;



#endif