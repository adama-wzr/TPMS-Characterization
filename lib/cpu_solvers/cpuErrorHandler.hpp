#ifndef _CPU_ERR

#define _CPU_ERR

#include <iostream>
#include <data_structures.hpp>

int gpuHandler(options *opts,
               meshInfo *mesh,
               saveInfo *save,
               float *Concentration,
               float *CoeffMatrix,
               float *RHS)
{
    printf("------------------------------------------------------\n\n");
    printf("    User Selected GPU Solver: no CUDA Found\n");
    printf("    Please set GPU Flag to false! Exiting....\n");
    printf("\n\n------------------------------------------------------\n\n");
    return 1;
}

#endif
