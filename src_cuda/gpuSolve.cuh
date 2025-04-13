#ifndef _GPU_SOLVE
#define _GPU_SOLVE

#include <lib/data_structures.hpp>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <GPU_kernels.cuh>
#include <cuda_structs.cuh>

#define MAX_GPU 32

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

int initGPU_3DSOR(float **d_Coeff,
                  float **d_RHS,
                  float **d_Conc,
                  float **d_ConcTemp,
                  meshInfo *mesh)
{
    /*
        Function initGPU_3DSOR:
        Inputs:
            - float pointer to d_Coeff, storing coeff matrix in GPU
            - float pointer to d_RHS, storing RHS vector in GPU
            - float pointer to d_Conc, the concentration array in GPU memory
            - float pointer to d_ConcTemp, where the concentration array will
                be modified in GPU memory
            - pointer to meshInfo, holding general information about the mesh.
        Outputs:
            - None.

        The function will allocate the sufficient space for the arrays needed for
        the Standard Over-Relaxed Jacobi Method. It also initializes the arrays.
        Error calls are returned if it fails.
    */

    // Set device

    CHECK_CUDA(cudaSetDevice(0));

    // Allocate space

    CHECK_CUDA(cudaMalloc((void **)&(*d_Coeff), mesh->nElements * sizeof(float) * 7));
    CHECK_CUDA(cudaMalloc((void **)&(*d_RHS), mesh->nElements * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void **)&(*d_Conc), mesh->nElements * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void **)&(*d_ConcTemp), mesh->nElements * sizeof(float)));

    // Set buffers

    CHECK_CUDA(cudaMemset((*d_Coeff), 0, mesh->nElements * sizeof(float) * 7));
    CHECK_CUDA(cudaMemset((*d_RHS), 0, mesh->nElements * sizeof(float)));
    CHECK_CUDA(cudaMemset((*d_Conc), 0, mesh->nElements * sizeof(float)));
    CHECK_CUDA(cudaMemset((*d_ConcTemp), 0, mesh->nElements * sizeof(float)));

    return 0;
}

int JI3D_SOR(float *Coeff,
             float *RHS,
             float *Concentration,
             float *d_Coeff,
             float *d_RHS,
             float *d_Conc,
             float *d_ConcTemp,
             options *opts,
             meshInfo *mesh)
{
    /*
        Function JI3D_SOR:
        Inputs:
            - pointer to coefficient matrix array
            - pointer to RHS matrix array
            - pointer to Concentration distribution array
            - pointer to device coefficient matrix
            - pointer to device right-hand side array
            - pointer to device concentration array
            - pointer to device temporary concentration array storage
            - pointer to options struct
            - pointer to mesh struct
        Outputs:
            - None

        This function will manage the host-device interactions for the Jacobi Iteration method
        in 3D, with a standard over-relaxation applied. The function will manage data transfers,
        convergence criteria, and kernel coordination.
    */

    long int iterCount = 0;
    int threads_per_block = 128;
    int numBlocks = mesh->nElements / threads_per_block + 1;

    float pctChange = 1;
    int iterToCheck = 1000;

    // copy arrays into GPU

    CHECK_CUDA(cudaMemcpy(d_Conc, Concentration,
                          sizeof(float) * mesh->nElements, cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMemcpy(d_ConcTemp, Concentration,
                          sizeof(float) * mesh->nElements, cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMemcpy(d_RHS, RHS,
                          sizeof(float) * mesh->nElements, cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMemcpy(d_Coeff, Coeff,
                          sizeof(float) * mesh->nElements * 7, cudaMemcpyHostToDevice));

    // Create Array to store temp Conc

    float *TempConc = (float *)malloc(sizeof(float) * mesh->nElements);

    memcpy(TempConc, Concentration, sizeof(float) * mesh->nElements);

    // start the main loop

    while (iterCount < opts->MAX_ITER && pctChange > opts->ConvergeCriteria)
    {
        // call kernel
        if (opts->PB)
        {
            JI_SOR3D_kernelPB<<<numBlocks, threads_per_block>>>(d_Coeff, d_ConcTemp, d_RHS, d_Conc,
                                                                mesh->nElements, mesh->numCellsX, mesh->numCellsY, mesh->numCellsZ);
        }
        else
        {
            JI_SOR3D_kernel<<<numBlocks, threads_per_block>>>(d_Coeff, d_ConcTemp, d_RHS, d_Conc,
                                                              mesh->nElements, mesh->numCellsX, mesh->numCellsY);
        }

        CHECK_CUDA(cudaGetLastError());

        // check convergence
        if (iterCount % iterToCheck == 0 && iterCount != 0)
        {
            // copy array from device to host
            CHECK_CUDA(cudaMemcpy(Concentration, d_Conc, sizeof(float) * mesh->nElements, cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();
            // compare
            float sum = 0;
            long int count = 0;

            for (int i = 0; i < mesh->nElements; i++)
            {
                if (Concentration[i] != 0)
                {
                    sum += fabs((Concentration[i] - TempConc[i]) / Concentration[i]);
                    count++;
                }
            }
            // calculate the change
            pctChange = sum / count;
            // copy memory to temp conc
            memcpy(TempConc, Concentration, sizeof(float) * mesh->nElements);
        }

        if (iterCount % 10000 == 0 && opts->verbose == 1)
        {
            printf("Iter %ld, pct Change = %lf\n", iterCount, pctChange);
        }

        // update d_Conc = d_ConcTemp

        CHECK_CUDA(cudaMemcpy(d_ConcTemp, d_Conc, sizeof(float) * mesh->nElements, cudaMemcpyDeviceToDevice));

        // increment
        iterCount++;
    }

    // copy the solution

    CHECK_CUDA(cudaMemcpy(Concentration, d_ConcTemp,
                          sizeof(float) * mesh->nElements, cudaMemcpyDeviceToHost));

    // print success

    if (opts->verbose)
    {
        printf("Total iter = %ld, pct change = %lf\n", iterCount, pctChange);
    }

    // store info to print

    mesh->conv = pctChange;
    mesh->iterCount = iterCount;

    return 0;
}

int JI3D_SOR_multi(float *Coeff,
                   float *RHS,
                   float *Concentration,
                   MGPU_Solve *gpu,
                   options *opts,
                   meshInfo *mesh,
                   int nDevices)
{
    /*
        Function JI3D_SOR_multi:
        Inputs:
        - pointer to coefficient matrix array
        - pointer to RHS matrix array
        - pointer to Concentration distribution array
        - pointer to GPU solver struct
        - pointer to options struct
        - pointer to mesh struct
        - number of GPUs
        Outputs:
        - None

        This function will manage the host-device interactions for the Jacobi Iteration method
        in 3D, with a standard over-relaxation applied. The function will manage data transfers,
        convergence criteria, and kernel coordination.
    */

    long int iterCount = 0;
    int threads_per_block = 128;
    // int numBlocks = mesh->nElements / threads_per_block + 1;
    
    for(int i = 0; i < nDevices; i++)
    {
        gpu[i].nBlocks = gpu[i].dataN / threads_per_block + 1;
    }

    float pctChange = 1;
    int iterToCheck = 1000;

    // copy arrays into GPU (sync)

    for (int i = 0; i < nDevices; i++)
    {
        // Set device
        CHECK_CUDA(cudaSetDevice(i));

        // Copy arrays asynchronously
        CHECK_CUDA(
            cudaMemcpyAsync(gpu[i].d_Coeff, Coeff + (gpu[i].xOffset * 7), gpu[i].dataN * 7 * sizeof(float), 
                            cudaMemcpyHostToDevice, gpu[i].stream)
        );
        CHECK_CUDA(
            cudaMemcpyAsync(gpu[i].d_RHS, RHS + gpu[i].xOffset, gpu[i].dataN * sizeof(float), 
                            cudaMemcpyHostToDevice, gpu[i].stream)
        );
        CHECK_CUDA(
            cudaMemcpyAsync(gpu[i].d_Temp, Concentration + gpu[i].xOffset, gpu[i].dataN * sizeof(float), 
                            cudaMemcpyHostToDevice, gpu[i].stream)
        );
        CHECK_CUDA(
            cudaMemcpyAsync(gpu[i].d_X, Concentration, mesh->nElements * sizeof(float), 
                            cudaMemcpyHostToDevice, gpu[i].stream)
        );
    }

    // Create Array to store temp Conc

    float *TempConc = (float *)malloc(sizeof(float) * mesh->nElements);

    memcpy(TempConc, Concentration, sizeof(float) * mesh->nElements);

    // start the main loop

    while (iterCount < opts->MAX_ITER && pctChange > opts->ConvergeCriteria)
    {
        // call kernel
        if (opts->PB)
        {
            for(int i = 0; i < nDevices; i++)
            {
                // set device
                CHECK_CUDA(cudaSetDevice(i));
                // make sure all streams have finished
                cudaStreamSynchronize(gpu[i].stream);
                // launch kernel async
                JI_SOR3D_PB_multi<<<gpu[i].nBlocks, threads_per_block, 0, gpu[i].stream>>>(
                    gpu[i].d_Coeff, gpu[i].d_X, gpu[i].d_RHS, gpu[i].d_Temp,
                    mesh->nElements, gpu[i].xOffset, mesh->numCellsX, mesh->numCellsY,
                    mesh->numCellsZ, gpu[i].dataN
                );
            }
        }
        else
        {
            printf("Not Implemented Yet\n");
        }

        CHECK_CUDA(cudaGetLastError());

        // check convergence

        if (iterCount % iterToCheck == 0 && iterCount != 0)
        {
            for(int i = 0; i < nDevices; i++)
            {
                // set device
                CHECK_CUDA(cudaSetDevice(i));

                // Copy data asynchronously
                CHECK_CUDA(
                    cudaMemcpyAsync(Concentration + gpu[i].xOffset, gpu[i].d_Temp, sizeof(float) * gpu[i].dataN, cudaMemcpyDeviceToHost, gpu[i].stream));
            }
            // compare
            float sum = 0;
            long int count = 0;

            for(int i = 0; i < nDevices; i++)
            {
                // Wait for all processes to finish
                cudaStreamSynchronize(gpu[i].stream);
            }

            for (int i = 0; i < mesh->nElements; i++)
            {
                if (Concentration[i] != 0)
                {
                    sum += fabs((Concentration[i] - TempConc[i]) / Concentration[i]);
                    count++;
                }
            }
            // calculate the change
            pctChange = sum / count;
            // copy memory to temp conc
            memcpy(TempConc, Concentration, sizeof(float) * mesh->nElements);
        }

        if (iterCount % 10000 == 0 && opts->verbose == 1)
        {
            printf("Iter %ld, pct Change = %lf\n", iterCount, pctChange);
        }

        // update d_Conc = d_ConcTemp

        for (int i = 0; i < nDevices; i++)
        {
            // Set device
            CHECK_CUDA(cudaSetDevice(i));
            // make sure all streams have finished
            cudaStreamSynchronize(gpu[i].stream);
            // copy necessary data
            CHECK_CUDA(
                cudaMemcpyAsync(Concentration + gpu[i].xOffset, gpu[i].d_Temp, sizeof(float) * gpu[i].dataN, cudaMemcpyDeviceToHost, gpu[i].stream)
            );
        }

        // copy of Concentration to device

        for (int i = 0; i < nDevices; i++)
        {
            // Set device
            CHECK_CUDA(cudaSetDevice(i));
            // make sure all streams have finished
            cudaStreamSynchronize(gpu[i].stream);
            // copy necessary data
            CHECK_CUDA(
                cudaMemcpyAsync(gpu[i].d_X, Concentration, sizeof(float)*mesh->nElements, cudaMemcpyHostToDevice, gpu[i].stream)
            );
        }

        // increment
        iterCount++;
    }

    // copy the solution

    for(int i = 0; i<nDevices; i++)
    {
        CHECK_CUDA(cudaSetDevice(i));
        // Wait for all processes to finish
        cudaStreamSynchronize(gpu[i].stream);
        // copy result to host
        CHECK_CUDA(
            cudaMemcpy(Concentration + gpu[i].xOffset, gpu[i].d_Temp, sizeof(float)*gpu[i].dataN, cudaMemcpyDeviceToHost)
        );
    }

    // print success

    if (opts->verbose)
    {
        printf("Total iter = %ld, pct change = %lf\n", iterCount, pctChange);
    }

    // store info to print

    mesh->conv = pctChange;
    mesh->iterCount = iterCount;

    // Memory management
    free(TempConc);

    return 0;
}

int unInitGPU_SOR(float **d_Coeff,
                  float **d_RHS,
                  float **d_Conc,
                  float **d_ConcTemp)
{
    /*
        Function unInitGPU_SOR:
        Inputs:
            - float pointer to d_Coeff, storing coeff matrix in GPU
            - float pointer to d_RHS, storing RHS vector in GPU
            - float pointer to d_Conc, the concentration array in GPU memory
            - float pointer to d_ConcTemp, where the concentration array will
                be modified in GPU memory
        Outputs:
            - None.

        The function will free space in device memory.
    */

    CHECK_CUDA(cudaFree((*d_Coeff)));
    CHECK_CUDA(cudaFree((*d_RHS)));
    CHECK_CUDA(cudaFree((*d_Conc)));
    CHECK_CUDA(cudaFree((*d_ConcTemp)));

    return 0;
}

int gpuHandler(options *opts,
               meshInfo *mesh,
               saveInfo *save,
               float *Concentration,
               float *CoeffMatrix,
               float *RHS)
{
    /*
        Function gpuHandler:
        Inputs:
            - pointer to options struct
            - pointer to mesh struct
            - pointer to save struct
            - pointer to Concentration matrix
            - pointer to coefficient matrix
            - pointer to RHS
        Outputs:
            - None

        Function will handle GPU soluiton, which includes reading user options,
        calling correct solver, returning final concentration solution, and
        handling memory allocations.
    */

    // Figure out how many GPUs to use

    int nDevices;

    cudaGetDeviceCount(&nDevices);

    if (nDevices == 0)
    {
        printf("Error, CUDA could not find a suitable device.\n");
        printf("Exiting now!\n");
        return 1;
    }
    else if (nDevices < opts->nGPU)
    {
        printf("Requested %d GPUs, only %d found. Proceeding...\n", opts->nGPU, nDevices);
    }
    else if (nDevices > opts->nGPU)
    {
        nDevices = opts->nGPU;
    }

    /*
        Currently only support one GPU, so that's a little useless
    */

    if (nDevices == 1) // single device execution
    {
        // pointers to GPU arrays

        float *d_Coeff = NULL;
        float *d_RHS = NULL;
        float *d_Conc = NULL;
        float *d_ConcTemp = NULL;

        // Intialize GPU Arrays

        initGPU_3DSOR(&d_Coeff, &d_RHS, &d_Conc, &d_ConcTemp, mesh);

        // Solve

        JI3D_SOR(CoeffMatrix, RHS, Concentration, d_Coeff,
                 d_RHS, d_Conc, d_ConcTemp, opts, mesh);

        // un-initialize
        unInitGPU_SOR(&d_Coeff, &d_RHS, &d_Conc, &d_ConcTemp);
    }
    else
    {
        if (opts->verbose)
            printf("MUlti-GPU Solve: %d Devices Detected\n", nDevices);
        // Multi-GPU Execution
        MGPU_Solve *gpu = (MGPU_Solve *)malloc(sizeof(MGPU_Solve) * nDevices);

        // Divide the data for each GPU

        for(int i = 0; i < nDevices; i++)
        {
            gpu[i].dataN = mesh->nElements / nDevices;
        }

        // take into account "odd" data sizes
        for (int i = 0; i < mesh->nElements % nDevices; i++)
        {
            gpu[i].dataN++;
        }

        // Get offsets
        long int offset = 0;
        for (int i = 0; i < nDevices; i++)
        {
            if (i == 0)
            {
                gpu[i].xOffset = 0;
                offset = gpu[i].dataN;
            }
            else
            {
                gpu[i].xOffset = offset;
                offset += gpu[i].dataN;
            }
        }

        // Create Streams and allocate memory
        for (int i = 0; i < nDevices; i++)
        {
            // set device
            CHECK_CUDA(cudaSetDevice(i));
            // Create Stream
            CHECK_CUDA(cudaStreamCreate(&gpu[i].stream));
            // Allocate the necessary arrays
            CHECK_CUDA(
                cudaMalloc((void **)&gpu[i].d_Coeff, gpu[i].dataN * 7 * sizeof(float)));
            CHECK_CUDA(
                cudaMalloc((void **)&gpu[i].d_RHS, gpu[i].dataN * sizeof(float)));
            CHECK_CUDA(
                cudaMalloc((void **)&gpu[i].d_Temp, gpu[i].dataN * sizeof(float)));
            CHECK_CUDA(
                cudaMalloc((void **)&gpu[i].d_X, mesh->nElements * sizeof(float)));
        }

        // Allocate page-locked memory

        float *h_Coeff, *h_RHS, *h_Conc;

        CHECK_CUDA(
            cudaMallocHost((void **)&h_Coeff, mesh->nElements * 7 * sizeof(float))
        );
        CHECK_CUDA(
            cudaMallocHost((void **)&h_RHS, mesh->nElements * sizeof(float))
        );
        CHECK_CUDA(
            cudaMallocHost((void **)&h_Conc, mesh->nElements * sizeof(float))
        );
        
        // Copy arrays to page-locked memory
        memcpy(h_Coeff, CoeffMatrix, sizeof(float) * mesh->nElements * 7);
        memcpy(h_RHS, RHS, sizeof(float) * mesh->nElements);
        memcpy(h_Conc, Concentration, sizeof(float) * mesh->nElements);

        // Call Solver

        JI3D_SOR_multi(h_Coeff, h_RHS, h_Conc, gpu, opts, mesh, nDevices);

        // GPU Memory Management

        for (int i = 0; i < nDevices; i++)
        {
            // Set device
            CHECK_CUDA(cudaSetDevice(i));
            // free gpu array pointers
            CHECK_CUDA(
                cudaFree(gpu[i].d_Coeff));
            CHECK_CUDA(
                cudaFree(gpu[i].d_RHS));
            CHECK_CUDA(
                cudaFree(gpu[i].d_X));
            CHECK_CUDA(
                cudaFree(gpu[i].d_Temp));
            // Stream Destroy
            CHECK_CUDA(
                cudaStreamDestroy(gpu[i].stream));
        }

        // Free the cuda struct
        free(gpu);
    }

    return 0;
}

#endif