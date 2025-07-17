/*

Handle everything related to size-distributions.

Last modified 04/03/2025
Andre Adam
*/

#ifndef _SD
#define _SD

#include <data_structures.hpp>
#include <stdlib.h>
#include <omp.h>

int pass12_Global(bool *target_arr, float *EDT, int j, int k, int width, int height, int depth, int primaryPhase)
{
    /*
        Function pass12_Global:
        Inputs:
            - pointer to target_arr, where structure is held
            - pointer to EDT, where the EDT is held
            - j and k give column and slice numbers for row scanning
            - width, height, and depth are the dimensions in number of pixels.
            - the primary phase is the phase relative to which we calculate the distances.
        Outputs:
            - None.
        This function calculates Pass 1 and 2 according to Meijster algorithm.
    */
    int stride = width;
    int offset = k * (width * height) + j;

    if (target_arr[offset] == primaryPhase)
        EDT[offset] = 0;
    else
        EDT[offset] = height + width + depth;

    // scan 1: forward
    for (int i = 1; i < height; i++)
    {
        if (target_arr[offset + i * stride] == primaryPhase)
            EDT[offset + i * stride] = 0;
        else
            EDT[offset + i * stride] = 1 + EDT[offset + (i - 1) * stride];
    }
    // scan 2: backward
    for (int i = height - 2; i >= 0; i--)
    {
        if (EDT[offset + (i + 1) * stride] < EDT[offset + i * stride])
            EDT[offset + i * stride] = 1 + EDT[offset + (i + 1) * stride];
    }
    return 0;
}

int pass34_Global(float *EDT, float *EDT_temp, int scanLength, int stride, int *s, int *t)
{
    /*
        Function pass34_Global:
        Inputs:
            - pointer to EDT, where the EDT is held
            - pointer to EDT_temp, where the old static EDT is held.
            - scan length gives number of pixels to scan
            - stride gives the actual memory stride for scanning single directions
            - pointers to s and t are necessary according to Meijster algorithm
        Outputs:
            - None.
        This function calculates Pass 3 and 4 according to Meijster algorithm. This is repeated (5 and 6)
        for 3D structures.
    */
    for (int i = 0; i < scanLength; i++)
        EDT_temp[i] = EDT[i * stride];

    int q = 0;
    double w = 0;
    s[0] = 0;
    t[0] = 0;

    // forward scan

    for (int u = 1; u < scanLength; u++)
    {
        while (q >= 0 && ((pow((t[q] - s[q]), 2) + pow(EDT_temp[s[q]], 2)) >=
                          (pow((t[q] - u), 2) + pow(EDT_temp[u], 2))))
            q--;

        if (q < 0)
        {
            q = 0;
            s[0] = u;
            t[0] = 0;
        }
        else
        {
            w = 1 + trunc(
                        (pow(u, 2) - pow(s[q], 2) + pow(EDT_temp[u], 2) - pow(EDT_temp[s[q]], 2)) / (2 * (double)(u - s[q])));
            if (w < scanLength)
            {
                q++;
                s[q] = u;
                t[q] = (int)w;
            }
        }
    }

    // backward update
    for (int u = scanLength - 1; u >= 0; u--)
    {
        EDT[u * stride] = sqrt(pow(u - s[q], 2) + pow(EDT_temp[s[q]], 2));
        if (u == t[q])
            q--;
    }

    return 0;
}

void pMeijster3D(bool *targetArray,
                 float *targetEDT,
                 meshInfo *structureInfo,
                 int primaryPhase)
{
    /*
    Function pMeijster3D:
    Inputs:
    - pointer to targetArray, where the structure is held
    - pointer to targetEDT, where the EDT will go.
    - pointer to sizeInfo sruct
    - primaryPhase dictates the phase which the EDT will be calculated in relation to.
    Outputs:
    - None.
    This function calculates the EDT in 3D using parallel computing.
    */

    // To make the periodic BCs we need an array that's twice as big

    int height, width, depth;

    height = structureInfo->numCellsY * 2;
    width = structureInfo->numCellsX * 2;
    depth = structureInfo->numCellsZ * 2;
    long int nElements = structureInfo->nElements * 8;

    // make the new array

    bool *B_ext = (bool *)malloc(sizeof(bool) * nElements);
    memset(B_ext, 0, nElements * sizeof(bool));

    // store offsets

    int realWidth, realHeight, realDepth;
    realWidth = width / 2;
    realHeight = height / 2;
    realDepth = depth / 2;

    // Store old array into the new one

    for (long int index = 0; index < nElements; index++)
    {
        // break down index
        int row, col, slice;

        slice = index/(height * width);
        row = (index - slice * height * width)/width;
        col = index - slice * height * width - row * width;

        // get real index (offset by 50% for periodic boundary)
        int realCol, realRow, realSlice;

        realCol = col - realWidth/2;
        realRow = row - realHeight/2;
        realSlice = slice - realDepth/2;

        // fix out-of-bounds index

        if(realCol < 0)
        {
            realCol = realWidth + realCol;
        }
        else if (realCol > realWidth - 1)
        {
            realCol = realCol - realWidth;
        }

        if (realRow < 0)
        {
            realRow = realHeight + realRow;
        }
        else if( realRow > realHeight - 1)
        {
            realRow  = realRow - realHeight;
        }

        if(realSlice < 0)
        {
            realSlice = realDepth + realSlice;
        }
        else if(realSlice > realDepth - 1)
        {
            realSlice = realSlice - realDepth;
        }

        // Finally assing the new structure

        if(targetArray[realSlice * realWidth * realHeight + realRow * realWidth + realCol] == 1)
        {
            B_ext[index] = 1;
        }

    }

    float *B_EDT = (float *)malloc(sizeof(float) * nElements);
    memset(B_EDT, 0, sizeof(float) * nElements);

#pragma omp parallel
    {
        // initialize s and t locally
        int *s = (int *)malloc(sizeof(int) * width * height);
        int *t = (int *)malloc(sizeof(int) * width * height);
        memset(s, 0, sizeof(int) * width * height);
        memset(t, 0, sizeof(int) * width * height);

        int LL = (height > width) ? height : width;
        LL = (depth > LL) ? depth : LL;

        float *tempEDT = (float *)malloc(sizeof(float) * LL);

        // phase 1

        for (int j = 0; j < width; j++)
        {
#pragma omp for schedule(auto)
            for (int k = 0; k < depth; k++)
            {
                pass12_Global(B_ext, B_EDT, j, k, width, height, depth, primaryPhase);
            }
        }

        // phase 2

        for (int k = 0; k < depth; k++)
        { // all slices
#pragma omp for schedule(auto)
            for (int i = 0; i < height; i++)
            { // all rows
                size_t offset = k * width * height + i * width;
                pass34_Global(B_EDT + offset, tempEDT, width, 1, s, t);
            }
        }

        // maybe not necessary, but clean arrays anyways
        memset(s, 0, sizeof(int) * width * height);
        memset(t, 0, sizeof(int) * width * height);

        for (int j = 0; j < width; j++)
        {
#pragma omp for schedule(auto)
            for (int i = 0; i < height; i++)
            {
                size_t offset = i * width + j;
                pass34_Global(B_EDT + offset, tempEDT, depth, height * width, s, t);
            }
        }
        // memory management inside parallel loop
        free(s);
        free(t);
        free(tempEDT);
    }

    // deconstruct B_EDT into target EDT

    for(int index = 0; index < nElements/8; index++)
    {
        // get index for real array
        int col, row, slice;
        slice = index / (realHeight * realWidth);
        row = (index - slice * realHeight * realWidth) / realWidth;
        col = index - slice * realHeight * realWidth - row * realWidth;

        // tranform it into the index for the extended array
        int extCol, extRow, extSlice;
        extCol = col + realWidth/2;
        extRow = row + realHeight/2;
        extSlice = slice + realDepth/2;

        // Transfer EDT
        targetEDT[index] = B_EDT[extSlice * width * height + extRow * width + extCol];
    }

    // Memory management
    free(B_EDT);
    free(B_ext);

    return;
}

int partSD_3D(options* opts, meshInfo *mesh, saveInfo *save, char *P, int POI)
{
    /*
        partSD_3D:
        - 
    */

    if(opts->verbose)
        printf("Particle-Size Distribution:\n");
    
    // set num threads
    omp_set_num_threads(opts->nThreads);

    // Loop variables

    long int p_sum, d_sum, e_sum;

    e_sum = 1;      // initialized for while loop

    // Index 0 = p, index 1 = d, index 2 = e

    long int *PDE_sum = (long int *)malloc(sizeof(long int) * opts->maxR * 3);

    memset(PDE_sum, 0, sizeof(long int) * opts->maxR * 3);

    // Array for storirng radius labels

    int *R = (int *)malloc(sizeof(int) * mesh->nElements);

    for(int i = 0; i < mesh->nElements; i++)
        R[i] = -1;

    // arrays for holding the EDT

    float *EDT_D = (float *)malloc(sizeof(float) * mesh->nElements);
    float *EDT_E = (float *)malloc(sizeof(float) * mesh->nElements);
    
    // Arrays for holding the binary arrays, B, D, and E

    bool *E = (bool *)malloc(sizeof(bool) * mesh->nElements);
    bool *D = (bool *)malloc(sizeof(bool) * mesh->nElements);
    bool *B = (bool *)malloc(sizeof(bool) * mesh->nElements);

    memset(E, 0, sizeof(bool) * mesh->nElements);
    memset(D, 0, sizeof(bool) * mesh->nElements);
    memset(B, 0, sizeof(bool) * mesh->nElements);

    // if P[i] = POI, B[i] = 1
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < mesh->nElements; i++)
    {
        if(P[i] == POI)
            B[i] = 1;
    }

    // EDT for dilation only calculated once
    pMeijster3D(B, EDT_D, mesh, 0);
    int radius = 1;

    // Main Loop

    while (e_sum != 0 && radius <= opts->maxR)
    {
        // copy P into D (probably not necessary)

        memcpy(D, B, sizeof(bool) * mesh->nElements);
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < mesh->nElements; i++)
        {
            if (pow(EDT_D[i], 2) <= (float)radius * radius)
                D[i] = 0;
        }

        // Copy D into E

        memcpy(E, D, sizeof(bool) * mesh->nElements);

        // Meijster in D

        pMeijster3D(D, EDT_E, mesh, 1);

// Update E
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < mesh->nElements; i++)
        {
            if (pow(EDT_E[i], 2) <= (float)radius * radius)
                E[i] = 1;
        }

        // Quantify changes

        e_sum = 0;
        d_sum = 0;
        p_sum = 0;

#pragma omp parallel for reduction(+ : p_sum, d_sum, e_sum)
        for (int i = 0; i < mesh->nElements; i++)
        {
            p_sum += B[i];
            d_sum += D[i];
            e_sum += E[i];

            if (B[i] - E[i] == 1 && R[i] == -1)
                R[i] = radius;
        }

        PDE_sum[(radius - 1 ) * 3 + 0] = p_sum;
        PDE_sum[(radius - 1 ) * 3 + 1] = d_sum;
        PDE_sum[(radius - 1 ) * 3 + 2] = e_sum;

        // print verbose

        if (opts->verbose)
            printf("R = %d, P = %ld, E = %ld, D = %ld\n", radius, p_sum, e_sum, d_sum);

        // increment radius
        radius++;
    }

    int lastR = radius;

    // calculate partSD and print to output file

    long int sum_removed = 0;
    double *partRemoved = (double *)malloc(sizeof(double) * lastR);

    // get particles removed at R = 1

    partRemoved[0] = PDE_sum[0 * 3 + 0] - PDE_sum[0 * 3 + 2];
    sum_removed += (int)partRemoved[0];

    for (int i = 1; i < lastR; i++)
    {
        partRemoved[i] = PDE_sum[(i - 1) * 3 + 2] - PDE_sum[i * 3 + 2];
        sum_removed += (int)partRemoved[i];
    }

    FILE *partSD_OUT = fopen(opts->partSDOut, "w+");

    fprintf(partSD_OUT, "r,p(r)\n");
    for (int i = 0; i < (lastR ); i++)
    {
        fprintf(partSD_OUT, "%d,%lf\n", i + 1, (double)partRemoved[i] / sum_removed);
    }

    fclose(partSD_OUT);

    // Calculate D50
    float p, D50;
    D50 = 0;
    for(int i = 0; i < lastR; i++)
    {
        p = (double)partRemoved[i] / sum_removed;
        D50 += p*(2*(i + 1));
    }

    save->part50 = (float)D50/mesh->numCellsX;

    // char filename[100];

    // sprintf(filename, "testR.csv");

    // // save labels

    // saveLabels3D(R, mesh, filename);

    // memory management

    free(EDT_D);
    free(EDT_E);

    free(B);
    free(E);
    free(D);
    free(R);
    return 0;
}

int poreSD_3D(options* opts, meshInfo *mesh, saveInfo *save, char *P, int POI)
{
    /*
        poreSD_3D:
        - 
    */

    if(opts->verbose)
        printf("Pore-Size Distribution:\n");

    // Loop variables

    // set num threads

    omp_set_num_threads(opts->nThreads);

    long int p_sum, d_sum, e_sum;

    e_sum = 1;      // initialized for while loop

    // Index 0 = p, index 1 = d, index 2 = e

    long int *PDE_sum = (long int *)malloc(sizeof(long int) * opts->maxR * 3);

    memset(PDE_sum, 0, sizeof(long int) * opts->maxR * 3);

    // Array for storirng radius labels

    int *R = (int *)malloc(sizeof(int) * mesh->nElements);

    for(int i = 0; i < mesh->nElements; i++)
        R[i] = -1;

    // arrays for holding the EDT

    float *EDT_D = (float *)malloc(sizeof(float) * mesh->nElements);
    float *EDT_E = (float *)malloc(sizeof(float) * mesh->nElements);
    
    // Arrays for holding the binary arrays, B, D, and E

    bool *E = (bool *)malloc(sizeof(bool) * mesh->nElements);
    bool *D = (bool *)malloc(sizeof(bool) * mesh->nElements);
    bool *B = (bool *)malloc(sizeof(bool) * mesh->nElements);

    memset(E, 0, sizeof(bool) * mesh->nElements);
    memset(D, 0, sizeof(bool) * mesh->nElements);
    memset(B, 0, sizeof(bool) * mesh->nElements);

    // if P[i] = POI, B[i] = 1
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < mesh->nElements; i++)
    {
        if(P[i] == POI)
            B[i] = 1;
    }

    // EDT for dilation only calculated once
    pMeijster3D(B, EDT_D, mesh, 0);
    int radius = 1;

    // Main Loop

    while (e_sum != 0 && radius <= opts->maxR)
    {
        // copy P into D (probably not necessary)

        memcpy(D, B, sizeof(bool) * mesh->nElements);
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < mesh->nElements; i++)
        {
            if (pow(EDT_D[i], 2) <= (float)radius * radius)
                D[i] = 0;
        }

        // Copy D into E

        memcpy(E, D, sizeof(bool) * mesh->nElements);

        // Meijster in D

        pMeijster3D(D, EDT_E, mesh, 1);

// Update E
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < mesh->nElements; i++)
        {
            if (pow(EDT_E[i], 2) <= (float)radius * radius)
                E[i] = 1;
        }

        // Quantify changes

        e_sum = 0;
        d_sum = 0;
        p_sum = 0;

#pragma omp parallel for reduction(+ : p_sum, d_sum, e_sum)
        for (int i = 0; i < mesh->nElements; i++)
        {
            p_sum += B[i];
            d_sum += D[i];
            e_sum += E[i];

            if (B[i] - E[i] == 1 && R[i] == -1)
                R[i] = radius;
        }

        PDE_sum[(radius - 1 ) * 3 + 0] = p_sum;
        PDE_sum[(radius - 1 ) * 3 + 1] = d_sum;
        PDE_sum[(radius - 1 ) * 3 + 2] = e_sum;

        // print verbose

        if (opts->verbose)
            printf("R = %d, P = %ld, E = %ld, D = %ld\n", radius, p_sum, e_sum, d_sum);

        // increment radius
        radius++;
    }

    int lastR = radius;

    // calculate partSD and print to output file

    long int sum_removed = 0;
    double *partRemoved = (double *)malloc(sizeof(double) * lastR);

    // get particles removed at R = 1

    partRemoved[0] = PDE_sum[0 * 3 + 0] - PDE_sum[0 * 3 + 2];
    sum_removed += (int)partRemoved[0];

    for (int i = 1; i < lastR; i++)
    {
        partRemoved[i] = PDE_sum[(i - 1) * 3 + 2] - PDE_sum[i * 3 + 2];
        sum_removed += (int)partRemoved[i];
    }

    FILE *poreSDOut = fopen(opts->poreSDOut, "w+");

    fprintf(poreSDOut, "r,p(r)\n");
    for (int i = 0; i < (lastR ); i++)
    {
        fprintf(poreSDOut, "%d,%lf\n", i + 1, (double)partRemoved[i] / sum_removed);
    }

    fclose(poreSDOut);

    // Calculate D50
    float p, D50;
    D50 = 0;
    for(int i = 0; i < lastR; i++)
    {
        p = (double)partRemoved[i] / sum_removed;
        D50 += p*(2*(i + 1));
    }

    save->pore50 = (float)D50/mesh->numCellsX;

    // char filename[100];

    // sprintf(filename, "test_pore_R.csv");

    // // save labels

    // saveLabels3D(R, mesh, filename);


    // memory management

    free(EDT_D);
    free(EDT_E);

    free(B);
    free(E);
    free(D);
    free(R);
    return 0;
}


#endif