#ifndef _HELPER
#define _HELPER

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <stdbool.h>
#include <fstream>
#include <cfloat>
#include <set>
#include <tuple>
#include <string>
#include <iostream>
#include "cuda_runtime.h"
#include "cuda.h"
#include <omp.h>
#define PI 3.14159265

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

/*

    Data Structure Definitions:

*/

// Handling user input

typedef struct
{
    float CLeft;            // Concentration of trace species in left boundary
    float CRight;           // Concentration of trace species in right boundary
    long int MAX_ITER;      // Max iterations
    float ConvergeCriteria; // Convergence Criteria
    int printOut;           // Flag to print output or not
    char *outputFilename;   // Output filename
    int verbose;            // verbose flag
    int nVoxels;            // height in number of pixels
    int nThreads;           // number of threads
    int useGPU;             // Use GPU or not?
    int nGPU;               // number of GPUs
    int TPMS_Type;          // type of TPMS structure
    float isoValues;        // iso-bounds "b" for TPMS thickness
    int Tau;                // run tortuosity or not
    int PB;                 // periodic BC
    int poreSD;             // pore-size distribution
    char *poreSDOut;        // pore size distribution file name
    int partSD;             // particle size distribution
    char *partSDOut;        // particle size distribution file name
    int maxR;               // maximum scan radius;
} options;

// Mesh related information

typedef struct
{
    int numCellsX;
    int numCellsY;
    int numCellsZ;
    long int nElements;
    float dx;
    float dy;
    float dz;
    long int iterCount;
    float conv;
    float VF; // solid volume fraction
} meshInfo;

// Array holding output data

typedef struct
{
    float porosity;
    float ePorosity;
    float Deff_TH_MAX;
    float Deff;
    float Tau;
    float SA;
    long int nElements;
    float pore50;
    float part50;
    int nVoxel;
} saveInfo;


// Define coords for Flood Fill

typedef std::tuple<int, int, int> coord;


/*

    GPU Kernels

*/

// 3D GPU SOR

__global__ void JI_SOR3D_kernel(
    float *A,
    float *x,
    float *b,
    float *xNew,
    long int nElements,
    int nCols,
    int nRows)
{
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

/*

    Handling User Input:

*/

void printOpts(options *opts)
{
    /*
        printOpts Function:
        Inputs:
            - pointer to options struct
        Outputs:
            - none.

        Function will simply print some of the user entered options
        to the command line.
    */

    printf("--------------------------------------\n\n");
    printf("TPMS Simulation\n");
    printf("Current selected options:\n\n");
    printf("--------------------------------------\n");
    printf("TPMS Type = %d\n", opts->TPMS_Type);
    printf("Dimension size = %d voxels\n", opts->nVoxels);
    printf("Number of Threads = %d\n", opts->nThreads);

    if (opts->Tau)
        printf("Running Tortuosity\n");

    if(opts->PB)
        printf("BC: Periodic\n");
    else
        printf("BC: No Flux [default]\n");
    


    return;
}

void readInputGeneral(char *filename, options *opts)
{
    /*
        readInputGeneral Function:
        Inputs:
            - FileName: pointer to where the input file name is stored.
            - struct options: pass a struct with the options.
        Outputs: None

        Function reads the input file and stores the options in the opts struct.
    */

    // initiate necessary variables for input reading
    std::string myText;

    char tempC[1000];
    double tempD;
    char tempFilenames[1000];
    std::ifstream InputFile(filename);

    opts->outputFilename = (char *)malloc(1000 * sizeof(char));
    opts->partSDOut = (char *)malloc(1000 * sizeof(char));
    opts->poreSDOut = (char *)malloc(1000 * sizeof(char));

    // default values

    opts->TPMS_Type = 1;
    opts->Tau = 0;
    opts->PB = 0;
    opts->poreSD = 0;
    opts->partSD = 0;
    opts->CLeft = 0;
    opts->CRight = 1;
    opts->maxR = 1000;
    opts->nThreads = 1;

    // read options

    while (std::getline(InputFile, myText))
    {
        sscanf(myText.c_str(), "%s %lf", tempC, &tempD);
        if (strcmp(tempC, "OutputName:") == 0)
        {
            sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
            strcpy(opts->outputFilename, tempFilenames);
        }
        else if (strcmp(tempC, "Verbose:") == 0)
        {
            opts->verbose = (int)tempD;
        }
        else if (strcmp(tempC, "nVoxels:") == 0)
        {
            opts->nVoxels = (int)tempD;
        }
        else if (strcmp(tempC, "isoValue:") == 0)
        {
            opts->isoValues = (float)tempD;
        }
        else if (strcmp(tempC, "TPMS:") == 0)
        {
            opts->TPMS_Type = (int)tempD;
        }
        else if (strcmp(tempC, "Tau:") == 0)
        {
            opts->Tau = (int)tempD;
        }
        else if (strcmp(tempC, "poreSD:") == 0)
        {
            opts->poreSD = (int)tempD;
        }
        else if (strcmp(tempC, "partSD:") == 0)
        {
            opts->partSD = (int)tempD;
        }
        else if (strcmp(tempC, "pb:") == 0)
        {
            opts->PB = (int)tempD;
        }
        else if (strcmp(tempC, "nThreads:") == 0)
        {
            opts->nThreads = (int)tempD;
        }
        else if (strcmp(tempC, "Convergence:") == 0)
        {
            opts->ConvergeCriteria = tempD;
        }
        else if (strcmp(tempC, "MaxIter:") == 0)
        {
            opts->MAX_ITER = (long int) tempD;
        }
        else if (strcmp(tempC, "partSDOut:") == 0)
        {
            sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
            strcpy(opts->partSDOut, tempFilenames);
        }
        else if (strcmp(tempC, "poreSDOut:") == 0)
        {
            sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
            strcpy(opts->poreSDOut, tempFilenames);
        }
    }

    return;
}

/*

    Handling output:

*/

void outSA(options *opts, saveInfo *save)
{
    /*
        outSA Function:
        Inputs:
            - pointer to user options struct
            - pointer to saveInfo struct
        Output:
            - None

        Function will print the output from running the code to a output file
        taking in consideration the user options.
    */

    bool headerFlag = true;

    // Check if file exists

    if (FILE *TEST = fopen(opts->outputFilename, "r"))
    {
        fclose(TEST);
        headerFlag = false;
    }

    // Open file

    FILE *OUT = fopen(opts->outputFilename, "a+");

    if (headerFlag)
    {
        fprintf(OUT, "Structure,nElement,iso,pore,SA\n");
    }

    // print output from inputs

    fprintf(OUT, "%d,%d,%1.3e,%1.3e,%1.3e\n", opts->TPMS_Type, save->nVoxel,
            opts->isoValues, save->porosity, save->SA);

    fclose(OUT);
    return;
}

void outputGeneral(options *opts, saveInfo *save)
{
    /*
        outputGeneral Function:
        Inputs:
            - pointer to user options struct
            - pointer to saveInfo struct
        Output:
            - None

        Function will print the output from running the code to a output file
        taking in consideration the user options. Handles general output.
    */

    bool headerFlag = true;

    // Check if file exists

    if (FILE *TEST = fopen(opts->outputFilename, "r"))
    {
        fclose(TEST);
        headerFlag = false;
    }

    // Open file

    FILE *OUT = fopen(opts->outputFilename, "a+");

    if (headerFlag)
    {
        fprintf(OUT, "Structure,nElement,iso,pore,SA,"
                     "DeffMax,Deff,Tau,pore50,part50\n");
    }

    // print output from inputs

    fprintf(OUT, "%d,%d,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e\n",
        opts->TPMS_Type, save->nVoxel, opts->isoValues,
        save->porosity, save->SA, save->Deff_TH_MAX,
        save->Deff, save->Tau, save->pore50, save->part50);

    fclose(OUT);
    return;
}

/*

    TPMS Definitions:

*/

int C_Y(char *P, float b, int numCellsX, meshInfo *info)
{
    /*
        Function C_Y3D:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - int numCellsX will define the average dx size
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the C(Y) structure is defined by the function

        f(x,y,z) = -sin(x)sin(y)sin(z)+sin(2x)sin(y)+sin(2y)sin(z) +
                    sin(2z)sin(x)-cos(x)cos(y)cos(z)+sin(2x)cos(z) +
                    sin(2y)cos(x) + sin(2z)cos(y)
    */
    // Populate the mesh array
    info->numCellsX = numCellsX;
    info->numCellsY = numCellsX;
    info->numCellsZ = numCellsX;

    info->nElements = numCellsX * numCellsX * numCellsX;

    // define step dx

    float dx = (float)2.0 * PI / (float)numCellsX;

    // We will also calculate the porosity based on this discretization

    size_t count = 0;

    // Begin main loop

    for (long int i = 0; i < info->nElements; i++)
    {
        // break down i into
        int slice = i / (numCellsX * numCellsX);
        int row = (i - slice * numCellsX * numCellsX) / numCellsX;
        int col = i - slice * numCellsX * numCellsX - row * numCellsX;

        // Now get x,y,z numbers
        float x = -PI + dx / 2.0 + (float)col * dx;
        float y = -PI + dx / 2.0 + (float)row * dx;
        float z = -PI + dx / 2.0 + (float)slice * dx;

        // calculate f

        float f = -sin(x) * sin(y) * sin(z) + sin(2 * x) * sin(y) + sin(2 * y) * sin(z) +
                  sin(2 * z) * sin(x) - cos(x) * cos(y) * cos(z) + sin(2 * x) * cos(z) +
                  sin(2 * y) * cos(x) + sin(2 * z) * cos(y);

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->VF = (float)count / (float)info->nElements;

    return 0;
}

int SchwarzP_3D(char *P, float b, int numCellsX, meshInfo *info)
{
    /*
        Function SchwarzP_3D:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - int numCellsX will define the average dx size
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the Schwarz-P structure is defined by the function

        f(x,y,z) = cos x + cos y + cos z

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    info->numCellsX = numCellsX;
    info->numCellsY = numCellsX;
    info->numCellsZ = numCellsX;

    info->nElements = numCellsX * numCellsX * numCellsX;

    // Define dx

    float dx = (float)2.0 * PI / (float)numCellsX;

    // we will also calculate SVF based on this discretization

    size_t count = 0;

    // Begin main loop

    for (long int i = 0; i < info->nElements; i++)
    {
        // Break down i into the integer indexes
        int slice = i / (numCellsX * numCellsX);
        int row = (i - slice * numCellsX * numCellsX) / numCellsX;
        int col = i - slice * numCellsX * numCellsX - row * numCellsX;
        // Now get x,y,z numbers
        float x = -PI + dx / 2.0 + (float)col * dx;
        float y = -PI + dx / 2.0 + (float)row * dx;
        float z = -PI + dx / 2.0 + (float)slice * dx;

        float f = cos(x) + cos(y) + cos(z);

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->VF = (float)count / (float)info->nElements;

    return 0;
}

int genStruct(options *opts, meshInfo *info, char *P)
{
    /*
        Function getStruct:
        Inputs:
            - pointer to user options struct
            - pointer to mesh struct
            - pointer to char array where struct will be stored.
        Outputs:
            - none.

        The funtion will simply call other functions to generate the TPMS
        structure based on user input.
    */

    if (opts->TPMS_Type == 1)
    {
        if (opts->verbose)
            printf("Schwarz-P\n");
        SchwarzP_3D(P, opts->isoValues, info->numCellsX, info);
    }
    else if (opts->TPMS_Type == 7)
    {
        if (opts->verbose)
            printf("C(Y)\n");
        C_Y(P, opts->isoValues, info->numCellsX, info);
    }
    else
    {
        return 1;
    }
    return 0;
}

/*

    Surface Area Calculation:

*/

void SA(char *P, meshInfo *info, saveInfo *save)
{
    /*
        Fuction SA:
        Inputs:
            - pointer to P array holding the structure
            - pointer to meshInfo struc
        Outputs:
            - None.

        Calculates the surface area of the domain in units of voxel squared. The
        surface area is calculated by voxel-face-interface counting.
    */

    size_t interfaceCount = 0;

    // Copy struct variables locally
    int nRows, nCols, nSlices;

    nSlices = info->numCellsZ;
    nCols = info->numCellsX;
    nRows = info->numCellsY;

    // Loop over the whole structure
    int slice, row, col;

    for (long int i = 0; i < info->nElements; i++)
    {
        // if P is not solid, continue
        if (P[i] != 1)
            continue;
        // Break down i into the integer indexes
        slice = i / (nCols * nRows);
        row = (i - slice * nRows * nCols) / nCols;
        col = i - slice * nRows * nCols - row * nCols;

        // Test all cases

        // cols

        if (col != 0)
        {
            if (P[i] != P[i - 1])
                interfaceCount++;
        }

        if (col != nCols - 1)
        {
            if (P[i] != P[i + 1])
                interfaceCount++;
        }

        // rows

        if (row != 0)
        {
            if (P[i] != P[i - nCols])
                interfaceCount++;
        }

        if (row != nRows - 1)
        {
            if (P[i] != P[i + nCols])
                interfaceCount++;
        }

        // slices

        if (slice != 0)
        {
            if (P[i] != P[i - nRows * nCols])
                interfaceCount++;
        }

        if (slice != nSlices - 1)
        {
            if (P[i] != P[i + nRows * nCols])
                interfaceCount++;
        }
    }

    save->SA = (float)interfaceCount / (float)info->nElements;
    return;
}


/*

    Tortuosity related Functions

*/

void SetDC_Tau(float *DC, char *P, meshInfo* mesh)
{
    /*
        Function SetDC_Tau:
        
        Inputs:
            - pointer DC, array holding diffusion coefficients
            - pointer to P, array holding the structure.
            - pointer to struct holding the mesh info array.
        Outputs:
            - None.
        
        Function will populate the DC array with diffusion coefficients.
    */

    for (int i = 0; i < mesh->nElements; i++)
    {
        if (P[i] == 0)
            DC[i] = 1.0;
    }

    return;
}

int FloodFill3D(meshInfo *mesh, float *DC)
{
    /*
        FloodFill3D function:
        Inputs:
            - pointer to mesh struct
            - pointer to array with DC's
        Outputs:
            - None

        The function will search the domain, and will find non-participating media.
        The non-participating media will have the diffusion coefficient set to 0, and
        will receive the "wall" treatment.
    */

    char *Domain = (char *)malloc(mesh->nElements * sizeof(char));

    // Initialize all the impermeable matter in the domain:

    long int count = 0;

    for (long int index = 0; index < mesh->nElements; index++)
    {
        if (DC[index] == 0)
        {
            Domain[index] = 1;
        }
        else
        {
            Domain[index] = -1;
            count++;
        }
    }

    // Find Fluid in both boundaries, add to open list

    std::set<coord> cList;

    int left = 0;
    int right = mesh->numCellsX - 1;

    for (int row = 0; row < mesh->numCellsY; row++)
    {
        for (int slice = 0; slice < mesh->numCellsZ; slice++)
        {
            long int indexL = slice * mesh->numCellsX * mesh->numCellsY + row * mesh->numCellsX + left;
            long int indexR = slice * mesh->numCellsX * mesh->numCellsY + row * mesh->numCellsX + right;
            // set left
            if (Domain[indexL] == -1)
            {
                Domain[indexL] = 0;
                cList.insert(std::make_tuple(left, row, slice));
            }
            // set right
            if (Domain[indexR] == -1)
            {
                Domain[indexR] = 0;
                cList.insert(std::tuple(right, row, slice));
            }
        }
    }

    // Search Full Domain

    while (!cList.empty())
    {
        // pop first item on the list
        coord pop = *cList.begin();

        // remove from open list
        cList.erase(cList.begin());

        // get coordinates from the list
        int col = std::get<0>(pop);
        int row = std::get<1>(pop);
        int slice = std::get<2>(pop);

        /*
            We need to check North, South, East, West, Back, and Front for more fluid:

            North = col + 0, row - 1, slice + 0
            South = col + 0, row + 1, slice + 0
            East  = col + 1, row + 0, slice + 0
            West  = col - 1, row + 0, slice + 0
            Front = col + 0, row + 0, slice - 1
            Back  = col + 0, row + 0, slice + 1

            Note that diagonals are not considered a connection.
            This code assumes no periodic boundary conditions (currently).
        */

        int tempRow, tempCol, tempSlice;
        long int tempIndex;

        // North

        tempCol = col;
        tempSlice = slice;

        if (row != 0)
        {
            tempRow = row - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            tempRow = row + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // Front

        tempCol = col;
        tempRow = row;

        if (slice != 0)
        {
            tempSlice = slice - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            tempSlice = slice + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // West

        tempRow = row;
        tempSlice = slice;

        if (col != 0)
        {
            tempCol = col - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // East

        if (col != mesh->numCellsX - 1)
        {
            tempCol = col + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }
        // repeat until cList is empty
    }

    // Every flag that is still -1 means a non-participating media

    long int npCount = 0;

    for (int index = 0; index < mesh->nElements; index++)
    {
        if (Domain[index] != -1)
        {
            continue;
        }
        else
            DC[index] = 0;
        npCount++;
    }

    // memory management
    free(Domain);

    return 0;
}

int Disc3D_Tau(options *opts,
               meshInfo *mesh,
               float *DC,
               float *CoeffMatrix,
               float *RHS)
{
    /*
        Function Disc3D_Tau:
        Inputs:
            - pointer to options data structure
            - pointer to mesh data structure
            - pointer to float array DC holding diffusion coefficients
            - pointer to float array CoeffMatrix Coefficient Matrix
            - pointer to float array RHS holding right-hand side of discretized system.
        Output:
            - none.

        Function creates a discretization for a simulation of tortuosity. It will populate the
        Coefficient Matrix array and the RHS array (where BC's are held).
    */

    // Set necessary variables

    int nCols, nRows;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;

    float dx, dy, dz;
    dx = mesh->dx;
    dy = mesh->dy;
    dz = mesh->dz;

    int row, col, slice;
    float dw, de, ds, dn, db, df;

    for (long int i = 0; i < mesh->nElements; i++)
    {
        // dissolve index into rows and cols
        slice = i / (nRows * nCols);
        row = (i - slice * nRows * nCols) / nCols;
        col = (i - slice * nRows * nCols - row * nCols);

        // make sure CoeffMatrix and RHS are zero

        RHS[i] = 0;
        for (int k = 0; k < 7; k++)
        {
            CoeffMatrix[i * 7 + k] = 0;
        }

        /*
            Correct for non-participating media, analogous to
            pressure-decoupled solid velocity correction:
            https://doi.org/10.1016/j.ijheatmasstransfer.2009.12.057
        */

        if (DC[i] == 0)
        {
            // 1 * phi = 0;
            CoeffMatrix[i * 7 + 0] = 1;
            RHS[i] = 0;
            continue;
        }

        // Participating fluid

        /*
            Indexing for coeff marix:

            0 : P       i
            1 : W       i - 1
            2 : E       i + 1
            3 : S       i + nCols
            4 : N       i - nCols
            5 : B       i + nRows * nCols
            6 : F       i - nRows * nCols
        */

        // West

        if (col == 0)
        {
            // Left boundary
            dw = DC[i];
            RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
        }
        else if (DC[i - 1] != 0)
        {
            // West is participating media
            dw = DC[i];
            CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
        }

        // East

        if (col == mesh->numCellsX - 1)
        {
            // Right boundary
            de = DC[i];
            RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
        }
        else if (DC[i + 1] != 0)
        {
            // East is participating media
            de = DC[i];
            CoeffMatrix[i * 7 + 2] = de * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / dx;
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            if (DC[i + nCols] != 0)
            {
                // Participating South
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
        }

        // North

        if (row != 0)
        {
            if (DC[i - nCols] != 0)
            {
                // Participating North
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            if (DC[i + nCols * nRows] != 0)
            {
                // Participating Back
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        }

        // Front

        if (slice != 0)
        {
            if (DC[i - nCols * nRows] != 0)
            {
                // Participating Front
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
        }

    } // end for

    return 0;
}


int Disc3D_TauPB(options *opts,
               meshInfo *mesh,
               float *DC,
               float *CoeffMatrix,
               float *RHS)
{
    /*
        Function Disc3D_Tau:
        Inputs:
            - pointer to options data structure
            - pointer to mesh data structure
            - pointer to float array DC holding diffusion coefficients
            - pointer to float array CoeffMatrix Coefficient Matrix
            - pointer to float array RHS holding right-hand side of discretized system.
        Output:
            - none.

        Function creates a discretization for a simulation of tortuosity. It will populate the
        Coefficient Matrix array and the RHS array (where BC's are held). This version is
        for periodic BCs
    */

    // Set necessary variables

    int nCols, nRows, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    float dx, dy, dz;
    dx = mesh->dx;
    dy = mesh->dy;
    dz = mesh->dz;

    int row, col, slice;
    float dw, de, ds, dn, db, df;

    for (long int i = 0; i < mesh->nElements; i++)
    {
        // dissolve index into rows and cols
        slice = i / (nRows * nCols);
        row = (i - slice * nRows * nCols) / nCols;
        col = (i - slice * nRows * nCols - row * nCols);

        // make sure CoeffMatrix and RHS are zero

        RHS[i] = 0;
        for (int k = 0; k < 7; k++)
        {
            CoeffMatrix[i * 7 + k] = 0;
        }

        /*
            Correct for non-participating media, analogous to
            pressure-decoupled solid velocity correction:
            https://doi.org/10.1016/j.ijheatmasstransfer.2009.12.057
        */

        if (DC[i] == 0)
        {
            // 1 * phi = 0;
            CoeffMatrix[i * 7 + 0] = 1;
            RHS[i] = 0;
            continue;
        }

        // Participating fluid

        /*
            Indexing for coeff marix:

            0 : P       i
            1 : W       i - 1
            2 : E       i + 1
            3 : S       i + nCols
            4 : N       i - nCols
            5 : B       i + nRows * nCols
            6 : F       i - nRows * nCols
        */

        // West

        if (col == 0)
        {
            // Left boundary
            dw = DC[i];
            RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
        }
        else if (DC[i - 1] != 0)
        {
            // West is participating media
            dw = DC[i];
            CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
        }

        // East

        if (col == mesh->numCellsX - 1)
        {
            // Right boundary
            de = DC[i];
            RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
        }
        else if (DC[i + 1] != 0)
        {
            // East is participating media
            de = DC[i];
            CoeffMatrix[i * 7 + 2] = de * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / dx;
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            if (DC[i + nCols] != 0)
            {
                // Participating South
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
        } else
        {
            if (DC[slice * nCols * nRows + col] != 0)
            {
                // Periodioc BC
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / dy;
            }
        }

        // North

        if (row != 0)
        {
            if (DC[i - nCols] != 0)
            {
                // Participating North
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
        } else
        {
            if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] != 0)
            {
                // Periodic North
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            if (DC[i + nCols * nRows] != 0)
            {
                // Participating Back
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        } else
        {
            if (DC[row * nCols + col] != 0)
            {
                // Periodic Back
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        }

        // Front

        if (slice != 0)
        {
            if (DC[i - nCols * nRows] != 0)
            {
                // Participating Front
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
        }
        else
        {
            if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] != 0)
            {
                // Periodic Back
                db = DC[i];
                CoeffMatrix[i * 7 + 6] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        }
    }

    return 0;
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



int TauSim3D(options *opts, meshInfo *mesh, saveInfo *info, char *P)
{
    /*
        Function Tau3D_Sim:
        
        Inputs:
            - pointer to struct options
            - pointer to mesh struct
            - pointer to save struct
            - pointer to array holding the structure, P.
        Outputs:
            - None.
        
        Function will setup and run a tortuosity simulation based on the
        user entered options. All releant information is saved to struct.
    */
    
    mesh->dx = (float) 1.0 / mesh->numCellsX;
    mesh->dy = (float) 1.0 / mesh->numCellsY;
    mesh->dz = (float) 1.0 / mesh->numCellsZ; 

    // declare and define DC in the main flow channel

    float *DC = (float *)malloc(sizeof(float) * mesh->nElements);

    memset(DC, 0 , mesh->nElements * sizeof(float));

    // Populate the array based on the structure

    SetDC_Tau(DC, P, mesh);

    // Find participating media, this will remove cutoff channels

    printf("Flood Fill\n");

    FloodFill3D(mesh, DC);

    // allocate the arrays for simulation

    float *CoeffMatrix = (float *)malloc(mesh->nElements * 7 * sizeof(float));
    float *RHS = (float *)malloc(mesh->nElements * sizeof(float));
    float *Concentration = (float *)malloc(mesh->nElements * sizeof(float));

    // initialize memory

    memset(CoeffMatrix, 0, sizeof(float) * 7 * mesh->nElements);
    memset(RHS, 0, mesh->nElements * sizeof(float));
    memset(Concentration, 0, mesh->nElements * sizeof(float));

    // Linear initialize the memory

    for (long int i = 0; i < mesh->nElements; i++)
    {
        if (DC[i] == 0)
            continue;
        int slice = i / (mesh->numCellsX * mesh->numCellsY);
        int row = (i - slice * mesh->numCellsX * mesh->numCellsY)/mesh->numCellsX;
        int col = i - slice * mesh->numCellsX * mesh->numCellsY - row * mesh->numCellsX;
        Concentration[i] = ((float)col / mesh->numCellsX) * (opts->CRight - opts->CLeft) + opts->CLeft;
    }

    // Discretize

    printf("Discretize\n");

    if (opts->PB)
    {
        Disc3D_TauPB(opts, mesh, DC, CoeffMatrix, RHS);
    } else
    {
        Disc3D_Tau(opts, mesh, DC, CoeffMatrix, RHS);
    }

    

    // Solve

    int nDevices;
    
    cudaGetDeviceCount(&nDevices);

    if (nDevices < 1)
    {
        printf("No CUDA-capable GPU Detected! Exiting....\n");
        return 1;
    }

    // pointers to GPU arrays

    float *d_Coeff = NULL;
    float *d_RHS = NULL;
    float *d_Conc = NULL;
    float *d_ConcTemp = NULL;

    // initialize the GPU arrays

    initGPU_3DSOR(&d_Coeff, &d_RHS, &d_Conc, &d_ConcTemp, mesh);

    // Solve

    JI3D_SOR(CoeffMatrix, RHS, Concentration, d_Coeff,
                 d_RHS, d_Conc, d_ConcTemp, opts, mesh);

    // un-initialize

    unInitGPU_SOR(&d_Coeff, &d_RHS, &d_Conc, &d_ConcTemp);

    // Calculate Tortuosity

    double Q1 = 0;
    double Q2 = 0;
    int right = mesh->numCellsX - 1;
    int left = 0;

    for (int k = 0; k < mesh->numCellsZ; k++)
    {
        for (int i = 0; i < mesh->numCellsY; i++)
        {
            long int indexL = k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + left;
            long int indexR = k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + right;
            Q1 += DC[indexL] * (Concentration[indexL] - opts->CLeft) / (mesh->dx / 2);
            Q2 += DC[indexR] * (opts->CRight - Concentration[indexR]) / (mesh->dx / 2);
        }
    }

    double qAvg = (Q1 + Q2) / (2.0 * mesh->numCellsY * mesh->numCellsZ);

    info->Deff_TH_MAX = info->porosity * 1.0;
    info->Deff = qAvg / (opts->CRight - opts->CLeft);
    info->Tau = info->Deff_TH_MAX / info->Deff;

    if (opts->verbose == 1)
    {
        printf("VF = %1.3lf, DeffMax = %1.3e, Deff = %1.3e, Tau = %1.3e\n",
               info->porosity, info->Deff_TH_MAX, info->Deff, info->Tau);
    }

    // FILE *OUT;

    // OUT = fopen("Test_CY.csv", "w");
    // fprintf(OUT, "x,y,z,c\n");
    // for (int i = 0; i < mesh->numCellsY; i++)
    // {
    //     for (int j = 0; j < mesh->numCellsX; j++)
    //     {
    //         for (int k = 0; k < mesh->numCellsZ; k++)
    //         {
    //             fprintf(OUT, "%d,%d,%d,%1.3lf\n", j, i, k, Concentration[k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + j]);
    //         }
    //     }
    // }

    // fclose(OUT);

    return 0;
}

/*

    Size distributions related functions:

*/

int saveLabels3D(int *R,
                 meshInfo *structureInfo,
                 char *filename)
{
    /*
    Function saveLabels3D:
    Inputs:
    - pointer to R labels
    - pointer to sizeInfo struct
    - pointer to output filename
    Output:
    - None

    Function will save labels for this 3D simulation.
    */
    // read data structure
    int height, width;

    height = structureInfo->numCellsY;
    width = structureInfo->numCellsX;

    long int nElements = structureInfo->nElements;

    // Open File

    FILE *Particle;

    Particle = fopen(filename, "a+");

    fprintf(Particle, "x,y,z,R\n");

    int slice, row, col;

    // save everything

    for (int i = 0; i < nElements; i++)
    {
        if (R[i] != -1)
        {
            slice = i / (height * width);
            row = (i - slice * height * width) / width;
            col = (i - slice * height * width - row * width);
            fprintf(Particle, "%d,%d,%d,%d\n", col, row, slice, (int)R[i]);
        }
    }

    // close file

    fclose(Particle);

    return 0;
}

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