/*

This file is dedicated to data structure definitions.

Last modified 03/27/2025
Andre Adam

*/

#ifndef _STRUCTS
#define _STRUCTS

#include <set>
#include <tuple>
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
    int runSA;              // calculate surface area flag
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

#endif