/*

This file is dedicated to data structure definitions.

Last modified 07/18/2025
Andre Adam

*/

#ifndef _STRUCTS
#define _STRUCTS

#include <set>
#include <tuple>

// Sub-Domain Information

typedef struct
{
    bool FC;                // 0 for no-flow, 1 for flow
    long int nElements;     // total pixels in volume
    float VF;               // local volume fraction
    float SA;               // surface area
    float Tau;              // tortuosity
    float pore50;           // avg pore size
}subDinfo;

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
    int subOut;             // Per-sub domain output flag
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
    float SVF; // solid volume fraction
    float porosity;
    int nChannels;
    bool congruent;     // 0 == false, 1 == true (currently unused)
    int nFC;            // number of 'fully-connected' subdomains
    subDinfo *sdInfo;   // pointer to array of structs of sub-domain info
} meshInfo;

// Array holding output data

typedef struct
{
    float porosity;
    float SVF;
    float Deff_TH_MAX;
    float Deff;
    float Tau;
    float TauSolid;
    float SA;
    long int nElements;
    float pore50;
    float part50;
    int nVoxel;
} saveInfo;

// Define coords for Flood Fill

typedef std::tuple<int, int, int> coord;

// Declare the type of function pointers for TPMS definitions

typedef int (*tpms_ptr)(char*, float, meshInfo*);


/*

Initializers:

*/

void InitSave(saveInfo *save)
{
    /*
        This function simply initializes a declared save function.
    */

    save->porosity = 0.0f;
    save->SVF = 0.0f;
    save->Tau = 0.0f;
    save->TauSolid = 0.0f;
    save->SA = 0.0f;
    save->Deff = 0.0f;
    save->Deff_TH_MAX = 0.0f;
    save->nElements = 1;
    save->pore50 = 0.0f;
    save->pore50 = 0.0f;
    save->nVoxel = 1;

    
    return;
}

#endif