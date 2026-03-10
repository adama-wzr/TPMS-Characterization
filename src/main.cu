/*

Main file will handle program execution.

***

This is an exact copy of main.cpp

The reason for this is because I can't write cuda
code properly to compile with c++ nor do I know
how to do the separate compilation and link the
executables.

Therefore, this is the solution I found :)

***

Last modified 07/29/2025
Andre Adam
*/

#include <lib/TPMS_definitions.hpp>
#include <lib/data_structures.hpp>
#include <lib/usrInput.hpp>
#include <lib/TPMS_helpers.hpp>
#include <lib/surfaceArea.hpp>
#include <lib/TauSim.hpp>
#include <lib/sizeDistributions.hpp>
#include <lib/sfSim.hpp>
#include <lib/output.hpp>
#include <subDomainFF.hpp>
#include <math.h>



int main(int argc, char **argv)
{
    // performance
    fflush(stdout);
    // Declare data structure with general user options
    options opts;
    meshInfo mesh;
    saveInfo save;

    // Initialize save array

    InitSave(&save);

    // Other useful flags

    bool errorFlag = 0;

    // Read user input

    char inputFilename[40];

    sprintf(inputFilename, "input.txt");

    if (readInputGeneral(inputFilename, &opts) == 1)
    {
        // Return if no file exists
        printf("No Input File Found! Exiting....\n");
        return 1;
    }

    if (opts.verbose)
        printOptsGeneral(&opts);

    // Error Check Inputs

    if(errorCheckInput(&opts) == 1)
        return 1;

    // Create TPMS

    char *P;

    TPMS_Init(&P, &opts, &mesh);

    // Update the save array
    save.nVoxel = opts.nVoxels;
    save.porosity = mesh.porosity;
    save.SVF = mesh.SVF;
    
    // Get N-Channels

    char *subDomains = (char *)malloc(sizeof(char) * mesh.nElements);

    subDomainFF(&mesh, P, subDomains);

    // allocated space for sub-domain data

    mesh.sdInfo = (subDinfo *)malloc(sizeof(subDinfo) * mesh.nChannels);

    // Fully-Connected or not?

    subDomainFC(&mesh, subDomains);

    /*
    
        Surface Area (per channel?)

    */

    if (opts.runSA)
        SA(P, &mesh, &save, &opts);
    
    if (opts.runSA && opts.subOut)
        SA_sub(&opts, &mesh, subDomains);
        
    
    /*
    
        Tortuosity:
    
    */

    // pore-space tortuosity
    if (opts.Tau_f)
        errorFlag = TauSim3D(&opts, &mesh, &save, P, subDomains, 0);

    // solid space tortuosity
    if (opts.Tau_s)
        errorFlag = TauSim3D(&opts, &mesh, &save, P, subDomains, 1);

    if(errorFlag)
        return 1;
    
    // Shape Factors

    if (opts.runSF)
    {
        SF_Sim3D(&opts, &mesh, &save, P, subDomains);
    }

    /*
    
        Size - Distributions:
    
    */

    // Calculate size distributions

    if(opts.partSD)
        partSD_3D(&opts, &mesh, &save, P, 1);
    else
        save.part50 = 0;
    
    if (opts.poreSD)
        poreSD_3D(&opts, &mesh, &save, P, subDomains, 0);
    else
        save.pore50 = 0;

    // Shape Factors is last because it depends on others

    if (opts.runSF)
    {
        SF_Sim3D(&opts, &mesh, &save, P, subDomains);
    }

    outputGeneral(&opts, &save, &mesh);

    return 0;
}