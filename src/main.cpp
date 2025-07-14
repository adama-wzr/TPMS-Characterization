/*

Main file will handle program execution.

Last modified 07/12/2025
Andre Adam
*/

#include <lib/TPMS_definitions.hpp>
#include <lib/data_structures.hpp>
#include <lib/usrInput.hpp>
#include <lib/TPMS_helpers.hpp>
#include <lib/surfaceArea.hpp>
#include <lib/TauSim.hpp>
#include <lib/sizeDistributions.hpp>
#include <lib/output.hpp>
#include <subDomainFF.hpp>
#include <math.h>


int main(int argc, char **argv)
{
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

    readInputGeneral(inputFilename, &opts);

    if (opts.verbose)
        printOptsGeneral(&opts);

    // Error Check Inputs

    /*
        Not yet implemented
    */

    // Create TPMS

    char *P;

    TPMS_Init(&P, &opts, &mesh);

    // Update the save array
    save.nVoxel = opts.nVoxels;
    save.porosity = mesh.porosity;
    save.SVF = mesh.SVF;

    // Calculate Surface Area, if applicable
    if (opts.runSA)
        SA(P, &mesh, &save);
    
    // Get N-Channels

    char *subDomains = (char *)malloc(sizeof(char) * mesh.nElements);

    subDomainFF(&mesh, P, subDomains);
    
    // calculate Tortuosity, if applicable

    // pore-space tortuosity
    if (opts.Tau)
        errorFlag = TauSim3D(&opts, &mesh, &save, P, 0);

    // solid space tortuosity
    if (opts.Tau)
        errorFlag = TauSim3D(&opts, &mesh, &save, P, 1);

    if(errorFlag)
        return 1;

    // Calculate size distributions

    omp_set_num_threads(opts.nThreads);

    if(opts.partSD)
        partSD_3D(&opts, &mesh, &save, P, 1);
    else
        save.part50 = 0;
    
    if (opts.poreSD)
        poreSD_3D(&opts, &mesh, &save, P, 0);
    else
        save.pore50 = 0;

    outputGeneral(&opts, &save);

    return 0;
}