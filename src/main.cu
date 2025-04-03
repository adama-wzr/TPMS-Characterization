/*

Main file will handle program execution.

Last modified 04/03/2025
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
#include <math.h>


int main(int argc, char **argv)
{
    // Declare data structure with general user options
    options opts;
    meshInfo mesh;
    saveInfo save;

    // Other useful flags

    bool errorFlag = 0;

    // Read user input

    char inputFilename[40];

    sprintf(inputFilename, "input.txt");

    readInputGeneral(inputFilename, &opts);

    if (opts.verbose)
        printOptsGeneral(&opts);

    // Error Check Inputs



    // Create TPMS

    char *P;

    TPMS_Init(&P, &opts, &mesh);

    // Calculate Surface Area, if applicable
    if (opts.runSA)
        SA(P, &mesh, &save);
    
    
    // calculate Tortuosity, if applicable

    if (opts.Tau)
        errorFlag = TauSim3D(&opts, &mesh, &save, P);

    
    if(errorFlag)
        return 1;

    // Calculate size distributions

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