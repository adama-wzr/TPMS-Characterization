/*

User Input-parsing

Last modified 04/02/2025
Andre Adam

*/

#ifndef _TPMS_HELPER
#define _TPMS_HELPER

#include <stdlib.h>
#include <math.h>

#include <lib/constants.hpp>
#include <lib/data_structures.hpp>
#include <lib/TPMS_definitions.hpp>

/*
    Lookup table for TPMS Names:
*/

static const char *TPMS_Names[] = {
    "Gyroid", "Diamond", "Primitive",
    "IWP", "Neovius", "C(Y)",
    "Lidinoid", "OCTO", "FRD",
    "S", "P+C(P)", "Split-P",
    "F", "C(D)", "G'",
    "G'2", "D'", "K",
    "C(S)", "Y", "+-Y",
    "C(+-Y)", "C(I2-Y)", "W",
    "Q*", "C(G)", "Slotted-P"
};

/*

    TPMS Definitions:

*/

int TPMS_Init(char *P, options *opts, meshInfo *mesh)
{
    /*
        TPMS_Init Function:
        Inputs:
            - pointer to P (empty)
            - pointer to opts struct
            - pointer to mesh struct
        Outputs:
            - None
        Function populates mesh struct with user options. Then allocates space
        and creates the user defined TPMS structure on array P.
    */
    // Populate mesh
    mesh->numCellsX = opts->nVoxels;
    mesh->numCellsY = opts->nVoxels;
    mesh->numCellsZ = opts->nVoxels;
    mesh->nElements = (long int)pow(opts->nVoxels, 3);

    float dx = (float)2.0 * PI / (float)mesh->numCellsX;
    float dy = (float)2.0 * PI / (float)mesh->numCellsX;
    float dz = (float)2.0 * PI / (float)mesh->numCellsX;

    // Allocate and initialize P

    P = (char *)malloc(sizeof(char) * mesh->nElements);
    
    memset(P, 0, sizeof(char) * mesh->nElements);

    // Call function to generate TPMS according to user input

    TPMS_Functions[opts->TPMS_Type - 1](P, opts->isoValues, mesh);
    
    return 0;
}

#endif