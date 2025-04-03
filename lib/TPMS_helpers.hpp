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
    Lookup table for TPMS Pinch:
*/

static const float TPMS_Pinch[] = {
    1.414,
    0.999,
    1.0,
    3.0,
    0.749,
    1.835,
    0.299,
    0.055,
    1.65,
    0.794,
    0.265,
    0.599,
    0.0,
    0.158,
    0.319,
    2.999,
    0.16,
    0.299,
    1.407,
    0.451,
    0.378,
    1.653,
    0.999,
    0.0,
    2.249,
    3.419,
    2.999
};

/*
    Lookup table for TPMS Crit:
*/

static const float TPMS_Crit[] = {
    1.5,
    1.414,
    3.0,
    4.998,
    13.0,
    3.535,
    2.7,
    3.249,
    10.0,
    1.414,
    1.75,
    1.8,
    1.0,
    4.25,
    1.82,
    7.5,
    2.41,
    1.6,
    6.455,
    4.95,
    2.828,
    2.0,
    3.0,
    4.0,
    4.0,
    10.5,
    12
};

/*

    TPMS Definitions:

*/

int TPMS_Init(char **P, options *opts, meshInfo *mesh)
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

    *P = (char *)malloc(sizeof(char) * mesh->nElements);
    
    memset(*P, 0, sizeof(char) * mesh->nElements);

    // Call function to generate TPMS according to user input

    TPMS_Functions[opts->TPMS_Type - 1](*P, opts->isoValues, mesh);
    
    return 0;
}

#endif