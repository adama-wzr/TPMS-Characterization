
/*

TPMS Structure definitions


Last modified: 04/03/2025

Silven Stallard
*/


#ifndef _TPMS
#define _TPMS

#include <math.h>
#include <cfloat>
#include <data_structures.hpp>
#include <constants.hpp>

/*

    TPMS Definitions:

*/

// Index 1

int Gyroid(char *P, float b, meshInfo *info)
{
    /*
        Function Gyroid:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the Gyroid structure is defined by the function

        f(x,y,z) = cos(x)sin(y) + cos(y)sin(z) + cos(z)sin(x)

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

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

        float f = cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x);

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 2

int TPMS_D(char *P, float b, meshInfo *info)
{
    /*
        Function TPMS_D:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None
        
        Function will populate the P array based on the iso-bounds, b. In short,
        the D structure is defined by the function

        f(x,y,z) = sin(x)sin(y)sin(z) + sin(x)cos(y)cos(z) + 
                   cos(x)sin(y)cos(z) + cos(x)cos(y)sin(z)

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

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

        float f = sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + 
                  cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z);

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 3

int SchwarzP(char *P, float b, meshInfo *info)
{
    /*
        Function SchwarzP:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
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
    int numCellsX = info->numCellsX;

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

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 4

int IWP(char *P, float b, meshInfo *info)
{
    /*
        Function IWP:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None
        
        Function will populate the P array based on the iso-bounds, b. In short,
        the IWP structure is defined by the function

        f(x,y,z) = 2(cos(x)cos(y) + cos(y)cos(z) + cos(z)cos(x)) - (cos(2x) + cos(2y) + cos(2z))

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

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

        float f = 2 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) -
                  (cos(2 * x) + cos(2 * y) + cos(2 * z));

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 5

int Neovius(char *P, float b, meshInfo *info)
{
    /*
        Function Neovius:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None
        Function will populate the P array based on the iso-bounds, b. In short,
        the Neovius structure is defined by the function

        f(x,y,z) = 3(cos(x) + cos(y) + cos(z)) + 4cos(x)cos(y)cos(z)

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

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

        float f = 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z);

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 6

int C_Y(char *P, float b, meshInfo *info)
{
    /*
        Function C_Y:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the C(Y) structure is defined by the function

        f(x,y,z) = -sin(x)sin(y)sin(z)+sin(2x)sin(y)+sin(2y)sin(z) +
                    sin(2z)sin(x)-cos(x)cos(y)cos(z)+sin(2x)cos(z) +
                    sin(2y)cos(x) + sin(2z)cos(y)
    */
    int numCellsX = info->numCellsX;

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

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 7

int Lidinoid(char *P, float b, meshInfo *info)
{
    /*
        Function Lidinoid:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the Lidinoid structure is defined by the function

        f(x,y,z) = (sin(2x)cos(y)sin(z) + sin(x)sin(2y)cos(z) + cos(x)sin(y)sin(2z)) -
                  (cos(2x)cos(2y) + cos(2y)cos(2z) + cos(2z)cos(2x)) + 0.3

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

    // define step dx

    float dx = (float)2.0 * PI / (float)numCellsX;

    // We will also calculate the porosity based on this discretization

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

        float f = (sin(2 * x) * cos(y) * sin(z) + sin(x) * sin(2 * y) * cos(z) + cos(x) * sin(y) * sin(2 * z)) -
                  (cos(2 * x)*cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)) + 0.3;

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

// Index 8

int OCTO(char *P, float b, meshInfo *info)
{
    /*
        Function OCTO:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the OCTO structure is defined by the function

        f(x,y,z) = 0.6(cos(x)cos(y) + cos(y)cos(z) + cos(z)cos(x)) -
                   0.4(cos(x) + cos(y) + cos(z)) + 0.25

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

    // define step dx

    float dx = (float)2.0 * PI / (float)numCellsX;

    // We will also calculate the porosity based on this discretization

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

        float f = 0.6 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) -
                  0.4 * (cos(x) + cos(y) + cos(z)) + 0.25;

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}


// Index 9

int FRD(char *P, float b, meshInfo *info)
{
    /*
        Function FRD:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - pointer to meshInfo struct containing mesh parameters
        Outputs:
            - None

        Function will populate the P array based on the iso-bounds, b. In short,
        the FRD structure is defined by the function

        f(x,y,z) = 8cos(x)cos(y)cos(z) + cos(2x)cos(2y)cos(2z) - 
                  (cos(2x)cos(2y) + cos(2y)cos(2z) + cos(2z)cos(2x))

        It can be given a real thickness by defining the iso-values for which the structure
        is defined, i.e.

        P[x,y,z] = 1 if -b < f(x,y,z) < b
        P[x,y,z] = 0 otherwise

        This function will populate P according to the criteria above.
        P is the used as basis for the PSD analysis.
    */
    int numCellsX = info->numCellsX;

    // define step dx

    float dx = (float)2.0 * PI / (float)numCellsX;

    // We will also calculate the porosity based on this discretization

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

        float f = 8 * cos(x) * cos(y) * cos(z) + cos(2 * x) * cos(2 * y) * cos(2 * z) - 
                  (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x));

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->SVF = (float)count / (float)info->nElements;
    info->porosity = 1.0f - info->SVF;

    return 0;
}

/*
    TPMS Function Lookup Table:
*/

tpms_ptr TPMS_Functions[] = {
    Gyroid,
    TPMS_D,
    SchwarzP,
    IWP,
    Neovius,
    C_Y,
    Lidinoid,
    OCTO,
    FRD
};

#endif
