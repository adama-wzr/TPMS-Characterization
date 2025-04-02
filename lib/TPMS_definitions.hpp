
/*

TPMS Structure definitions

Last modified: 04/02/2025

Silven Stallard
*/


#ifndef _TPMS
#define _TPMS

#include <math.h>
#include <cfloat>
#include <data_structures.hpp>
#include <constants.hpp>



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

int Gyroid(char *P, float b, int numCellsX, meshInfo *info)
{
    /*
        Function Gyroid:
        Inputs:
            - pointer to P, will store binary phase information
            - float b gives the iso-bounds for solid structure definition
            - int numCellsX will define the average dx size
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

        float f = cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x);

        if (f > -b && f < b)
        {
            P[i] = 1;
            count++;
        }
    }

    info->VF = (float)count / (float)info->nElements;

    return 0;
}

#endif
