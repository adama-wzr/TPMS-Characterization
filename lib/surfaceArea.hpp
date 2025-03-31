/*

This file is dedicated to surface area functions.

Last modified 03/31/2025
Andre Adam

*/
#ifndef _SURF_AREA
#define _SURF_AREA

#include <data_structures.hpp>

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


#endif