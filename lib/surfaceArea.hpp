/*

This file is dedicated to surface area functions.

Last modified 03/31/2025
Andre Adam

*/
#ifndef _SURF_AREA
#define _SURF_AREA

#include <data_structures.hpp>

void SA(char *P, meshInfo *mesh, saveInfo *save, options *opts)
{
    /*
        Fuction SA:
        Inputs:
            - pointer to P array holding the structure
            - pointer to meshInfo struct
            - pointer to saveInfo struct
            - pointer to opts struct
        Outputs:
            - None.

        Calculates the surface area of the domain in units of voxel squared. The
        surface area is calculated by voxel-face-interface counting.
        Will calculate surface area with periodic boundaries (y and z) 
        if user enters the periodic boundary (pb) flag.
    */

    size_t interfaceCount = 0;

    // Copy struct variables locally
    int nRows, nCols, nSlices;

    nSlices = mesh->numCellsZ;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;

    // Loop over the whole structure
    int slice, row, col;

    for (long int i = 0; i < mesh->nElements; i++)
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
        else if(opts->PB)
        {
            // periodic (North) boundary
            if (P[i] != P[slice * nRows * nCols + (nRows - 1) * nCols + col])
                interfaceCount++;
        }

        if (row != nRows - 1)
        {
            if (P[i] != P[i + nCols])
                interfaceCount++;
        }
        else if(opts->PB)
        {
            // periodic (South) boundary
            if (P[i] != P[slice * nRows * nCols + col])
                interfaceCount++;
        }

        // slices

        if (slice != 0)
        {
            if (P[i] != P[i - nRows * nCols])
                interfaceCount++;
        }
        else if(opts->PB)
        {
            // periodic (Front) boundary
            if (P[i] != P[(nSlices - 1) * nRows * nCols + row * nCols + col])
                interfaceCount++;
        }

        if (slice != nSlices - 1)
        {
            if (P[i] != P[i + nRows * nCols])
                interfaceCount++;
        }
        else if(opts->PB)
        {
            // periodic (Back) boundary
            if (P[i] != P[row * nCols + col])
                interfaceCount++;
        }
    }

    save->SA = (float)interfaceCount / (float)mesh->nElements;

    return;
}

void SA_sub(options *opts, meshInfo *mesh, char *subD)
{
    /*
        Function SA_sub:
        Inputs:
            - pointer to options struct
            - pointer to meshInfo struct
            - pointer to sub-domains array
        Output:
            - None
        
        Function will calculate the specific surface area
        on a sub-domain basis.
    */

    size_t interfaceCount = 0;

    // Copy struct variables locally
    int nRows, nCols, nSlices;

    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    // Loop over the whole structure
    int slice, row, col;

    for(int nSub = 1; nSub <= mesh->nChannels; nSub++)
    {
        // if not fully-connected, not interested
        if (mesh->sdInfo[nSub-1].FC == 0)
            continue;

        // zero counts
        interfaceCount = 0;
        
        // scan full domain
        for(int i = 0; i < mesh->nElements; i++)
        {
            // continue if not correct phase
            if (subD[i] != nSub)
                continue;

            // Break down i into the integer indexes
            slice = i / (nCols * nRows);
            row = (i - slice * nRows * nCols) / nCols;
            col = i - slice * nRows * nCols - row * nCols;

            // search immediate neighbors
            
            // cols
            if (col != 0)
            {
                if (subD[i] != subD[i - 1])
                    interfaceCount++;
            }

            if (col != nCols - 1)
            {
                if (subD[i] != subD[i + 1])
                    interfaceCount++;
            }
            // rows

            if (row != 0)
            {
                if (subD[i] != subD[i - nCols])
                    interfaceCount++;
            }
            else if(opts->PB)
            {
                // periodic (North) boundary
                if (subD[i] != subD[slice * nRows * nCols + (nRows - 1) * nCols + col])
                    interfaceCount++;
            }

            if (row != nRows - 1)
            {
                if (subD[i] != subD[i + nCols])
                    interfaceCount++;
            }
            else if(opts->PB)
            {
                // periodic (South) boundary
                if (subD[i] != subD[slice * nRows * nCols + col])
                    interfaceCount++;
            }

            // slices

            if (slice != 0)
            {
                if (subD[i] != subD[i - nRows * nCols])
                    interfaceCount++;
            }
            else if(opts->PB)
            {
                // periodic (Front) boundary
                if (subD[i] != subD[(nSlices - 1) * nRows * nCols + row * nCols + col])
                    interfaceCount++;
            }

            if (slice != nSlices - 1)
            {
                if (subD[i] != subD[i + nRows * nCols])
                    interfaceCount++;
            }
            else if(opts->PB)
            {
                // periodic (Back) boundary
                if (subD[i] != subD[row * nCols + col])
                    interfaceCount++;
            }
        }
        // done searching
        mesh->sdInfo[nSub - 1].SA = (float) interfaceCount / mesh->sdInfo[nSub - 1].nElements;
        if(opts->verbose)
            printf("Channel = %d, SSA = %1.3e\n", nSub - 1, mesh->sdInfo[nSub - 1].SA );
    }
    return;
}


#endif