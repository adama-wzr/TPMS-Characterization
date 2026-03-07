/*
Discretization function for Shape Factor simulation for TPMS.

3/5/2026
Last modified by Silven Stallard

Once this is working, I will merge this function on the sfSim.hpp file.
I will assume this is working for now.

- Andre
*/

#ifndef _DISC_SF

#define _DISC_SF

#include <data_structures.hpp>

int Disc3D_SF_PB(options *opts,
               meshInfo *mesh,
               float *DC,
               float *CoeffMatrix,
               float *RHS)
{
    /*
        Function Disc3D_SF_PB:
        Inputs:
            - pointer to options data structure
            - pointer to mesh data structure
            - pointer to float array DC holding diffusion coefficients
            - pointer to float array CoeffMatrix Coefficient Matrix
            - pointer to float array RHS holding right-hand side of discretized system.
        Output:
            - none.

        Function creates a discretization for a simulation of shape factor. It will populate the
        Coefficient Matrix array and the RHS array (where BC's are held). This version is
        for periodic BCs. Each element boundary has three BC options: 

        1. Participating media
        2. First channel BC (T=T_high)
        3. Second Channel BC (T=T_low)
    */

    // Set necessary variables

    int nCols, nRows, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    float dx, dy, dz;
    dx = mesh->dx;
    dy = mesh->dy;
    dz = mesh->dz;

    int row, col, slice;
    float dw, de, ds, dn, db, df;

    for (long int i = 0; i < mesh->nElements; i++)
    {
        // dissolve index into rows and cols
        slice = i / (nRows * nCols);
        row = (i - slice * nRows * nCols) / nCols;
        col = (i - slice * nRows * nCols - row * nCols);

        // make sure CoeffMatrix and RHS are zero

        RHS[i] = 0;
        for (int k = 0; k < 7; k++)
        {
            CoeffMatrix[i * 7 + k] = 0;
        }

        /*
            Correct for non-participating media, analogous to
            pressure-decoupled solid velocity correction:
            https://doi.org/10.1016/j.ijheatmasstransfer.2009.12.057
        */

        if (DC[i] != 1)
        {
            // 1 * phi = 0;
            CoeffMatrix[i * 7 + 0] = 1;
            RHS[i] = 0;
            continue;
        }

        // Participating fluid

        /*
            Indexing for coeff marix:

            0 : P       i
            1 : W       i - 1
            2 : E       i + 1
            3 : S       i + nCols
            4 : N       i - nCols
            5 : B       i + nRows * nCols
            6 : F       i - nRows * nCols
        */

        // West

        if (col == 0)
        {
            // Periodic Boundary
            if (DC[i - 1 + nCols] == 1) 
            {
            // Left boundary, participating media
            dw = DC[i];
            CoeffMatrix[i * 7 + 2] = dw * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
            }
            else if (DC[i - 1 + nCols] == 2)
            {
            // Left boundary, first channel
            dw = DC[i];
            RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
            else if (DC[i - 1 + nCols] == 3)
            {
            // Left boundary, second channel
            dw = DC[i];
            RHS[i] -= opts->CRight * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
        } else  
        {
            // Non-Boundary Neighbor
            if (DC[i - 1] == 1)
            {
            // West is participating media
            dw = DC[i];
            CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
            }
            else if (DC[i - 1] == 2)
            {
            // West is first channel
            dw = DC[i];
            RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
            else if (DC[i - 1] == 3)
            {
            // West is second channel
            dw = DC[i];
            RHS[i] -= opts->CRight * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
        }

        // East

        if (col == mesh->numCellsX - 1)
        {
            // Periodic Boundary
            if (DC[i + 1 - nCols] == 1)
            {
            // Right boundary, participating media
            de = DC[i];
            CoeffMatrix[i * 7 + 2] = de * (dy * dz) / (dx);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx);
            }
            else if (DC[i + 1 - nCols] == 2)
            {
            // Right boundary, first channel
            de = DC[i];
            RHS[i] -= opts->CLeft * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
            if (DC[i + 1 - nCols] == 3)
            {
            // Right boundary, second channel
            de = DC[i];
            RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
        } else 
        {
            // Non-Boundary Neighbor
            if (DC[i + 1] == 1)
            {
            // East, participating media
            de = DC[i];
            CoeffMatrix[i * 7 + 2] = de * (dy * dz) / (dx);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx);
            }
            else if (DC[i + 1] == 2)
            {
            // East, first channel
            de = DC[i];
            RHS[i] -= opts->CLeft * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
            if (DC[i + 1] == 3)
            {
            // East, second channel
            de = DC[i];
            RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            // Non-Boundary Neighbor
            if (DC[i + nCols] == 1)
            {
                // South, participating
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
            else if (DC[i + nCols] == 2)
            {
                // South, first channel
                ds = DC[i];
                RHS[i] -= opts->CLeft * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / (dy / 2);
            }
            else if (DC[i + nCols] == 3)
            {
                // South, second channel
                ds = DC[i];
                RHS[i] -= opts->CRight * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / (dy / 2);
            }
        } else
        {
            // Periodic Boundary
            if (DC[slice * nCols * nRows + col] == 1) //Periodic Boundary
            {
                // South, participating
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / dy;
            }
            else if (DC[slice * nCols * nRows + col] == 2)
            {
                // South, first channel
                ds = DC[i];
                RHS[i] -= opts->CLeft * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / (dy / 2);
            }
            else if (DC[slice * nCols * nRows + col] == 3)
            {
                // South, second channel
                ds = DC[i];
                RHS[i] -= opts->CRight * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / (dy / 2);
            }
        }

        // North

        if (row != 0)
        {
            // Non-Boundary Neighbor
            if (DC[i - nCols] == 1)
            {
                // North, participating media
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
            else if (DC[i - nCols] == 2)
            {
                // North, first channel
                dn = DC[i];
                RHS[i] -= opts->CLeft * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            }
            else if (DC[i - nCols] == 3)
            {
                // North, second channel
                dn = DC[i];
                RHS[i] -= opts->CRight * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            }
        } else
        {
            // Periodic Boundary
            if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 1) 
            {
                //North, participating media
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            } else if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 2)
            {
                // North, first channel
                dn = DC[i];
                RHS[i] -= opts->CLeft * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            } else if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 3)
            {
                // North, second channel
                dn = DC[i];
                RHS[i] -= opts->CRight * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            // Non-Boundary Neighbor
            if (DC[i + nCols * nRows] == 1)
            {
                // Back, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (DC[i + nCols * nRows] == 2)
            {
                // Back, first channel
                db = DC[i];
                RHS[i] -= opts->CLeft * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
            else if (DC[i + nCols * nRows] == 3)
            {
                // Back, second channel
                db = DC[i];
                RHS[i] -= opts->CRight * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
        } else
        {
            //Periodic Boundary
            if (DC[row * nCols + col] == 1)
            {
                // Back, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (DC[row * nCols + col] == 2)
            {
                // Periodic Back
                db = DC[i];
                RHS[i] -= opts->CLeft * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
            else if (DC[row * nCols + col] == 3)
            {
                // Periodic Back
                db = DC[i];
                RHS[i] -= opts->CRight * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
        }

        // Front

        if (slice != 0)
        {
            // Non-Boundary Neighbor
            if (DC[i - nCols * nRows] == 1)
            {
                // Front, participating media
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
            else if (DC[i - nCols * nRows] == 2)
            {
                // Front, first channel
                df = DC[i];
                RHS[i] -= opts->CLeft * df * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (dz / 2);
            }
            else if (DC[i - nCols * nRows] == 3)
            {
                // Front, second channel
                df = DC[i];
                RHS[i] -= opts->CRight * df * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (dz / 2);
            }
        }
        else
        {
            // Periodic Boundary
            if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 1)
            {
                // Front, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 6] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                // Front, first channel
                df = DC[i];
                RHS[i] -= opts->CLeft * df * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (dz / 2);
            }
            else if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                // Front, second channel
                df = DC[i];
                RHS[i] -= opts->CRight * df * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (dz / 2);
            }
        }
    }

    return 0;
}

#endif
