/*

Shape Factor Simulation.

Last modified 03/07/2026
Silven Stallard
Andre Adam
*/

#ifndef _SF_SIM
#define _SF_SIM

#include <math.h>
#include <stdlib.h>
#include <data_structures.hpp>
#include <output.hpp>
#include <cpu_solvers/cpuSolvers.hpp>

#ifdef USE_CUDA
    #include <cuda_solvers/gpuSolve.cu>
#endif

#ifndef USE_CUDA
    #include <cpu_solvers/cpuErrorHandler.hpp>
#endif

void SetDC_Tau(float *DC, char *P, meshInfo* mesh, int POI)
{
    /*
        Function SetDC_Tau:
        
        Inputs:
            - pointer DC, array holding diffusion coefficients
            - pointer to P, array holding the structure.
            - pointer to struct holding the mesh info array.
            - integer POI (phase of interest)
        Outputs: 
            - None.
        
        Function will populate the DC array with diffusion coefficients.
    */

    for (int i = 0; i < mesh->nElements; i++)
    {
        if (P[i] == POI)
            DC[i] = 1.0;
    }

    return;
}

int FloodFill3D(meshInfo *mesh, float *DC)
{
    /*
        FloodFill3D function:
        Inputs:
            - pointer to mesh struct
            - pointer to array with DC's
        Outputs:
            - None

        The function will search the domain, and will find non-participating media.
        The non-participating media will have the diffusion coefficient set to 0, and
        will receive the "wall" treatment.
    */

    char *Domain = (char *)malloc(mesh->nElements * sizeof(char));

    // Initialize all the impermeable matter in the domain:

    long int count = 0;

    for (long int index = 0; index < mesh->nElements; index++)
    {
        if (DC[index] == 0)
        {
            Domain[index] = 1;
        }
        else
        {
            Domain[index] = -1;
            count++;
        }
    }

    // Find Fluid in both boundaries, add to open list

    std::set<coord> cList;

    int left = 0;
    int right = mesh->numCellsX - 1;

    for (int row = 0; row < mesh->numCellsY; row++)
    {
        for (int slice = 0; slice < mesh->numCellsZ; slice++)
        {
            long int indexL = slice * mesh->numCellsX * mesh->numCellsY + row * mesh->numCellsX + left;
            long int indexR = slice * mesh->numCellsX * mesh->numCellsY + row * mesh->numCellsX + right;
            // set left
            if (Domain[indexL] == -1)
            {
                Domain[indexL] = 0;
                cList.insert(std::make_tuple(left, row, slice));
            }
            // set right
            if (Domain[indexR] == -1)
            {
                Domain[indexR] = 0;
                cList.insert(std::tuple(right, row, slice));
            }
        }
    }

    // Search Full Domain

    while (!cList.empty())
    {
        // pop first item on the list
        coord pop = *cList.begin();

        // remove from open list
        cList.erase(cList.begin());

        // get coordinates from the list
        int col = std::get<0>(pop);
        int row = std::get<1>(pop);
        int slice = std::get<2>(pop);

        /*
            We need to check North, South, East, West, Back, and Front for more fluid:

            North = col + 0, row - 1, slice + 0
            South = col + 0, row + 1, slice + 0
            East  = col + 1, row + 0, slice + 0
            West  = col - 1, row + 0, slice + 0
            Front = col + 0, row + 0, slice - 1
            Back  = col + 0, row + 0, slice + 1

            Note that diagonals are not considered a connection.
            This code assumes no periodic boundary conditions (currently).
        */

        int tempRow, tempCol, tempSlice;
        long int tempIndex;

        // North

        tempCol = col;
        tempSlice = slice;

        if (row != 0)
        {
            tempRow = row - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            tempRow = row + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // Front

        tempCol = col;
        tempRow = row;

        if (slice != 0)
        {
            tempSlice = slice - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            tempSlice = slice + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // West

        tempRow = row;
        tempSlice = slice;

        if (col != 0)
        {
            tempCol = col - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // East

        if (col != mesh->numCellsX - 1)
        {
            tempCol = col + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }
        // repeat until cList is empty
    }

    // Every flag that is still -1 means a non-participating media

    long int npCount = 0;

    for (int index = 0; index < mesh->nElements; index++)
    {
        if (Domain[index] != -1)
        {
            continue;
        }
        else
            DC[index] = 0;
        npCount++;
    }

    // memory management
    free(Domain);

    return 0;
}

int FloodFill3D_PB(meshInfo *mesh, float *DC)
{
    /*
        FloodFill3D_PB function:
        Inputs:
            - pointer to mesh struct
            - pointer to array with DC's
        Outputs:
            - None

        The function will search the domain, and will find non-participating media.
        The non-participating media will have the diffusion coefficient set to 0, and
        will receive the "wall" treatment. This one includes periodic BCs
    */

    char *Domain = (char *)malloc(mesh->nElements * sizeof(char));

    // Initialize all the impermeable matter in the domain:

    long int count = 0;

    for (long int index = 0; index < mesh->nElements; index++)
    {
        if (DC[index] == 0)
        {
            Domain[index] = 1;
        }
        else
        {
            Domain[index] = -1;
            count++;
        }
    }

    // Find Fluid in both boundaries, add to open list

    std::set<coord> cList;

    int left = 0;
    int right = mesh->numCellsX - 1;

    for (int row = 0; row < mesh->numCellsY; row++)
    {
        for (int slice = 0; slice < mesh->numCellsZ; slice++)
        {
            long int indexL = slice * mesh->numCellsX * mesh->numCellsY + row * mesh->numCellsX + left;
            long int indexR = slice * mesh->numCellsX * mesh->numCellsY + row * mesh->numCellsX + right;
            // set left
            if (Domain[indexL] == -1)
            {
                Domain[indexL] = 0;
                cList.insert(std::make_tuple(left, row, slice));
            }
            // set right
            if (Domain[indexR] == -1)
            {
                Domain[indexR] = 0;
                cList.insert(std::tuple(right, row, slice));
            }
        }
    }

    // Search Full Domain

    while (!cList.empty())
    {
        // pop first item on the list
        coord pop = *cList.begin();

        // remove from open list
        cList.erase(cList.begin());

        // get coordinates from the list
        int col = std::get<0>(pop);
        int row = std::get<1>(pop);
        int slice = std::get<2>(pop);

        /*
            We need to check North, South, East, West, Back, and Front for more fluid:

            North = col + 0, row - 1, slice + 0
            South = col + 0, row + 1, slice + 0
            East  = col + 1, row + 0, slice + 0
            West  = col - 1, row + 0, slice + 0
            Front = col + 0, row + 0, slice - 1
            Back  = col + 0, row + 0, slice + 1

            Note that diagonals are not considered a connection.
            This code assumes periodic boundary conditions.
        */

        int tempRow, tempCol, tempSlice;
        long int tempIndex;

        // North

        tempCol = col;
        tempSlice = slice;

        if (row != 0)
        {
            tempRow = row - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        } else
        {
            // Periodic
            tempRow = mesh->numCellsY - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            tempRow = row + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        } else
        {
            // Periodic
            tempRow = 0;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // Front

        tempCol = col;
        tempRow = row;

        if (slice != 0)
        {
            tempSlice = slice - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        } else
        {
            // Periodic
            tempSlice = mesh->numCellsZ - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            tempSlice = slice + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }
        else
        {
            // Periodic
            tempSlice = 0;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // West

        tempRow = row;
        tempSlice = slice;

        if (col != 0)
        {
            tempCol = col - 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }

        // East

        if (col != mesh->numCellsX - 1)
        {
            tempCol = col + 1;
            tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            if (Domain[tempIndex] == -1)
            {
                Domain[tempIndex] = 0;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
        }
        // repeat until cList is empty
    }

    // Every flag that is still -1 means a non-participating media

    long int npCount = 0;

    for (int index = 0; index < mesh->nElements; index++)
    {
        if (Domain[index] != -1)
        {
            continue;
        }
        else
            DC[index] = 0;
        npCount++;
    }

    // memory management
    free(Domain);
    return 0;
}

int Disc3D_Tau(options *opts,
               meshInfo *mesh,
               float *DC,
               float *CoeffMatrix,
               float *RHS)
{
    /*
        Function Disc3D_Tau:
        Inputs:
            - pointer to options data structure
            - pointer to mesh data structure
            - pointer to float array DC holding diffusion coefficients
            - pointer to float array CoeffMatrix Coefficient Matrix
            - pointer to float array RHS holding right-hand side of discretized system.
        Output:
            - none.

        Function creates a discretization for a simulation of tortuosity. It will populate the
        Coefficient Matrix array and the RHS array (where BC's are held).
    */

    // Set necessary variables

    int nCols, nRows;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;

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

        if (DC[i] == 0)
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
            // Left boundary
            dw = DC[i];
            RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
        }
        else if (DC[i - 1] != 0)
        {
            // West is participating media
            dw = DC[i];
            CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
        }

        // East

        if (col == mesh->numCellsX - 1)
        {
            // Right boundary
            de = DC[i];
            RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
        }
        else if (DC[i + 1] != 0)
        {
            // East is participating media
            de = DC[i];
            CoeffMatrix[i * 7 + 2] = de * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / dx;
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            if (DC[i + nCols] != 0)
            {
                // Participating South
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
        }

        // North

        if (row != 0)
        {
            if (DC[i - nCols] != 0)
            {
                // Participating North
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            if (DC[i + nCols * nRows] != 0)
            {
                // Participating Back
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        }

        // Front

        if (slice != 0)
        {
            if (DC[i - nCols * nRows] != 0)
            {
                // Participating Front
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
        }

    } // end for

    return 0;
}

int Disc3D_TauPB(options *opts,
               meshInfo *mesh,
               float *DC,
               float *CoeffMatrix,
               float *RHS)
{
    /*
        Function Disc3D_Tau:
        Inputs:
            - pointer to options data structure
            - pointer to mesh data structure
            - pointer to float array DC holding diffusion coefficients
            - pointer to float array CoeffMatrix Coefficient Matrix
            - pointer to float array RHS holding right-hand side of discretized system.
        Output:
            - none.

        Function creates a discretization for a simulation of tortuosity. It will populate the
        Coefficient Matrix array and the RHS array (where BC's are held). This version is
        for periodic BCs
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

        if (DC[i] == 0)
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
            // Left boundary
            dw = DC[i];
            RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
        }
        else if (DC[i - 1] != 0)
        {
            // West is participating media
            dw = DC[i];
            CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
        }

        // East

        if (col == mesh->numCellsX - 1)
        {
            // Right boundary
            de = DC[i];
            RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
        }
        else if (DC[i + 1] != 0)
        {
            // East is participating media
            de = DC[i];
            CoeffMatrix[i * 7 + 2] = de * (dy * dz) / dx;
            CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / dx;
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            if (DC[i + nCols] != 0)
            {
                // Participating South
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
        } else
        {
            if (DC[slice * nCols * nRows + col] != 0)
            {
                // Periodioc BC
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / dy;
            }
        }

        // North

        if (row != 0)
        {
            if (DC[i - nCols] != 0)
            {
                // Participating North
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
        } else
        {
            if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] != 0)
            {
                // Periodic North
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            if (DC[i + nCols * nRows] != 0)
            {
                // Participating Back
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        } else
        {
            if (DC[row * nCols + col] != 0)
            {
                // Periodic Back
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        }

        // Front

        if (slice != 0)
        {
            if (DC[i - nCols * nRows] != 0)
            {
                // Participating Front
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
        }
        else
        {
            if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] != 0)
            {
                // Periodic Back
                db = DC[i];
                CoeffMatrix[i * 7 + 6] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
        }
    }

    return 0;
}

int TauSim3D(options *opts, meshInfo *mesh, saveInfo *save, char *P, char *subDomain, int POI)
{
    /*
        Function Tau3D_Sim:
        
        Inputs:
            - pointer to struct options
            - pointer to mesh struct
            - pointer to save struct
            - pointer to array holding the structure, P.
            - pointer to subDomain array
            - integer POI, phase of intereset
        Outputs:
            - None.
        
        Function will setup and run a tortuosity simulation based on the
        user entered options. All releant information is saved to struct.
    */
    
    mesh->dx = (float) 1.0 /mesh->numCellsX;
    mesh->dy = (float) 1.0 /mesh->numCellsY;
    mesh->dz = (float) 1.0 /mesh->numCellsZ;

    // declare and define DC in the main flow channel

    float *DC = (float *)malloc(sizeof(float) * mesh->nElements);

    memset(DC, 0 , mesh->nElements * sizeof(float));

    // Populate the array based on the structure

    SetDC_Tau(DC, P, mesh, POI);

    // Find participating media, this will remove cutoff channels

    if (opts->verbose)
        printf("Flood Fill\n");

    if (opts->PB)
        FloodFill3D_PB(mesh, DC);
    else
        FloodFill3D(mesh, DC);

    // allocate the arrays for simulation

    float *CoeffMatrix = (float *)malloc(mesh->nElements * 7 * sizeof(float));
    float *RHS = (float *)malloc(mesh->nElements * sizeof(float));
    float *Concentration = (float *)malloc(mesh->nElements * sizeof(float));

    // initialize memory

    memset(CoeffMatrix, 0, sizeof(float) * 7 * mesh->nElements);
    memset(RHS, 0, mesh->nElements * sizeof(float));
    memset(Concentration, 0, mesh->nElements * sizeof(float));

    // Linear initialize the concentration

    for (long int i = 0; i < mesh->nElements; i++)
    {
        if (DC[i] == 0)
            continue;
        int slice = i / (mesh->numCellsX * mesh->numCellsY);
        int row = (i - slice * mesh->numCellsX * mesh->numCellsY)/mesh->numCellsX;
        int col = i - slice * mesh->numCellsX * mesh->numCellsY - row * mesh->numCellsX;
        Concentration[i] = ((float)col / mesh->numCellsX) * (opts->CRight - opts->CLeft) + opts->CLeft;
    }

    // Discretize

    if(opts->verbose)
        printf("Discretizing\n");

    if (opts->PB)
    {
        Disc3D_TauPB(opts, mesh, DC, CoeffMatrix, RHS);
    } else
    {
        Disc3D_Tau(opts, mesh, DC, CoeffMatrix, RHS);
    }

    // Solve

    bool errorFlag = 0;

    if(opts->useGPU)
    {
        errorFlag = gpuHandler(opts, mesh, save, Concentration, CoeffMatrix, RHS);
        if(errorFlag)
            return 1;
    }
    else
    {
        // cpuSolve
        pGS3D_handle(opts, mesh, save, CoeffMatrix, RHS, Concentration);
    }

    // Calculate Tortuosity

    double Q1 = 0;
    double Q2 = 0;
    int right = mesh->numCellsX - 1;
    int left = 0;

    for (int k = 0; k < mesh->numCellsZ; k++)
    {
        for (int i = 0; i < mesh->numCellsY; i++)
        {
            long int indexL = k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + left;
            long int indexR = k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + right;
            Q1 += DC[indexL] * (Concentration[indexL] - opts->CLeft) / (mesh->dx / 2);
            Q2 += DC[indexR] * (opts->CRight - Concentration[indexR]) / (mesh->dx / 2);
        }
    }

    double qAvg = (Q1 + Q2) / (2.0 * mesh->numCellsY * mesh->numCellsZ);

    if (POI == 0)
    {
        // get effective porosity
        size_t count = 0;
        for(long int i = 0; i < mesh->nElements; i++)
        {
            if (DC[i] != 0)
                count++;
        }

        save->ePore = (float) count / mesh->nElements;

        // Calculate tortuosity

        save->Deff_TH_MAX = save->ePore;
        save->Deff = qAvg / (opts->CRight - opts->CLeft);
        save->Tau = save->Deff_TH_MAX / save->Deff;

        if (opts->verbose == 1)
        {
            printf("VF = %1.3lf, DeffMax = %1.3e, Deff = %1.3e, Tau = %1.3e\n",
                   save->porosity, save->Deff_TH_MAX, save->Deff, save->Tau);
        }
        // print CMAP if needed
        if(opts->CMAP)
        {
            char out_end[] = "_TauF.csv";
            char filename[200];
            strcpy(filename, opts->CMAP_Name);
            strncat(filename, out_end, 100);
            printCMAP(opts, mesh, filename, Concentration, P, 0);
        }
        if(opts->subOut)
        {
            for (int sub = 1; sub <= mesh->nChannels; sub++)
            {
                // skip if not fully connected
                if (mesh->sdInfo[sub - 1].FC == 0)
                    continue;

                // calculate local fluxes
                Q1 = 0;
                Q2 = 0;
                for (int k = 0; k < mesh->numCellsZ; k++)
                {
                    for (int i = 0; i < mesh->numCellsY; i++)
                    {
                        long int indexL = k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + left;
                        long int indexR = k * mesh->numCellsX * mesh->numCellsY + i * mesh->numCellsX + right;
                        if (subDomain[indexL] == sub)
                            Q1 += DC[indexL] * (Concentration[indexL] - opts->CLeft) / (mesh->dx / 2);
                        if (subDomain[indexR] == sub)
                            Q2 += DC[indexR] * (opts->CRight - Concentration[indexR]) / (mesh->dx / 2);
                    }
                }
                // calculate avg flux and tau
                qAvg = (Q1 + Q2) / (2.0 * mesh->numCellsY * mesh->numCellsZ);
                float D_TH_MAX = mesh->sdInfo[sub - 1].VF;
                float Deff = qAvg / (opts->CRight - opts->CLeft);
                mesh->sdInfo[sub - 1].Tau = D_TH_MAX / Deff;
                // print
                if (opts->verbose)
                    printf("sub = %d, VF = %1.3f, Tau = %1.3f\n", sub, mesh->sdInfo[sub - 1].VF, mesh->sdInfo[sub - 1].Tau);
                // print CMAP if needed
                if (opts->CMAP)
                {
                    char out_end[100];
                    sprintf(out_end, "TauSub%d.csv", sub);
                    char filename[200];
                    strcpy(filename, opts->CMAP_Name);
                    strncat(filename, out_end, 100);
                    printCMAP(opts, mesh, filename, Concentration, subDomain, sub);
                }
            }
        }
    }
    else if (POI == 1)
    {
        save->Deff_TH_MAX = save->SVF;
        save->Deff = qAvg / (opts->CRight - opts->CLeft);
        save->TauSolid = save->Deff_TH_MAX / save->Deff;

        if (opts->verbose == 1)
        {
            printf("VF = %1.3lf, DeffMax = %1.3e, Deff = %1.3e, Tau = %1.3e\n",
                   save->SVF, save->Deff_TH_MAX, save->Deff, save->TauSolid);
        }
        // print CMAP if needed
        if (opts->CMAP)
        {
            char out_end[] = "_TauS.csv";
            char filename[200];
            strcpy(filename, opts->CMAP_Name);
            strncat(filename, out_end, 100);
            printCMAP(opts, mesh, filename, Concentration, P, 1);
        }
    }

    return 0;
}

#endif
