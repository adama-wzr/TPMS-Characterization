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

/*

Legacy/Debug Functions

*/

void saveDC_SF(float *DC, meshInfo* mesh)
{
    /*
        Function saveDC_SF:

        This is a function to export the information on the discretization
        for the shape factor code.

        I was originally created as a debug option, and I will just leave it
        here if anyone wants to use.
    */

    FILE *TEST = fopen("PhaseInfo_SF.csv", "w+");

    fprintf(TEST, "x,y,z,P\n");

    int row, col, slice;

    for(int i = 0; i < mesh->nElements; i++)
    {
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX * mesh->numCellsY) / mesh->numCellsX;
        col = (i - slice * mesh->numCellsY * mesh->numCellsX - row * mesh->numCellsX);
        fprintf(TEST, "%d,%d,%d,%1.0f\n", col, row, slice, DC[i]);
    }

    fclose(TEST);

    return;
}

void saveTemp_SF(float *Concentration, meshInfo* mesh)
{
    /*
        Function saveTemp_SF:

        This is a function to save temperatures in the shape factor simulation.

        Function not called anywhere, just here for debugging.
    */

    // save to see temperatures calculated

    FILE *TEST = fopen("TempInfo_SF.csv", "w+");

    fprintf(TEST, "x,y,z,T\n");

    int row, col, slice;

    for(int i = 0; i < mesh->nElements; i++)
    {
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX * mesh->numCellsY) / mesh->numCellsX;
        col = (i - slice * mesh->numCellsY * mesh->numCellsX - row * mesh->numCellsX);
        fprintf(TEST, "%d,%d,%d,%1.3f\n", col, row, slice, Concentration[i]);
    }

    fclose(TEST);

    return;
}

/*

    Discretization Functions

*/

void SetDC_SF(float *DC, char *subDomain, meshInfo* mesh)
{
    /*
        Function SetDC_Tau:
        
        Inputs:
            - pointer DC, array holding diffusion coefficients
            - pointer to subDomain, array holding the structure's subdomain information.
            - pointer to struct holding the mesh info array.
            - integer POI (phase of interest)
        Outputs: 
            - None.
        
        Function will populate the DC array with correct index for Shape Factor simulation.

        1. Participating media
        2. First channel BC (T=T_high)
        3. Second Channel BC (T=T_low)
    */

    for (int i = 0; i < mesh->nElements; i++)
    {
        if (subDomain[i] == 0)
            DC[i] = 1;
        else if(subDomain[i] == 1)
            DC[i] = 2;
        else if(subDomain[i] == 2)
            DC[i] = 3;
    }

    return;
}

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

/*

    Simulation Control Function

*/

int SF_Sim3D(options *opts, meshInfo *mesh, saveInfo *save, char *P, char *subDomain)
{
    /*
        Function SF_Sim3D:
        
        Inputs:
            - pointer to struct options
            - pointer to mesh struct
            - pointer to save struct
            - pointer to array holding the structure, P.
            - pointer to subDomain array
        Outputs:
            - None.
        
        Function will setup and run a tortuosity simulation based on the
        user entered options. All releant information is saved to struct.
    */
    
    mesh->dx = (float) 1.0 /mesh->numCellsX;
    mesh->dy = (float) 1.0 /mesh->numCellsY;
    mesh->dz = (float) 1.0 /mesh->numCellsZ;

    int nRows, nCols, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    // declare and define DC in the main flow channel

    float *DC = (float *)malloc(sizeof(float) * mesh->nElements);

    memset(DC, 0 , mesh->nElements * sizeof(float));

    // Sub-Domain Info is absolutely necessary for this simulation

    if(mesh->nChannels != 2)
    {
        printf("Error Detected: nChannels = %d, nFC = %d\n", mesh->nChannels, mesh->nFC);
        printf("Returning.....");
        return 1;
        /*
        
        NEED TO ALSO CHECK IF THEY ARE FULLY CONNECTED!
        CURRENTLY NOT CHECKED!!!!!!

        */
    }

    // Populate the array based on the structure

    SetDC_SF(DC, subDomain, mesh);

    // allocate the arrays for simulation

    float *CoeffMatrix = (float *)malloc(mesh->nElements * 7 * sizeof(float));
    float *RHS = (float *)malloc(mesh->nElements * sizeof(float));
    float *Concentration = (float *)malloc(mesh->nElements * sizeof(float));

    // initialize memory

    memset(CoeffMatrix, 0, sizeof(float) * 7 * mesh->nElements);
    memset(RHS, 0, mesh->nElements * sizeof(float));

    // Initializing to CLeft, still converges fast
    memset(Concentration, 0, mesh->nElements * sizeof(float));
    for(int i = 0; i < mesh->nElements; i++)
    {
        if(DC[i] == 1)
            Concentration[i] = opts->CLeft;
    }

    // Discretize

    if(opts->verbose)
        printf("Discretizing\n");

    Disc3D_SF_PB(opts, mesh, DC, CoeffMatrix, RHS);

    // Solve

    bool errorFlag = 0;

    /*
        Code runs fast on CPU, I won't implement GPU support for now.
    */

    if(opts->useGPU)
    {
    //     errorFlag = gpuHandler(opts, mesh, save, Concentration, CoeffMatrix, RHS);
    //     if(errorFlag)
    //         return 1;
    // }
    // else
    // {
    //     // cpuSolve
    //     pGS3D_handle(opts, mesh, save, CoeffMatrix, RHS, Concentration);
        printf("GPU not currently Supported. Using CPU instead.\n");
    }

    // CPU solve
    pGS3D_SF_handle(opts, mesh, save, CoeffMatrix, RHS, Concentration);

    // Calculate SF

    // 1. Participating media
    // 2. First channel BC (T=T_high)
    // 3. Second Channel BC (T=T_low)

    float Q_21 = 0.0;
    float Q_13 = 0.0;

    int slice, row, col;

    for (long int i = 0; i < mesh->nElements; i++)
    {
        // if not participating media, ignore
        if (DC[i] != 1)
            continue;

        // get index components
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX * mesh->numCellsY) / mesh->numCellsX;
        col = (i - slice * mesh->numCellsY * mesh->numCellsX - row * mesh->numCellsX);

        // check all faces for neighbors

        if (col == 0)
        {
            // Periodic West
            if (DC[i + nCols - 1] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i + nCols - 1] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }

            // East
            if (DC[i + 1] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i + 1] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
        }
        else if (col == mesh->numCellsX - 1)
        {
            // West
            if (DC[i - 1] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i - 1] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }

            // Periodic East
            if (DC[i - (mesh->numCellsX - 1)] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i - (mesh->numCellsX - 1)] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
        }
        else
        {
            // West
            if (DC[i - 1] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i - 1] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }

            // East
            if (DC[i + 1] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i + 1] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
        }

        if (row == 0)
        {
            // Periodic North
            if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }

            // South
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
        }
        else if (row == nRows - 1)
        {
            // North
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }

            // Periodic South
            if (DC[slice * nCols * nRows + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
        }
        else
        {
            // North
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }

            // South
            if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
        }

        if (slice == 0)
        {
            // Periodic Front
            if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }

            // Back
            if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
        }
        else if (slice == nSlices - 1)
        {
            // Front
            if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }

            // Periodic Back
            if (DC[row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[row * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
        }
        else
        {
            // Front
            if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }

            // Back
            if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Concentration[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Concentration[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
        }
    } // end of the for loop

    // Calculate the shape factor

    float S = Q_21/(opts->CLeft - opts->CRight);

    printf("S = %1.3e, Q_21 = %1.3e, Q_13 = %1.3e, Residual = %1.3e\n", S, Q_21, Q_13, fabs(Q_21 - Q_13));
    

    return 0;
}

#endif
