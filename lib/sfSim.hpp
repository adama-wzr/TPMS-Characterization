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
#include <cfloat>
#include <data_structures.hpp>
#include <constants.hpp>
#include <output.hpp>
#include <cpu_solvers/cpuSolvers.hpp>
#include <surfaceArea.hpp>
#include <marchingCubes.hpp>
#include <sizeDistributions.hpp>
#include <omp.h>
#include "TPMS_definitions.hpp"

#ifdef USE_CUDA
    #include <cuda_solvers/gpuSolve.cu>
#endif

#ifndef USE_CUDA
    #include <cpu_solvers/cpuErrorHandler.hpp>
#endif

/*

Legacy/Debug Functions

*/

void fixSD_info(options *opts, meshInfo *mesh, char *subDomain)
{
    /*
     * Add function description later if it works.
     *
     *
     * */
    // memset(subDomain, 0, sizeof(char)*mesh->nElements);

    float dx = mesh->dx;

    int numCellsX = mesh->numCellsX;

    for(long int i = 0; i < mesh->nElements; i++)
    {
        // Break down i into the integer indexes
        int slice = i / (numCellsX * numCellsX);
        int row = (i - slice * numCellsX * numCellsX) / numCellsX;
        int col = i - slice * numCellsX * numCellsX - row * numCellsX;
        // Now get x,y,z numbers
        float x = -PI + dx / 2.0 + (float)col * dx;
        float y = -PI + dx / 2.0 + (float)row * dx;
        float z = -PI + dx / 2.0 + (float)slice * dx;

        float f = TPMS_F[opts->TPMS_Type - 1](x,y,z);

        if(f < opts->isoValues && f > -opts->isoValues)
        {
            subDomain[i] = 0;
        } else if(f >= opts->isoValues)
        {
            if (mesh->sdInfo[subDomain[i] - 1].FC)
                subDomain[i] = 1;
            else
                subDomain[i] = 3;
        }else if(f <= -opts->isoValues)
        {
            if (mesh->sdInfo[subDomain[i] - 1].FC)
                subDomain[i] = 2;
            else
                subDomain[i] = 3;
        }
        else{
            printf("Error on subDomain correction!\n");
        }
    }

    return;
}

void DC_loc_debug(options *opts, meshInfo *mesh, float *DC)
{
    /*
     *  Function DC_loc_debug:
     *  Inputs:
     *      - options opts
     *      - mesh
     *      - pointer to DC
     *  Outputs:
     *      - none
     *
     *  Function will calculate the exepected isoValues
     *  of the DC[i] == 1. These values are printed
     *  to a file.
     *
     * */

    FILE *TEST = fopen("DC_loc_debug.csv", "w+");
    fprintf(TEST, "x,y,z,Iso\n");

    long int solidCt = 0;
    long int outOfBound_Ct = 0;

    int slice, row, col;

    for(long int i = 0; i < mesh->nElements; i++)
    {
        if(fabs(DC[i] - 1) > 0.5)
            continue;
       
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX * mesh->numCellsY) / mesh->numCellsX;
        col = (i - slice * mesh->numCellsY * mesh->numCellsX - row * mesh->numCellsX);

        // Now get x,y,z numbers
        float x = -PI + mesh->dx / 2.0 + (float)col * mesh->dx;
        float y = -PI + mesh->dx / 2.0 + (float)row * mesh->dx;
        float z = -PI + mesh->dx / 2.0 + (float)slice * mesh->dx;

        float f_val = TPMS_F[opts->TPMS_Type - 1](x,y,z);

        solidCt++;

        if(f_val > opts->isoValues || f_val < (-opts->isoValues))
        {
           outOfBound_Ct++;
           fprintf(TEST, "%1.3f,%1.3f,%1.3f,%1.3f\n", x, y, z, f_val); 
        } 
    }

    printf("Solid Total: %ld, OOB Total: %ld\n", solidCt, outOfBound_Ct);

    fclose(TEST);

    return;
}

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

void saveTemp_SF(float *DC, float *Temperature, meshInfo* mesh)
{
    /*
        Function saveTemp_SF:
        Inputs:
            - DC (array w/ diffusion coefficients)
            - Temperature (array with temperature distributions. I know, the name...)
            - meshInfo struct
        Outputs:
            - none

        This is a function to save temperatures in the shape factor simulation. Outputs
        are saved to sfTemp.csv
    */

    // save to see temperatures calculated

    FILE *TEST = fopen("sfTemp.csv", "w+");

    fprintf(TEST, "x,y,z,T\n");

    int row, col, slice;

    for(int i = 0; i < mesh->nElements; i++)
    {
        if(fabs(DC[i]-1) > 0.5)
            continue;
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX * mesh->numCellsY) / mesh->numCellsX;
        col = (i - slice * mesh->numCellsY * mesh->numCellsX - row * mesh->numCellsX);
        fprintf(TEST, "%d,%d,%d,%1.3f\n", col, row, slice, Temperature[i]);
    }

    fclose(TEST);

    return;
}

void debug_BoundariesIB(options *opts, meshInfo *mesh, IBM_Correct *IB, float *DC, float *Temperature)
{
    /*
        Function debug_BoundariesIB:
        Inputs:
            - pointer to opts
            - pointer to mesh
            - pointer to IBs
            - pointer to DC
            - pointer to Temperatures
        Outputs:
            - none.
        
        Function will map out boundaries and print the IB distances, and their respective heat flux
    */

    int slice, row, col;

    int nRows, nCols, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    bool isInt;

    float qw, qe, qn, qs, qb, qf;

    FILE *TEST = fopen("test_IB.csv", "w+");
    fprintf(TEST, "x,y,z,dw,de,ds,dn,db,df,qw,qe,qs,qn,qb,qf\n");

    // main loop:

    for (long int i = 0; i < mesh->nElements; i++)
    {
        // if not participating media, ignore

        if (fabs(DC[i]-1) > 0.5)
            continue;
        
        // reset variables
        isInt = false;

        qw = 0;
        qe = 0;
        qs = 0;
        qn = 0;
        qb = 0;
        qf = 0;

        // get index components
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX * mesh->numCellsY) / mesh->numCellsX;
        col = (i - slice * mesh->numCellsY * mesh->numCellsX - row * mesh->numCellsX);

        // check all faces for neighbors

        if (col == 0)
        {
            // Periodic West
            if (fabs(DC[i + nCols - 1] - 2) < 0.5)
            {
                qw = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
                isInt = true;
            }
            else if ( fabs(DC[i + nCols - 1] - 3) < 0.5)
            {
                qw = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
                isInt = true;
            }

            // East
            if (fabs(DC[i + 1] - 2) < 0.5)
            {
                qe = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
                isInt = true;
            }
            else if (fabs(DC[i + 1]-3) < 0.5)
            {
                qe = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
                isInt = true;
            }
        }
        else if (col == mesh->numCellsX - 1)
        {
            // West
            if (fabs(DC[i - 1] - 2) < 0.5)
            {
                qw = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
                isInt = true;
            }
            else if (fabs(DC[i - 1] - 3) < 3)
            {
                qw = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
                isInt = true;
            }

            // Periodic East
            if (fabs(DC[i - (mesh->numCellsX - 1)] - 2) < 0.5)
            {
                qe = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
                isInt = true;
            }
            else if (fabs(DC[i - (mesh->numCellsX - 1)] - 3) < 0.5)
            {
                qe = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
                isInt = true;
            }
        }
        else
        {
            // West
            if (fabs(DC[i - 1] - 2) < 0.5)
            {
                qw = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
                isInt = true;
            }
            else if (fabs(DC[i - 1] - 3) < 0.5)
            {
                qw = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
                isInt = true;
            }

            // East
            if (fabs(DC[i + 1] - 2) < 0.5)
            {
                isInt = true;
                qe = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
            else if (fabs(DC[i + 1] - 3) < 0.5)
            {
                isInt = true;
                qe = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
        }

        if (row == 0)
        {
            // Periodic North
            if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qn = (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }
            else if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qn = (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }

            // South
            if (fabs(DC[slice * nCols * nRows + (row + 1) * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qs = (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
            else if (fabs(DC[slice * nCols * nRows + (row + 1) * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qs = (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
        }
        else if (row == nRows - 1)
        {
            // North
            if (fabs(DC[slice * nCols * nRows + (row - 1) * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qn = (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }
            else if (fabs(DC[slice * nCols * nRows + (row - 1) * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qn = (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }

            // Periodic South
            if (fabs(DC[slice * nCols * nRows + col] - 2) < 0.5)
            {
                isInt = true;
                qs = (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
            else if (fabs(DC[slice * nCols * nRows + col] - 3) < 0.5)
            {
                isInt = true;
                qs = (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
        }
        else
        {
            // North
            if (fabs(DC[slice * nCols * nRows + (row - 1) * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qn = (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }
            else if (fabs(DC[slice * nCols * nRows + (row - 1) * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qn = (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }

            // South
            if (fabs(DC[slice * nCols * nRows + (row + 1) * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qs = (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
            else if (fabs(DC[slice * nCols * nRows + (row + 1) * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qs = (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
        }

        if (slice == 0)
        {
            // Periodic Front
            if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qf = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }
            else if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qf = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }

            // Back
            if (fabs(DC[(slice + 1) * nRows * nCols + row * nCols + col] - 2)  < 0.5)
            {
                isInt = true;
                qb = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
            else if (fabs(DC[(slice + 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qb = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
        }
        else if (slice == nSlices - 1)
        {
            // Front
            if (fabs(DC[(slice - 1) * nRows * nCols + row * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qf = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }
            else if (fabs(DC[(slice - 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qf = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }

            // Periodic Back
            if (fabs(DC[row * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qb = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
            else if (fabs(DC[row * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qb = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
        }
        else
        {
            // Front
            if (fabs(DC[(slice - 1) * nRows * nCols + row * nCols + col] - 2) < 0.5)
            {
                isInt = true;
                qf = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }
            else if (fabs(DC[(slice - 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qf = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }

            // Back
            if (fabs(DC[(slice + 1) * nRows * nCols + row * nCols + col] - 2) < 0.5 )
            {
                isInt = true;
                qb = (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
            else if (fabs(DC[(slice + 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
            {
                isInt = true;
                qb = (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
        }


        if(isInt == true)
        {
            // print to file
            fprintf(TEST, "%d,%d,%d,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e\n",
            col, row, slice, IB[i].Dist[0], IB[i].Dist[1], IB[i].Dist[2], IB[i].Dist[3], IB[i].Dist[4], IB[i].Dist[5], qw, qe, qs, qn, qb, qf);
        }
    } // end for

    fclose(TEST);

    return;
}

/*

    Initializing Simulation Temperature Distribution:

*/

void TemperatureInit(float *Temperature, float *DC, meshInfo *mesh, options *opts)
{
    /*
    
        Function TemperatureInit:

        Inputs:
            - Pointer to temperature distribution array
            - pointer to DC array (boundary flags)
            - pointer to meshInfo struct
            - pointer to opts struct
        
        Outputs:
            - none
        
        Function will initialize the temperature as a linear distribution
        between CLeft and CRight, depending on the element
        distance from CLeft or CRight.
    
        Strategy: Calculate EDT for POI relative to CLeft and relative to CRight
        separately. Then, linearly interpolate what the temperature should be
        based on the difference between the distances.

        For distances DL and DR, assume the total distance to boundaries is DL + DR.
        Then, the temperature will be:

        T = opts->CLeft + (DL) / (DL + DR) * opts->CRight;
    */

    // Open two boolean arrays:

    bool *Bool_R = (bool *)malloc(sizeof(bool) * mesh->nElements);
    bool *Bool_L = (bool *)malloc(sizeof(bool) * mesh->nElements);

    // Two arrays to store distances

    float *EDT_R = (float *)malloc(sizeof(float) * mesh->nElements);
    float *EDT_L = (float *)malloc(sizeof(float) * mesh->nElements);

    // parallel programming stuff

    omp_set_num_threads(opts->nThreads);

    // Initialize

    memset(Bool_L, 0, sizeof(bool) * mesh->nElements);
    memset(Bool_R, 0, sizeof(bool) * mesh->nElements);

    memset(EDT_L, 0, sizeof(float) * mesh->nElements);
    memset(EDT_R, 0, sizeof(float) * mesh->nElements);

    // populate bool arrays

    #pragma omp parallel for schedule(auto)
    for(long int i = 0; i < mesh->nElements; i++)
    {
        if(fabs(DC[i] - 3) < 0.5)
            Bool_L[i] = 1;
        else if(fabs(DC[i] - 2) < 0.5)
            Bool_R[i] = 1;
    }

    // Calculate EDT according to Meijster's algorithm

    pMeijster3D(Bool_L, EDT_L, mesh, 1);

    pMeijster3D(Bool_R, EDT_R, mesh, 1);

    // Initialize Temperature array

    #pragma omp parallel for schedule(auto)
    for(long int i = 0; i < mesh->nElements; i++)
    {
        if(fabs(DC[i] - 1) < 0.5)
            Temperature[i] = 0.0 + EDT_L[i]/(EDT_R[i] + EDT_L[i]) * 1.0;
    }

    // Memory management
    free(EDT_R);
    free(EDT_L);
    free(Bool_R);
    free(Bool_L);

    return;
}

float bisectionIBM(options *opts, meshInfo *mesh, int col_idx, int row_idx, int slice_idx, int nb_col, int nb_row, int nb_slice, int dir_search)
{
    /*
        Function bisectionIBM:
        Inputs:
            - pointer to opts
            - pointer to mesh
            - idx of current location
            - idx of neighbor
            - direction of search
        Outputs:
            - corrected delta

        Function appplies the bisection method to find the exact
        location of the interface in relation to the center of the 
        two adjacent cells.
     */

    // setup variables

    Point interceptPoint;

    Point a, b, c;

    Point p1, p2;

    p1.x = (float)col_idx*mesh->dx + mesh->dx/2.0 - PI;
    p1.y = (float)row_idx*mesh->dy + mesh->dy/2.0 - PI;
    p1.z = (float)slice_idx*mesh->dz + mesh->dz/2.0 - PI;

    p2.x = (float)nb_col*mesh->dx + mesh->dx/2.0 - PI;
    p2.y = (float)nb_row*mesh->dy + mesh->dy/2.0 - PI;
    p2.z = (float)nb_slice*mesh->dz + mesh->dz/2.0 - PI;

    float p1_value = zeroFunction(opts, p1);
    float p2_value = zeroFunction(opts, p2);

    float tol = 1e-6;
    int maxIter = 1e4;
    int iteration = 0;

    float eps, c1;

    float delta = mesh->dx/2;

    // check if either p1 or p2 are below tol

    if(fabs(p1_value) <= tol)
        return mesh->dx/2000;
    else if(fabs(p2_value) <= tol)
        return mesh->dx;

    // if p1.v or p2.v is within eps, return

    if(p1_value < 0 && p2_value > 0)
    {
        a = p1;
        b = p2;
    }
    else if(p1_value > 0 && p2_value < 0)
    {
        a = p2;
        b = p1;
    }
    else
    {
        printf("Error! Intercept is out of range!\n");
        printf("Values v1 = %1.3e, v2 = %1.3e\n", p1_value, p2_value);
        return delta;
    }

    // find c

    c.x = (p1.x + p2.x)/2.0;
    c.y = (p1.y + p2.y)/2.0;
    c.z = (p1.z + p2.z)/2.0;

    eps = fabs(zeroFunction(opts, c));

    if(eps <= tol)
    {
        if(dir_search == 0)
        {
            c1 = fabs(c.x - p1.x);
        }
        else if(dir_search == 1)
        {
            c1 = fabs(c.y - p1.y);
        }
        else if(dir_search == 2)
        {
            c1 = fabs(c.z - p1.z);
        }

        if(c1 < mesh->dx/2000)
        {
            c1 = mesh->dx/2000;
        }
        return c1;
    }

    float cVal;

    // otherwise, iterate
    while (eps > tol && iteration < maxIter)
    {
        cVal = zeroFunction(opts, c);
        // assign new bounds
        if(cVal > 0)
            b = c;
        else if(cVal < 0)
            a = c;
        
        // find new c

        c.x = (a.x + b.x)/2;
        c.y = (a.y + b.y)/2;
        c.z = (a.z + b.z)/2;

        eps = fabs(zeroFunction(opts, c));

        iteration++;
    }
    // assign intercept point
    interceptPoint = c;

    // use know direction to find the delta

    if(dir_search == 0)
    {
        delta = fabs(interceptPoint.x - p1.x);
    }
    else if(dir_search == 1)
    {
        delta = fabs(interceptPoint.y - p1.y);
    }
    else if(dir_search == 2)
    {
        delta = fabs(interceptPoint.z - p1.z);
    }

    if(delta < mesh->dx/2000)
    {
        delta = mesh->dx/2000;
    }

    return delta;
}

float temp_opt(options *opts, meshInfo *mesh, int col_idx, int row_idx, int slice_idx, int nb_col, int nb_row, int nb_slice, int dir_search)
{
    /*
    
        This is just a placeholder function to test the code out.

        if dir = 0, x
        if dir = 1, y
        if dir = 2, z
    
    */

    float eps = 1.000;
    float ib_delta = 0.0;
    
    int steps = 1000;

    float step_size = 0;

    int lb, ub;

    if(dir_search == 0)
    {
        step_size = fabs(col_idx - nb_col)/steps;

        if(nb_col < col_idx)
        {
            ub = col_idx;
            lb = nb_col;
        }
        else if(nb_col > col_idx)
        {
            ub = nb_col;
            lb = col_idx;
        }
        else
        {
            printf("Big mistake, deltaX not calculated properly\n");
            return 1;
        }

        for(float i = lb; i <= ub; i+=step_size)
        {
            float x_coord = i*mesh->dx + mesh->dx/2 - PI;
            float y_coord = row_idx * mesh->dy + mesh->dy/2 - PI;
            float z_coord = slice_idx * mesh->dz + mesh->dz/2 - PI;

            float new_eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](x_coord, y_coord, z_coord));
            
            if(new_eps < eps)
            {
                eps = new_eps;
                ib_delta = i/steps * mesh->dx;
            }
        }
    }
    else if(dir_search == 1)
    {
        step_size = fabs(row_idx - nb_row)/steps;

        if(nb_row < row_idx)
        {
            ub = row_idx;
            lb = nb_row;
        }
        else if(nb_row > row_idx)
        {
            ub = nb_row;
            lb = row_idx;
        }
        else
        {
            printf("Big Mistake, deltaY not calculated properly\n");
            return 1;
        }

        for(float i = lb; i <= ub; i+=step_size)
        {
            float x_coord = col_idx * mesh->dx + mesh->dx/2 - PI;
            float y_coord = i*mesh->dy + mesh->dy/2 - PI;
            float z_coord = slice_idx * mesh->dz + mesh->dz/2 - PI;

            float new_eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](x_coord, y_coord, z_coord));

            if(new_eps < eps)
            {
                eps = new_eps;
                ib_delta = i/steps * mesh->dy;
            }
        }
    }
    else if(dir_search == 2)
    {
        step_size = fabs(slice_idx - nb_slice)/steps;

        if(nb_slice < slice_idx)
        {
            ub = slice_idx;
            lb = nb_slice;
        }
        else if(nb_slice > slice_idx)
        {
            ub = nb_slice;
            lb = slice_idx;
        }
        else
        {
            printf("Big Mistake, deltaY not calculated properly\n");
            return 1;
        }

        for(float i = lb; i <= ub; i+=step_size)
        {
            float x_coord = col_idx * mesh->dx + mesh->dx/2 - PI;
            float y_coord = row_idx*mesh->dy + mesh->dy/2 - PI;
            float z_coord = i * mesh->dz + mesh->dz/2 - PI;
            float new_eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](x_coord, y_coord, z_coord));

            if(new_eps < eps)
            {
                eps = new_eps;
                ib_delta = i/steps * mesh->dz;
            }
        }
    }

    // correct to avoid NaN's

    if(ib_delta < mesh->dx/2000)
    {
        ib_delta += mesh->dx/2000;
    }
    
    return ib_delta;
}

void correctIBs(float *DC, options *opts, meshInfo *mesh, IBM_Correct *IB)
{
    /*
        correctIBs Function:
        Inputs:
            - pointer to DC
            - pointer to meshInfo
            - pointer to IB table
        Outputs:
            - none
        Function will correct the distances to IBs where needed.
        If not needed, just populate IB with corresponding dx,
        dy, or dz.
    */

    int nCols, nRows, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;
    long int count = 0;
    long int bd_ct = 0;
    for(int i = 0; i < mesh->nElements; i++)
    {
        // initialize IB-Dists
        IB[i].Dist[0] = mesh->dx;
        IB[i].Dist[1] = mesh->dx;
        IB[i].Dist[2] = mesh->dy;
        IB[i].Dist[3] = mesh->dy;
        IB[i].Dist[4] = mesh->dz;
        IB[i].Dist[5] = mesh->dz;
        // check if fluid, no need to do anything
        if (fabs(DC[i]-1.0) > 0.5)
            continue;
        
        /*
            Solid voxel:
            - if boundary, correct dx, dy, dz
            - if not boundary, dx, dy, dz remain
                unchanged
        */

        int row, col, slice;
        slice = i / (nRows * nCols);
        row = (i - slice * nRows * nCols) / nCols;
        col = (i - slice * nRows * nCols - row * nCols);

        // West
        long int nbIdx = 0;
        if(col == 0)
        {
            // periodic west
            nbIdx = slice * nRows * nCols + row * nCols + (nCols - 1); 
            if ( fabs(DC[nbIdx] - DC[i]) > 0.5 )
            {
                bd_ct++;
                // Update dx_w accordingly
                float dx_w = bisectionIBM(opts, mesh, col, row, slice, -1, row, slice, 0); // special periodic dist -1
                IB[i].Dist[0] = dx_w;
                if(dx_w == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[0] = mesh->dx;
            }
        } else
        {
            // not periodic West
            nbIdx = i - 1;
            if (fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dx_w = bisectionIBM(opts, mesh, col, row, slice, (col-1), row, slice, 0);
                IB[i].Dist[0] = dx_w;
                if(dx_w == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[0] = mesh->dx;
            }
        }

        // East
        if(col == nCols - 1)
        {
            // periodic East
            nbIdx = slice * nRows * nCols + row * nCols + 0;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dx_e = bisectionIBM(opts, mesh, col, row, slice, nCols, row, slice, 0); // special periodic dist nCols
                IB[i].Dist[1] = dx_e;
                if(dx_e == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[1] = mesh->dx;
            }
        } else
        {
            // not periodic East
            nbIdx = i+1;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dx_e = bisectionIBM(opts, mesh, col, row, slice, (col+1), row, slice, 0);
                IB[i].Dist[1] = dx_e;
                if(dx_e == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[1] = mesh->dx;
            }
        }

        // South

        if(row == nRows - 1)
        {
            // Periodic South
            nbIdx = slice * nRows * nCols + (0) * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dy_s = bisectionIBM(opts, mesh, col, row, slice, col, nRows, slice, 1); // special periodic dist nRows
                IB[i].Dist[2] = dy_s;
                if(dy_s == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[2] = mesh->dy;
            }
        }
        else
        {
            // Not Periodic South
            nbIdx = slice * nRows * nCols + (row+1) * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dy_s = bisectionIBM(opts, mesh, col, row, slice, col, row + 1, slice, 1);
                IB[i].Dist[2] = dy_s;
                if(dy_s == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[2] = mesh->dy;
            }
        }

        // North

        if(row == 0)
        {
            // Periodic North
            nbIdx = slice * nRows * nCols + (nRows - 1) * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dy_n = bisectionIBM(opts, mesh, col, row, slice, col, -1, slice, 1); // special periodic dist -1
                IB[i].Dist[3] = dy_n;
                if(dy_n == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[3] = mesh->dy;
            }
        }
        else
        {
            // Not periodic North
            nbIdx = slice * nRows * nCols + (row - 1) * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dy_n = bisectionIBM(opts, mesh, col, row, slice, col, row - 1, slice, 1);
                IB[i].Dist[3] = dy_n;
                if(dy_n == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[3] = mesh->dy;
            }
        }

        // Back

        if(slice == nSlices - 1)
        {
            // Periodic Back
            nbIdx = (0) * nRows * nCols + row * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dz_b = bisectionIBM(opts, mesh, col, row, slice, col, row, nSlices, 2); // special periodic dist nSlices
                IB[i].Dist[4] = dz_b;
                if(dz_b == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[4] = mesh->dz;
            }
        }
        else
        {
            // Not Periodic Back
            nbIdx = (slice + 1) * nRows * nCols + row * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dz_b = bisectionIBM(opts, mesh, col, row, slice, col, row, slice + 1, 2);
                IB[i].Dist[4] = dz_b;
                if(dz_b == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[4] = mesh->dz;
            }
        }
    
        // Front

        if(slice == 0)
        {
            // Periodic Front
            nbIdx = (nSlices - 1) * nRows * nCols + row * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dz_f = bisectionIBM(opts, mesh, col, row, slice, col, row, -1, 2); // special periodic dist -1
                IB[i].Dist[5] = dz_f;
                if(dz_f == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[5] = mesh->dz;
            }
        }
        else
        {
            // Not Periodic Front
            nbIdx = (slice - 1) * nRows * nCols + row * nCols + col;
            if(fabs(DC[nbIdx] - DC[i]) > 0.5)
            {
                bd_ct++;
                float dz_f = bisectionIBM(opts, mesh, col, row, slice, col, row, slice - 1, 2);
                IB[i].Dist[5] = dz_f;
                if(dz_f == mesh->dx/2)
                    count++;
            }
            else
            {
                IB[i].Dist[5] = mesh->dz;
            }
        }
    
    } // end of for

    printf("Count Errors: %ld, Count total = %ld \n", count, bd_ct);
    return;
}

/*

    Discretization Functions

*/

void SetDC_SF(float *DC, char *subDomain, meshInfo *mesh)
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
        else if(subDomain[i] == 3)
        {
            DC[i] = 4;
        }
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

        if (DC[i] == 1 || DC[i] == 4)
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
            if (fabs(DC[i - 1 + nCols] - 1) < 0.5)
            {
                // Left boundary, participating media
                dw = DC[i];
                CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
            }
            else if (fabs(DC[i - 1 + nCols] - 2) < 0.5)
            {
                // Left boundary, first channel
                dw = DC[i];
                RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
            else if (fabs(DC[i - 1 + nCols] - 3) < 0.5)
            {
                // Left boundary, second channel
                dw = DC[i];
                RHS[i] -= opts->CRight * dw * (dy * dz) / (dx / 2);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
        }
        else
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i - 1] - 1) < 0.5)
            {
                // West is participating media
                dw = DC[i];
                CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
            }
            else if (fabs(DC[i - 1] - 2) < 0.5)
            {
                // West is first channel
                dw = DC[i];
                RHS[i] -= opts->CLeft * dw * (dy * dz) / (dx / 2);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (dx / 2);
            }
            else if (fabs(DC[i - 1] - 3) < 0.5)
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
            if (fabs(DC[i + 1 - nCols] - 1) < 0.5)
            {
                // Right boundary, participating media
                de = DC[i];
                CoeffMatrix[i * 7 + 2] = de * (dy * dz) / (dx);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx);
            }
            else if (fabs(DC[i + 1 - nCols] - 2) < 0.5)
            {
                // Right boundary, first channel
                de = DC[i];
                RHS[i] -= opts->CLeft * de * (dy * dz) / (dx / 2);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
            if (fabs(DC[i + 1 - nCols] - 3) < 0.5)
            {
                // Right boundary, second channel
                de = DC[i];
                RHS[i] -= opts->CRight * de * (dy * dz) / (dx / 2);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
        }
        else
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i + 1] - 1) < 0.5)
            {
                // East, participating media
                de = DC[i];
                CoeffMatrix[i * 7 + 2] = de * (dy * dz) / (dx);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx);
            }
            else if (fabs(DC[i + 1] - 2) < 0.5)
            {
                // East, first channel
                de = DC[i];
                RHS[i] -= opts->CLeft * de * (dy * dz) / (dx / 2);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx / 2);
            }
            if (fabs(DC[i + 1] - 3) < 0.5)
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
            if (fabs(DC[i + nCols] - 1) < 0.5)
            {
                // South, participating
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
            else if (fabs(DC[i + nCols] - 2) < 0.5)
            {
                // South, first channel
                ds = DC[i];
                RHS[i] -= opts->CLeft * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / (dy / 2);
            }
            else if (fabs(DC[i + nCols] - 3) < 0.5)
            {
                // South, second channel
                ds = DC[i];
                RHS[i] -= opts->CRight * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / (dy / 2);
            }
        } else
        {
            // Periodic Boundary
            if (fabs(DC[slice * nCols * nRows + col] - 1) < 0.5) //Periodic Boundary
            {
                // South, participating
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / dy;
            }
            else if (fabs(DC[slice * nCols * nRows + col] - 2) < 0.5)
            {
                // South, first channel
                ds = DC[i];
                RHS[i] -= opts->CLeft * ds * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / (dy / 2);
            }
            else if (fabs(DC[slice * nCols * nRows + col] - 3) < 0.5)
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
            if (fabs(DC[i - nCols] - 1) < 0.5)
            {
                // North, participating media
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
            else if (fabs(DC[i - nCols] - 2) < 0.5)
            {
                // North, first channel
                dn = DC[i];
                RHS[i] -= opts->CLeft * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            }
            else if (fabs(DC[i - nCols] - 3) < 0.5)
            {
                // North, second channel
                dn = DC[i];
                RHS[i] -= opts->CRight * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            }
        } else
        {
            // Periodic Boundary
            if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 1) < 0.5) 
            {
                //North, participating media
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            } else if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 2) < 0.5)
            {
                // North, first channel
                dn = DC[i];
                RHS[i] -= opts->CLeft * dn * (dx * dz) / (dy / 2);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (dy / 2);
            } else if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 3) < 0.5)
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
            if (fabs(DC[i + nCols * nRows] - 1) < 0.5)
            {
                // Back, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (fabs(DC[i + nCols * nRows] - 2) < 0.5)
            {
                // Back, first channel
                db = DC[i];
                RHS[i] -= opts->CLeft * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
            else if (fabs(DC[i + nCols * nRows] - 3) < 0.5)
            {
                // Back, second channel
                db = DC[i];
                RHS[i] -= opts->CRight * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
        } else
        {
            //Periodic Boundary
            if (fabs(DC[row * nCols + col] - 1) < 0.5 )
            {
                // Back, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (fabs(DC[row * nCols + col] - 2) < 0.5)
            {
                // Periodic Back
                db = DC[i];
                RHS[i] -= opts->CLeft * db * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (dz / 2);
            }
            else if (fabs(DC[row * nCols + col] - 3) < 0.5)
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
            if (fabs(DC[i - nCols * nRows] - 1) < 0.5)
            {
                // Front, participating media
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
            else if (fabs(DC[i - nCols * nRows] - 2) < 0.5)
            {
                // Front, first channel
                df = DC[i];
                RHS[i] -= opts->CLeft * df * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (dz / 2);
            }
            else if (fabs(DC[i - nCols * nRows] - 3) < 0.5)
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
            if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 1) < 0.5)
            {
                // Front, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 6] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 2) < 0.5)
            {
                // Front, first channel
                df = DC[i];
                RHS[i] -= opts->CLeft * df * (dx * dy) / (dz / 2);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (dz / 2);
            }
            else if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
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

int Disc3D_IB_SF(options *opts,
                meshInfo *mesh,
                float *DC,
                float *CoeffMatrix,
                float *RHS,
                IBM_Correct *IB)
{
    /*
        Function Disc3D_SF_PB:
        Inputs:
            - pointer to options data structure
            - pointer to mesh data structure
            - pointer to float array DC holding diffusion coefficients
            - pointer to float array CoeffMatrix Coefficient Matrix
            - pointer to float array RHS holding right-hand side of discretized system.
            - pointer to immersed boundary corrected distances
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

        if (DC[i] == 1 || DC[i] == 4 )
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
            if (fabs(DC[i - 1 + nCols] - 1) < 0.5)
            {
                // Left boundary, participating media
                dw = DC[i];
                CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
            }
            else if (fabs(DC[i - 1 + nCols] - 2) < 0.5 )
            {
                // Left boundary, first channel
                dw = DC[i];
                RHS[i] -= opts->CLeft * dw * (dy * dz) / (IB[i].Dist[0]);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (IB[i].Dist[0]);
            }
            else if (fabs(DC[i - 1 + nCols] - 3) < 0.5 )
            {
                // Left boundary, second channel
                dw = DC[i];
                RHS[i] -= opts->CRight * dw * (dy * dz) / (IB[i].Dist[0]);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (IB[i].Dist[0]);
            }
        }
        else
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i - 1] - 1) < 0.5 )
            {
                // West is participating media
                dw = DC[i];
                CoeffMatrix[i * 7 + 1] = dw * (dy * dz) / dx;
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / dx;
            }
            else if (fabs(DC[i - 1] - 2) < 0.5 )
            {
                // West is first channel
                dw = DC[i];
                RHS[i] -= opts->CLeft * dw * (dy * dz) / (IB[i].Dist[0]);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (IB[i].Dist[0]);
            }
            else if (fabs(DC[i - 1] - 3) < 0.5 )
            {
                // West is second channel
                dw = DC[i];
                RHS[i] -= opts->CRight * dw * (dy * dz) / (IB[i].Dist[0]);
                CoeffMatrix[i * 7 + 0] -= dw * (dy * dz) / (IB[i].Dist[0]);
            }
        }

        // East

        if (col == mesh->numCellsX - 1)
        {
            // Periodic Boundary
            if (fabs(DC[i + 1 - nCols] - 1) < 0.5 )
            {
                // Right boundary, participating media
                de = DC[i];
                CoeffMatrix[i * 7 + 2] = de * (dy * dz) / (dx);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx);
            }
            else if (fabs(DC[i + 1 - nCols] - 2) < 0.5 )
            {
                // Right boundary, first channel
                de = DC[i];
                RHS[i] -= opts->CLeft * de * (dy * dz) / (IB[i].Dist[1]);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (IB[i].Dist[1]);
            }
            if (fabs(DC[i + 1 - nCols] - 3) < 0.5 )
            {
                // Right boundary, second channel
                de = DC[i];
                RHS[i] -= opts->CRight * de * (dy * dz) / (IB[i].Dist[1]);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (IB[i].Dist[1]);
            }
        }
        else
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i + 1] - 1) < 0.5 )
            {
                // East, participating media
                de = DC[i];
                CoeffMatrix[i * 7 + 2] = de * (dy * dz) / (dx);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (dx);
            }
            else if (fabs(DC[i + 1] - 2) < 0.5)
            {
                // East, first channel
                de = DC[i];
                RHS[i] -= opts->CLeft * de * (dy * dz) / (IB[i].Dist[1]);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (IB[i].Dist[1]);
            }
            if (fabs(DC[i + 1] - 3) < 0.5 )
            {
                // East, second channel
                de = DC[i];
                RHS[i] -= opts->CRight * de * (dy * dz) / (IB[i].Dist[1]);
                CoeffMatrix[i * 7 + 0] -= de * (dy * dz) / (IB[i].Dist[1]);
            }
        }

        // South

        if (row != mesh->numCellsY - 1)
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i + nCols] - 1) < 0.5 )
            {
                // South, participating
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / dy;
            }
            else if (fabs(DC[i + nCols] - 2) < 0.5 )
            {
                // South, first channel
                ds = DC[i];
                RHS[i] -= opts->CLeft * ds * (dx * dz) / (IB[i].Dist[2]);
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / (IB[i].Dist[2]);
            }
            else if (fabs(DC[i + nCols] - 3) < 0.5)
            {
                // South, second channel
                ds = DC[i];
                RHS[i] -= opts->CRight * ds * (dx * dz) / (IB[i].Dist[2]);
                CoeffMatrix[i * 7 + 0] -= ds * (dx * dz) / (IB[i].Dist[2]);
            }
        } else
        {
            // Periodic Boundary
            if (fabs(DC[slice * nCols * nRows + col] - 1) < 0.5) //Periodic Boundary
            {
                // South, participating
                ds = DC[i];
                CoeffMatrix[i * 7 + 3] = ds * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / dy;
            }
            else if (fabs(DC[slice * nCols * nRows + col] - 2) < 0.5 )
            {
                // South, first channel
                ds = DC[i];
                RHS[i] -= opts->CLeft * ds * (dx * dz) / (IB[i].Dist[2]);
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / (IB[i].Dist[2]);
            }
            else if (fabs(DC[slice * nCols * nRows + col] - 3) < 0.5 )
            {
                // South, second channel
                ds = DC[i];
                RHS[i] -= opts->CRight * ds * (dx * dz) / (IB[i].Dist[2]);
                CoeffMatrix[i * 7 + 0] -= ds * (dx *dz) / (IB[i].Dist[2]);
            }
        }

        // North

        if (row != 0)
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i - nCols] - 1) < 0.5 )
            {
                // North, participating media
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            }
            else if (fabs(DC[i - nCols] - 2) < 0.5)
            {
                // North, first channel
                dn = DC[i];
                RHS[i] -= opts->CLeft * dn * (dx * dz) / (IB[i].Dist[3]);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (IB[i].Dist[3]);
            }
            else if (fabs(DC[i - nCols] - 3) < 0.5)
            {
                // North, second channel
                dn = DC[i];
                RHS[i] -= opts->CRight * dn * (dx * dz) / (IB[i].Dist[3]);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (IB[i].Dist[3]);
            }
        } else
        {
            // Periodic Boundary
            if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 1) < 0.5 ) 
            {
                //North, participating media
                dn = DC[i];
                CoeffMatrix[i * 7 + 4] = dn * (dx * dz) / dy;
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / dy;
            } else if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 2) < 0.5 )
            {
                // North, first channel
                dn = DC[i];
                RHS[i] -= opts->CLeft * dn * (dx * dz) / (IB[i].Dist[3]);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (IB[i].Dist[3]);
            } else if (fabs(DC[slice * nCols * nRows + (nRows - 1) * nCols + col] - 3) < 0.5 )
            {
                // North, second channel
                dn = DC[i];
                RHS[i] -= opts->CRight * dn * (dx * dz) / (IB[i].Dist[3]);
                CoeffMatrix[i * 7 + 0] -= dn * (dx * dz) / (IB[i].Dist[3]);
            }
        }

        // Back

        if (slice != mesh->numCellsZ - 1)
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i + nCols * nRows] - 1) < 0.5 )
            {
                // Back, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (fabs(DC[i + nCols * nRows] - 2) < 0.5 )
            {
                // Back, first channel
                db = DC[i];
                RHS[i] -= opts->CLeft * db * (dx * dy) / (IB[i].Dist[4]);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (IB[i].Dist[4]);
            }
            else if (fabs(DC[i + nCols * nRows] - 3) < 0.5 )
            {
                // Back, second channel
                db = DC[i];
                RHS[i] -= opts->CRight * db * (dx * dy) / (IB[i].Dist[4]);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (IB[i].Dist[4]);
            }
        } else
        {
            //Periodic Boundary
            if (fabs(DC[row * nCols + col] - 1) < 0.5 )
            {
                // Back, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 5] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (fabs(DC[row * nCols + col] - 2) < 0.5)
            {
                // Periodic Back
                db = DC[i];
                RHS[i] -= opts->CLeft * db * (dx * dy) / (IB[i].Dist[4]);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (IB[i].Dist[4]);
            }
            else if (fabs(DC[row * nCols + col] - 3) < 0.5)
            {
                // Periodic Back
                db = DC[i];
                RHS[i] -= opts->CRight * db * (dx * dy) / (IB[i].Dist[4]);
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / (IB[i].Dist[4]);
            }
        }

        // Front

        if (slice != 0)
        {
            // Non-Boundary Neighbor
            if (fabs(DC[i - nCols * nRows] - 1) < 0.5 )
            {
                // Front, participating media
                df = DC[i];
                CoeffMatrix[i * 7 + 6] = df * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / dz;
            }
            else if (fabs(DC[i - nCols * nRows] - 2) < 0.5 )
            {
                // Front, first channel
                df = DC[i];
                RHS[i] -= opts->CLeft * df * (dx * dy) / (IB[i].Dist[5]);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (IB[i].Dist[5]);
            }
            else if (fabs(DC[i - nCols * nRows] - 3) < 0.5)
            {
                // Front, second channel
                df = DC[i];
                RHS[i] -= opts->CRight * df * (dx * dy) / (IB[i].Dist[5]);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (IB[i].Dist[5]);
            }
        }
        else
        {
            // Periodic Boundary
            if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 1) < 0.5 )
            {
                // Front, participating media
                db = DC[i];
                CoeffMatrix[i * 7 + 6] = db * (dx * dy) / dz;
                CoeffMatrix[i * 7 + 0] -= db * (dx * dy) / dz;
            }
            else if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 2) < 0.5 )
            {
                // Front, first channel
                df = DC[i];
                RHS[i] -= opts->CLeft * df * (dx * dy) / (IB[i].Dist[5]);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (IB[i].Dist[5]);
            }
            else if (fabs(DC[(nSlices - 1) * nRows * nCols + row * nCols + col] - 3) < 0.5)
            {
                // Front, second channel
                df = DC[i];
                RHS[i] -= opts->CRight * df * (dx * dy) / (IB[i].Dist[5]);
                CoeffMatrix[i * 7 + 0] -= df * (dx * dy) / (IB[i].Dist[5]);
            }
        }
    }

    return 0;
}

/*

    Calculate Shape Factors:

*/

float Calc_Q_SF(options *opts, meshInfo *mesh, float *Temperature, float *DC)
{
    /*
        Function Calc_Q_SF:
        Inputs:
            - opts struct
            - mesh struct
            - Conc array, holding SF sim solutions
            - DC array, holding discretized sub-domain information
        Outputs:
            - (float) Q_avg: average heat flux in the domain
    */

    // Initialize variables
    float Q_Avg, Q_21, Q_13;

    Q_21 = 0;
    Q_13 = 0;

    int slice, row, col;

    int nRows, nCols, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    // main loop:

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
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i + nCols - 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }

            // East
            if (DC[i + 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i + 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
        }
        else if (col == mesh->numCellsX - 1)
        {
            // West
            if (DC[i - 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i - 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }

            // Periodic East
            if (DC[i - (mesh->numCellsX - 1)] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i - (mesh->numCellsX - 1)] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
        }
        else
        {
            // West
            if (DC[i - 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i - 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }

            // East
            if (DC[i + 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
            else if (DC[i + 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (mesh->dx / 2);
            }
        }

        if (row == 0)
        {
            // Periodic North
            if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }

            // South
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
        }
        else if (row == nRows - 1)
        {
            // North
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }

            // Periodic South
            if (DC[slice * nCols * nRows + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
        }
        else
        {
            // North
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }

            // South
            if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
            else if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (mesh->dy / 2);
            }
        }

        if (slice == 0)
        {
            // Periodic Front
            if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }

            // Back
            if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
        }
        else if (slice == nSlices - 1)
        {
            // Front
            if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }

            // Periodic Back
            if (DC[row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
        }
        else
        {
            // Front
            if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }

            // Back
            if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
            else if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (mesh->dz / 2);
            }
        }
    } // end of the for loop

    Q_Avg = (Q_13 + Q_21)/2;

    if (opts->verbose)
    {
        printf("Shape Factor Simulation Results:\n");
        printf("S = %1.3e, Q_21 = %1.3e, Q_13 = %1.3e, Pct. Error = %3.2f %%\n", Q_Avg/(opts->CLeft - opts->CRight), Q_21, Q_13, fabs((Q_21 - Q_13)/Q_21)*100);
    }

    return Q_Avg;
}


float Calc_Q_SF_IBs(options *opts, meshInfo *mesh, float *Temperature, float *DC, IBM_Correct *IB)
{
    /*
        Function Calc_Q_SF_IBs:
        Inputs:
            - opts struct
            - mesh struct
            - Conc array, holding SF sim solutions
            - DC array, holding discretized sub-domain information
            - Immersed boundary corrected distances
        Outputs:
            - (float) Q_avg: average heat flux in the domain
    */

    // Initialize variables
    float Q_Avg, Q_21, Q_13;

    Q_21 = 0;
    Q_13 = 0;

    int slice, row, col;

    int nRows, nCols, nSlices;
    nCols = mesh->numCellsX;
    nRows = mesh->numCellsY;
    nSlices = mesh->numCellsZ;

    // main loop:

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
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
            }
            else if (DC[i + nCols - 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
            }

            // East
            if (DC[i + 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
            else if (DC[i + 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
        }
        else if (col == mesh->numCellsX - 1)
        {
            // West
            if (DC[i - 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
            }
            else if (DC[i - 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
            }

            // Periodic East
            if (DC[i - (mesh->numCellsX - 1)] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
            else if (DC[i - (mesh->numCellsX - 1)] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
        }
        else
        {
            // West
            if (DC[i - 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
            }
            else if (DC[i - 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[0]);
            }

            // East
            if (DC[i + 1] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
            else if (DC[i + 1] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dz) / (IB[i].Dist[1]);
            }
        }

        if (row == 0)
        {
            // Periodic North
            if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }
            else if (DC[slice * nCols * nRows + (nRows - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }

            // South
            if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
            else if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
        }
        else if (row == nRows - 1)
        {
            // North
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }

            // Periodic South
            if (DC[slice * nCols * nRows + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
            else if (DC[slice * nCols * nRows + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
        }
        else
        {
            // North
            if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }
            else if (DC[slice * nCols * nRows + (row - 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[3]);
            }

            // South
            if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
            else if (DC[slice * nCols * nRows + (row + 1) * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dz * mesh->dx) / (IB[i].Dist[2]);
            }
        }

        if (slice == 0)
        {
            // Periodic Front
            if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }
            else if (DC[(nSlices - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }

            // Back
            if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
            else if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
        }
        else if (slice == nSlices - 1)
        {
            // Front
            if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }
            else if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }

            // Periodic Back
            if (DC[row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
            else if (DC[row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
        }
        else
        {
            // Front
            if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }
            else if (DC[(slice - 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[5]);
            }

            // Back
            if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 2)
            {
                Q_21 += (opts->CLeft - Temperature[i]) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
            else if (DC[(slice + 1) * nRows * nCols + row * nCols + col] == 3)
            {
                Q_13 += (Temperature[i] - opts->CRight) * (mesh->dy * mesh->dx) / (IB[i].Dist[4]);
            }
        }
    } // end of the for loop

    Q_Avg = (Q_13 + Q_21)/2;

    if (opts->verbose)
    {
        printf("Shape Factor Simulation Results:\n");
        printf("S = %1.3e, Q_21 = %1.3e, Q_13 = %1.3e, Pct. Error = %3.2f %%\n", Q_Avg/(opts->CLeft - opts->CRight), Q_21, Q_13, fabs((Q_21 - Q_13)/Q_21)*100);
    }

    return Q_Avg;
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

    // declare and define DC in the main flow channel

    float *DC = (float *)malloc(sizeof(float) * mesh->nElements);

    memset(DC, 0 , mesh->nElements * sizeof(float));

    // Sub-Domain Info is absolutely necessary for this simulation

    if(mesh->nFC < 2)
    {
        printf("Error Detected: nChannels = %d, nFC = %d\n", mesh->nChannels, mesh->nFC);
        printf("Returning.....");
        return 1;
    }

    fixSD_info(opts, mesh, subDomain);

    // Populate the array based on the structure

    SetDC_SF(DC, subDomain, mesh);

    //DC_loc_debug(opts, mesh, DC);

    // allocate the arrays for simulation

    float *CoeffMatrix = (float *)malloc(mesh->nElements * 7 * sizeof(float));
    float *RHS = (float *)malloc(mesh->nElements * sizeof(float));
    float *Temperature = (float *)malloc(mesh->nElements * sizeof(float));

    // initialize memory

    memset(CoeffMatrix, 0, sizeof(float) * 7 * mesh->nElements);
    memset(RHS, 0, mesh->nElements * sizeof(float));

    // Initializing temperature distribution
    memset(Temperature, 0, mesh->nElements * sizeof(float));

    // allocate the struct

    IBM_Correct *IBs = (IBM_Correct *)malloc(sizeof(IBM_Correct) * mesh->nElements);

    memset(IBs, 0, sizeof(float)*mesh->nElements*6);

    correctIBs(DC, opts, mesh, IBs);

    // debug

    // debug_BoundariesIB(opts, mesh, IBs, DC, Temperature);

    // Discretize

    Disc3D_IB_SF(opts, mesh, DC, CoeffMatrix, RHS, IBs);

    // Initialize temperature again

    if(opts->verbose)
        printf("Initializing Temperature\n");

    TemperatureInit(Temperature, DC, mesh, opts);

    // solve

    if(opts->useGPU)
    {
        printf("GPU not currently Supported. Using CPU instead.\n");
    }

    // CPU solve
    pGS3D_SF_handle(opts, mesh, save, CoeffMatrix, RHS, Temperature);

    // Calculate Shape Factor

    float Q_avg_IBs = Calc_Q_SF_IBs(opts, mesh, Temperature, DC, IBs);

    // Add debug function to print IBs and Q's for all boundaries

    // debug_BoundariesIB(opts, mesh, IBs, DC, Temperature);

    save->SF = Q_avg_IBs/(opts->CLeft - opts->CRight);
    
    // saveTemp

    if(opts->sfTMAP)
        saveTemp_SF(DC, Temperature, mesh);

    // memory management
    free(Temperature);
    free(CoeffMatrix);
    free(RHS);

    return 0;
}

#endif
