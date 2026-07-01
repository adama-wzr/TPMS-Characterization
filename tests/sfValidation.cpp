/*
 *
 * This is a test function to validate the shape factor code 
 * against some results that have analytic or theoretical solutions.
 *
 * I will put all the stuff that I need in this file, no .hpp
 *
 * However, as soon as I started this validation step, I realised
 * most of the shape-factor code needs to be slighly re-written, 
 * mainly due to the bounds of the problem not being +-C, like they 
 * are for the TPMS.
 *
 * On the other hand, we are currently constrained to +-C for the 
 * TPMS, so it might be worth it to convert opts->isoValues to a
 * tuple instead of a single float.
 *
 * Last modified: 06/24/2026
 * Andre Adam.
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <lib/data_structures.hpp> 
#include <usrInput.hpp>
#include <sizeDistributions.hpp>
#include <surfaceArea.hpp>
#include <marchingCubes.hpp>
#include <cpu_solvers/cpuSolvers.hpp>

/*
 *
 * Special Data Structures
 *
 * */

typedef struct{
    float SF_TH;
    float SF_sim;
    float area_MC;
    float area_Voxel;
} solutions_sf;

typedef struct{
    float L;
    float D1;       // D1 is "a" for square flow passage
    float D2;       // D2 is "b" for square flow passage
    float center;   // tipically the 4th input to th_case_f
} geometry;

typedef float (*th_case_f_ptr)(float, float, float, float);

/*
 *
 *      Implicit Functions:
 *
 * */

float sqPassage_F(float x, float y, float z, float w)
{
    float f = std::max(fabs(x),fabs(y)) - w/2.0;
    return f;
}

float cylTube_F(float x, float y, float z, float r)
{
    float f = x + y + z + r;
    return f;
}

/*
 *
 *      Lookup Tables
 *
 * */

static const char *TH_Cases[] =
{
    "Square Flow Passage", "Long Cylindrical Layer"
};

th_case_f_ptr TH_CaseF[] = {
    sqPassage_F,
    cylTube_F
};

/*
 *
 * Helper Functions (save/opts/print/generate)
 *
 * */

void printIndex(char max)
{
    printf("----------------------------------\n\n");
    printf("    Theoreical Case Index:\n\n");
    printf("------------------------------------\n\n");
    printf("Index      Name\n");
    printf("_____      ____\n");

    for(int i = 0; i < max; i++)
    {
        printf("%02d         %s\n", i+1, TH_Cases[i]);
    }

    printf("------------------------------------\n\n");

    return;
}

void generateFlowPassage(options *opts, meshInfo *mesh, geometry *geo, char *P)
{
    /*
     *  Function generateFlowPassage:
     *  
     *  Inputs:
     *      - pointer to options opts;
     *      - pointer to meshInfo mesh;
     *      - pointer to geometry geo;
     *      - pointer to char array holding structure;
     *  Outputs;
     *      - none
     *
     *  Function will create the flow passage according to the implicit equation
     *  and the user-entered information.
     *
     * */

    // define some stuff to make life easier
    float dx = mesh->dx;
    float dy = mesh->dy;
    float dz = mesh->dz;
    
    int nCols = mesh->numCellsX;
    int nRows = mesh->numCellsY;
    int nSlices = mesh->numCellsZ;

    float x, y, z;
    int col, row, slice;

    for(int i = 0; i < mesh->nElements; i++)
    {
        // get col, row, slice
        slice = i / (nCols * nRows);
        row = (i - slice * nCols * nRows) / nCols;
        col = i - slice * nRows * nCols - row * nCols;

        x = -PI + (float) col * dx + dx / 2.0;
        y = -PI + (float) row * dy + dy / 2.0;
        z = -PI + (float) slice * dz + dz / 2.0;

        float f = TH_CaseF[opts->TPMS_Type - 1](x,y,z, geo->center);

        if(f < opts->isoValues && f > -opts->isoValues)
        {
            P[i] = 1;
        }
    }

    return;
}

void saveTH_Geometry(meshInfo *mesh, char *P, std::string filename)
{
    FILE *SAVE = fopen(filename.c_str(), "w+");

    fprintf(SAVE, "x,y,z\n");

    for(int i = 0; i < mesh->nElements; i++)
    {
        if(P[i] == 1)
        {
            int slice = i / (mesh->numCellsX * mesh->numCellsY);
            int row = (i - slice * mesh->numCellsX * mesh->numCellsY)/mesh->numCellsX;
            int col = i - slice * mesh->numCellsX * mesh->numCellsY - row * mesh->numCellsX;
            fprintf(SAVE, "%d,%d,%d\n",col, row, slice);
        }
    }

    fclose(SAVE);

    return;
}

void fixSD_info(options *opts, meshInfo *mesh, geometry *geo, char *subDomain)
{
    /*
     * Add function description later if it works.
     *
     *
     * */
    memset(subDomain, 0, sizeof(char)*mesh->nElements);

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

        float f = TH_CaseF[opts->TPMS_Type - 1](x,y,z, (geo->D2 + geo->D1)/2.0);

        if(f < opts->isoValues && f > -opts->isoValues)
        {
            subDomain[i] = 0;
        } else if(f >= opts->isoValues)
        {
            subDomain[i] = 1;
        }else if(f <= -opts->isoValues)
        {
            subDomain[i] = 2;
        }
        else{
            printf("Error on subDomain correction!\n");
        }
    }

    return;
}

/*
 *
 *
 *
 * Debug Functions:
 *
 *
 *
 * */


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
 *
 *
 * Shape Factor Simulations:
 *
 *
 * */

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


float zeroFunctionTH(options *opts, Point p, float center)
{
    /*
        Function zeroFunctionTH:
        Inputs:
            - pointer to opts
            - point (x,y,z)
        Output:
            - returns r = C - |S_a(x,y,z)|, where 
            C is the isovalue and S_a is the TPMS
            function evaluated at px, py, and pz.
    */
    float result = opts->isoValues - fabs(TH_CaseF[opts->TPMS_Type - 1](p.x, p.y, p.z, center));

    return result;
}


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
    }

    return;
}

float bisectionIBM(options *opts, meshInfo *mesh, int col_idx, int row_idx, int slice_idx, int nb_col, int nb_row, int nb_slice, float center, int dir_search)
{
    /*
        Function bisectionIBM:
        Inputs:
            - pointer to opts
            - pointer to mesh
            - idx of current location
            - idx of neighbor
            - float center
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

    float p1_value = zeroFunctionTH(opts, p1, center);
    float p2_value = zeroFunctionTH(opts, p2, center);

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

    eps = fabs(zeroFunctionTH(opts, c, center));

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
        cVal = zeroFunctionTH(opts, c, center);
        // assign new bounds
        if(cVal > 0)
            b = c;
        else if(cVal < 0)
            a = c;
        
        // find new c

        c.x = (a.x + b.x)/2;
        c.y = (a.y + b.y)/2;
        c.z = (a.z + b.z)/2;

        eps = fabs(zeroFunctionTH(opts, c, center));

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


void correctIBs(float *DC, options *opts, meshInfo *mesh, geometry *geo, IBM_Correct *IB)
{
    /*
        correctIBs Function:
        Inputs:
            - pointer to DC
            - pointer to opts
            - pointer to meshInfo
            - geometry
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
    float center = geo->center;
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
                float dx_w = bisectionIBM(opts, mesh, col, row, slice, -1, row, slice, center, 0); // special periodic dist -1
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
                float dx_w = bisectionIBM(opts, mesh, col, row, slice, (col-1), row, slice, center, 0);
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
                float dx_e = bisectionIBM(opts, mesh, col, row, slice, nCols, row, slice, center, 0); // special periodic dist nCols
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
                float dx_e = bisectionIBM(opts, mesh, col, row, slice, (col+1), row, slice, center, 0);
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
                float dy_s = bisectionIBM(opts, mesh, col, row, slice, col, nRows, slice, center, 1); // special periodic dist nRows
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
                float dy_s = bisectionIBM(opts, mesh, col, row, slice, col, row + 1, slice, center, 1);
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
                float dy_n = bisectionIBM(opts, mesh, col, row, slice, col, -1, slice, center, 1); // special periodic dist -1
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
                float dy_n = bisectionIBM(opts, mesh, col, row, slice, col, row - 1, slice, center, 1);
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
                float dz_b = bisectionIBM(opts, mesh, col, row, slice, col, row, nSlices, center, 2); // special periodic dist nSlices
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
                float dz_b = bisectionIBM(opts, mesh, col, row, slice, col, row, slice + 1, center, 2);
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
                float dz_f = bisectionIBM(opts, mesh, col, row, slice, col, row, -1, center, 2); // special periodic dist -1
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
                float dz_f = bisectionIBM(opts, mesh, col, row, slice, col, row, slice - 1, center, 2);
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

        if (fabs(DC[i] - 1) > 0.5 )
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


int main()
{
    // declare data structure with general user options
    options opts;
    meshInfo mesh;
    saveInfo save;


    geometry geo;
    solutions_sf sol;

    // Initialize some stuff

    optionsInit(&opts);

    geo.D1 = 0.0;
    geo.D2 = 0.0;

    opts.CLeft = 0;
    opts.CRight = 1;

    opts.sfTMAP = 1;

    opts.MAX_ITER = 1e6;
    opts.ConvergeCriteria = 1e-6;

    printIndex(2);      // hardcoded number of theoretical cases

    // Choose the case and mesh size
 
    bool acceptableInput = false;

    unsigned int case_num;

    while(!acceptableInput)
    {
        printf("Please Enter the Case Number:\n");
        std::cin >> case_num;
        if(case_num == 1 || case_num == 2)
        {
            acceptableInput = true;
        }
        std::cin.clear();
    }

    // set TPMS_Type = case_num to facilitate things
    opts.TPMS_Type = case_num;

    acceptableInput = false;

    while(!acceptableInput)
    {
        printf("Please Enter Mesh Size:\n");
        std::cin >> mesh.numCellsX;
        if(mesh.numCellsX <= 25 || mesh.numCellsX >= 1000)
        {
            printf("Selected mesh size is not acceptable, try a number between 25 and 1000\n");
            std::cin.clear();
        }
        else
        {
            acceptableInput = true;
        }
    }

    // initialie everything that might be needed

    mesh.numCellsY = mesh.numCellsX;
    mesh.numCellsZ = mesh.numCellsX;

    mesh.dx = 2.0 * PI / mesh.numCellsX;
    mesh.dy = 2.0 * PI / mesh.numCellsY;
    mesh.dz = 2.0 * PI / mesh.numCellsZ;

    mesh.nElements = mesh.numCellsX * mesh.numCellsY * mesh.numCellsZ;

    char *P = (char *)malloc(sizeof(char) * mesh.nElements);
    memset(P, 0, sizeof(char) * mesh.nElements);

    float *DC = (float *)malloc(sizeof(float) * mesh.nElements);
    memset(DC, 0.0, sizeof(float) * mesh.nElements);

    // set omp parameters
    unsigned int nThreads = 1;

    acceptableInput = false;

    while(!acceptableInput)
    {
        printf("Enter Number of Threads:\n");
        std::cin >> nThreads;
        if(nThreads > 1)
        {
            acceptableInput = true;
        }
        else
        {
            acceptableInput = false;
            std::cin.clear();
        }
    }

    printf("Proceeding with %d threads\n", nThreads);

    opts.nThreads = nThreads;

    if(case_num == 1)
    {
        printf("------------------------------------------\n");
        printf("    Creating Square Flow Passage\n");
        printf("------------------------------------------\n");
        
        acceptableInput = false;

        while(!acceptableInput)
        {
            printf("Enter inner-side length (\"a\"):\n");
            std::cin >> geo.D1;
            if(geo.D1 < 0.2 || geo.D1 > 5.75)
            {
                printf("Error, try a value between 0.2 and 5.75\n");
                std::cin.clear();
            }
            else
            {
                acceptableInput = true;
            }
        }

        acceptableInput = false;
        while(!acceptableInput)
        {
            printf("Enter outer-side length (\"b\"):\n");
            std::cin >> geo.D2;
            if(geo.D2 < geo.D1 || geo.D2 > 5.75)
            {
                printf("Error, try a value between %1.2f and 5.75\n", geo.D1);
                std::cin.clear();
            }
            else
            {
                acceptableInput = true;
            }
        }

        // set isovalue
        opts.isoValues = (geo.D2 - geo.D1)/4.0;
        geo.center = (geo.D2 + geo.D1)/2.0;
        // generate structure

        generateFlowPassage(&opts, &mesh, &geo, P);

    }
    else if(case_num == 2)
    {
        printf("------------------------------------------\n");
        printf("    Creating Long Cylindrical Layer\n");
        printf("------------------------------------------\n");
    }

    // Save structure?
    acceptableInput = false;

    bool saveStruct = false;
    
    printf("Save Geometry? (0) No (1) Yes\n");
    std::cin >> saveStruct;

    std::cin.clear();

    // If user entered nonsense, then don't save

    std::string geometry_filename;

    if(saveStruct)
    {
        printf("Enter file name (without extension):\n");
        std::cin >> geometry_filename;
        std::string file_ext = geometry_filename + ".csv";
        std::cin.clear();

        saveTH_Geometry(&mesh, P, file_ext);
    }
    else
    {
        printf("Geometry won't be saved, moving on...\n");
    }

    // discretize and assign subDomains
    char *subDomain = (char *)malloc(sizeof(char) * mesh.nElements);
    memset(subDomain, 0, mesh.nElements);

    fixSD_info(&opts, &mesh, &geo, subDomain);

    // Set DC's based on the structure of subDomain

    SetDC_SF(DC, subDomain, &mesh);

    // allocate the arrays for simulation

    float *CoeffMatrix = (float *)malloc(mesh.nElements * 7 * sizeof(float));
    float *RHS = (float *)malloc(mesh.nElements * sizeof(float));
    float *Temperature = (float *)malloc(mesh.nElements * sizeof(float));

    // initialize memory

    memset(CoeffMatrix, 0, sizeof(float) * 7 * mesh.nElements);
    memset(RHS, 0, mesh.nElements * sizeof(float));

    // Initializing temperature distribution
    memset(Temperature, 0, mesh.nElements * sizeof(float));
    
    // allocate the struct

    IBM_Correct *IBs = (IBM_Correct *)malloc(sizeof(IBM_Correct) * mesh.nElements);

    memset(IBs, 0, sizeof(float)*mesh.nElements*6);

    correctIBs(DC, &opts, &mesh, &geo, IBs);

    // Discretize

    Disc3D_IB_SF(&opts, &mesh, DC, CoeffMatrix, RHS, IBs);

    // Initialize temperature again

    if(opts.verbose)
        printf("Initializing Temperature\n");

    TemperatureInit(Temperature, DC, &mesh, &opts);

    // CPU solve
    pGS3D_SF_handle(&opts, &mesh, &save, CoeffMatrix, RHS, Temperature);

    // Calculate Shape Factor

    float Q_avg_IBs = Calc_Q_SF_IBs(&opts, &mesh, Temperature, DC, IBs);

    save.SF = Q_avg_IBs/(opts.CLeft - opts.CRight);

    float TH_SF = 0.0;

    // hardcode theoretical cases because why not
    if(case_num == 1)
    {
        if(geo.D2/geo.D1 > 1.4)
        {
            TH_SF = 2*PI*(2*PI)/(0.93*std::log(0.948*geo.D2/geo.D1));
        }else
        {
            TH_SF = 2*PI*(2*PI)/(0.785*std::log(geo.D2/geo.D1));
        }
    }
 

    printf("Simulated: %1.3e, Theoretical: %1.3e\n", save.SF, TH_SF);


    // save Temp?
    if(opts.sfTMAP)
        saveTemp_SF(DC, Temperature, &mesh);


    // memory management
    free(Temperature);
    free(CoeffMatrix);
    free(RHS);


    return 0;
}
