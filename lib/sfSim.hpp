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
#include <Disc3D_SF_PB.hpp>

#ifdef USE_CUDA
    #include <cuda_solvers/gpuSolve.cu>
#endif

#ifndef USE_CUDA
    #include <cpu_solvers/cpuErrorHandler.hpp>
#endif

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
    memset(Concentration, 0, mesh->nElements * sizeof(float));

    // Linear initialize the concentration (not needed?)

    // for (long int i = 0; i < mesh->nElements; i++)
    // {
    //     if (DC[i] == 0)
    //         continue;
    //     int slice = i / (mesh->numCellsX * mesh->numCellsY);
    //     int row = (i - slice * mesh->numCellsX * mesh->numCellsY)/mesh->numCellsX;
    //     int col = i - slice * mesh->numCellsX * mesh->numCellsY - row * mesh->numCellsX;
    //     Concentration[i] = ((float)col / mesh->numCellsX) * (opts->CRight - opts->CLeft) + opts->CLeft;
    // }

    // Discretize

    if(opts->verbose)
        printf("Discretizing\n");

    Disc3D_SF_PB(opts, mesh, DC, CoeffMatrix, RHS);

    // Solve

    bool errorFlag = 0;

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

    return 1;

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

    int POI = 1;

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
