

#ifndef __OUT
#define __OUT

#include <data_structures.hpp>
#include <stdlib.h>


void printCMAP(options *opts, meshInfo *mesh, char *filename, float *C, char *P, int POI)
{
    /*
        printCMAP Function:
        Inputs:
            - pointer to user options
            - pointer to mesh info
            - pointer to filename (to save)
            - pointer to concentration distribution
            - pointer to phase-array
            - phase of interest for CMAP
        Outputs:
            - none.
        
        Function will create filename
    */
    // create file
    FILE *OUT = fopen(filename, "w+");

    fprintf(OUT,"x,y,z,C\n");

    int row, col, slice;

    // print

    for(int i = 0; i < mesh->nElements; i++)
    {
        if(P[i] != POI)
            continue;
        
        slice = i / (mesh->numCellsX * mesh->numCellsY);
        row = (i - slice * mesh->numCellsX*mesh->numCellsY)/mesh->numCellsX;
        col = i - slice * mesh->numCellsX*mesh->numCellsY - row * mesh->numCellsX;
        
        fprintf(OUT,"%d,%d,%d,%1.3e\n", col, row, slice, C[i]);
    }

    // close

    fclose(OUT);

    return;
}


void outputGeneral(options *opts, saveInfo *save, meshInfo *mesh)
{
    /*
        outputGeneral Function:
        Inputs:
            - pointer to user options struct
            - pointer to saveInfo struct
            - pointer to meshInfo struct
        Output:
            - None

        Function will print the output from running the code to a output file
        taking in consideration the user options. Handles general output.
    */

    bool headerFlag = true;

    // Check if file exists

    if (FILE *TEST = fopen(opts->outputFilename, "r"))
    {
        fclose(TEST);
        headerFlag = false;
    }

    // Open file

    FILE *OUT = fopen(opts->outputFilename, "a+");

    if (headerFlag && opts->subOut)
    {
        fprintf(OUT, "Structure,nElement,iso,pore,SVF,SA,"
                     "ePore,Tau,TauSolid,pore50,part50");
        
        int nChannel = 0;
        for(int nSub = 1; nSub <= mesh->nChannels; nSub++)
        {
            if(mesh->sdInfo[nSub - 1].FC)
            {
                nChannel++;
                fprintf(OUT, ",Tau%d,pore50_%d,SA%d,eVF%d", nChannel, nChannel, nChannel, nChannel);
            }
            else
                continue;
        }
        fprintf(OUT, "\n");
    }
    else if (headerFlag)
    {
        fprintf(OUT, "Structure,nElement,iso,pore,SVF,SA,"
                     "ePore,Tau,TauSolid,pore50,part50,SF\n");
    }

    // print output from inputs

    fprintf(OUT, "%d,%d,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e",
        opts->TPMS_Type, save->nVoxel, opts->isoValues,
        save->porosity, save->SVF, save->SA, save->ePore,
        save->Tau, save->TauSolid, save->pore50, save->part50, save->SF);
    
    if(opts->subOut)
    {
        for(int nSub = 1; nSub <= mesh->nChannels; nSub++)
        {
            if(mesh->sdInfo[nSub - 1].FC)
            {
                fprintf(OUT, ",%1.3e,%1.3e,%1.3e,%1.3e",
                    mesh->sdInfo[nSub - 1].Tau, mesh->sdInfo[nSub - 1].pore50,
                    mesh->sdInfo[nSub - 1].SA, mesh->sdInfo[nSub - 1].VF);
            }
        }
        fprintf(OUT,"\n");
    }
    else
    {
        fprintf(OUT,"\n");
    }

    fclose(OUT);
    return;
}

#endif