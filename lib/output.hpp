

#ifndef __OUT
#define __OUT

#include <data_structures.hpp>
#include <stdlib.h>

void outputGeneral(options *opts, saveInfo *save)
{
    /*
        outputGeneral Function:
        Inputs:
            - pointer to user options struct
            - pointer to saveInfo struct
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

    if (headerFlag)
    {
        fprintf(OUT, "Structure,nElement,iso,pore,SA,"
                     "DeffMax,Deff,Tau,pore50,part50\n");
    }

    // print output from inputs

    fprintf(OUT, "%d,%d,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e\n",
        opts->TPMS_Type, save->nVoxel, opts->isoValues,
        save->porosity, save->SA, save->Deff_TH_MAX,
        save->Deff, save->Tau, save->pore50, save->part50);

    fclose(OUT);
    return;
}

#endif