/*

User Input-parsing

Last modified 04/01/2025
Andre Adam

*/

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdbool.h>
#include <data_structures.hpp>

/*

Initialize opts

*/

void optionsInit(options *opts)
{
    /*
        optionsInit:
        Inputs:
            - pointer to options struct
    */

    // Initialize all flags to false

    opts->Tau = 0;
    opts->PB = 0;
    opts->runSA = 0;
    opts->poreSD = 0;
    opts->partSD = 0;
    opts->CLeft = 0;
    opts->CRight = 1;

    // Default Solver Options
    opts->nThreads = 1;
    opts->useGPU = 0;
    opts->nGPU = 1;

    // Other Defauls
    opts->maxR = 1000;

    // default TPMS_Type
    opts->TPMS_Type = 1;

    // Allocate Space for file names

    opts->outputFilename = (char *)malloc(200 * sizeof(char));
    opts->partSDOut = (char *)malloc(200 * sizeof(char));
    opts->poreSDOut = (char *)malloc(200 * sizeof(char));

    return;
}

/*

Un-Initialized Options:

*/

void unInitOptions(options *opts)
{
    /*
        unInitOptions:
        Inputs:
            - pointer to struct options
        Outputs:
            - None.
        Function will free the memory that was allocated by
        the options struct.
    */

    free(opts->outputFilename);
    free(opts->partSDOut);
    free(opts->poreSDOut);

    return;
}

/*

Print General options:

*/

void printOptsGeneral(options *opts)
{
    /*
        printOptsGeneral:
        Inputs:
            - pointer to usr inputs struct
        Outputs:
            - None

        Function simply prints the user entered options to the command line.
    */

    printf("----------------------------------\n\n");
    printf("         Selected Options:        \n\n");
    printf("-----------------------------------\n");

    printf("TPMS Type = %d\n", opts->TPMS_Type);
    printf("IsoValues = %1.3f\n", opts->isoValues);
    printf("Side Length = %d\n", opts->nVoxels);
    printf("Number of CPU Threads = %d\n", opts->nThreads);
    printf("Output Name = %s\n", opts->outputFilename);

    if (opts->Tau)
    {
        printf("--------------------------------\n");
        printf("Tortuosity Simulation Enabled\n");
        printf("GPU Flag = %d\n", opts->useGPU);
        if (opts->useGPU)
            printf("GPUs Requested = %d\n", opts->nGPU);

        if (opts->PB)
            printf("BC: Periodic\n");
        else
            printf("BC: No Flux\n");

        printf("Solver Options:\n");
        printf("SOR Jacobi Method\n"); // only option avaialable right now
        printf("Maximum Iterations = %ld\n", opts->MAX_ITER);
        printf("Convergence Type: Percent Change\n"); // also only options available
        printf("Convergence Cutoff = %1.3e\n", opts->ConvergeCriteria);
    }

    if (opts->poreSD)
    {
        printf("--------------------------------\n");
        printf("Pore-Size Distribution Enabled\n");
        printf("Output Pore-Size Distribution: %s\n", opts->poreSDOut);
        printf("Cutoff Radius = %d\n", opts->maxR);
    }

    if (opts->partSD)
    {
        printf("--------------------------------\n");
        printf("Particle-Size Distribution Enabled\n");
        printf("Output Particle-Size Distribution: %s\n", opts->partSDOut);
        printf("Cutoff Radius = %d\n", opts->maxR);
    }

    if (opts->runSA)
        printf("Surface Area Calculation Enabled\n");

    printf("----------------------------------\n\n");

    return;
}

/*

    Read Input General Case:

*/

void readInputGeneral(char *filename, options *opts)
{
    /*
        readInputGeneral Function:
        Inputs:
            - FileName: pointer to where the input file name is stored.
            - struct options: pass a struct with the options.
        Outputs: None

        Function reads the input file and stores the options in the opts struct.
    */

    // initiate necessary variables for input reading
    std::string myText;

    char tempC[1000];
    double tempD;
    char tempFilenames[1000];
    std::ifstream InputFile(filename);

    optionsInit(opts);

    // read options

    while (std::getline(InputFile, myText))
    {
        sscanf(myText.c_str(), "%s %lf", tempC, &tempD);
        if (strcmp(tempC, "OutputName:") == 0)
        {
            sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
            strcpy(opts->outputFilename, tempFilenames);
        }
        else if (strcmp(tempC, "Verbose:") == 0)
        {
            opts->verbose = (int)tempD;
        }
        else if (strcmp(tempC, "nVoxels:") == 0)
        {
            opts->nVoxels = (int)tempD;
        }
        else if (strcmp(tempC, "isoValue:") == 0)
        {
            opts->isoValues = (float)tempD;
        }
        else if (strcmp(tempC, "TPMS:") == 0)
        {
            opts->TPMS_Type = (int)tempD;
        }
        else if (strcmp(tempC, "Tau:") == 0)
        {
            opts->Tau = (int)tempD;
        }
        else if (strcmp(tempC, "poreSD:") == 0)
        {
            opts->poreSD = (int)tempD;
        }
        else if (strcmp(tempC, "partSD:") == 0)
        {
            opts->partSD = (int)tempD;
        }
        else if (strcmp(tempC, "pb:") == 0)
        {
            opts->PB = (int)tempD;
        }
        else if (strcmp(tempC, "nThreads:") == 0)
        {
            opts->nThreads = (int)tempD;
        }
        else if (strcmp(tempC, "Convergence:") == 0)
        {
            opts->ConvergeCriteria = tempD;
        }
        else if (strcmp(tempC, "MaxIter:") == 0)
        {
            opts->MAX_ITER = (long int)tempD;
        }
        else if (strcmp(tempC, "partSDOut:") == 0)
        {
            sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
            strcpy(opts->partSDOut, tempFilenames);
        }
        else if (strcmp(tempC, "poreSDOut:") == 0)
        {
            sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
            strcpy(opts->poreSDOut, tempFilenames);
        }
        else if (strcmp(tempC, "useGPU:") == 0)
        {
            opts->useGPU = (int)tempD;
        }
        else if (strcmp(tempC, "numGPU:") == 0)
        {
            opts->nGPU = (int)tempD;
        }
        else if (strcmp(tempC, "runSA:") == 0)
        {
            opts->runSA = (int)tempD;
        }
    }
    return;
}