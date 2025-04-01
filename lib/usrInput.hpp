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
    opts->poreSD = 0;
    opts->partSD = 0;
    opts->CLeft = 0;
    opts->CRight = 1;

    // Default Solver Options
    opts->nThreads = 1;
    opts->nGPU = 1;     // later need to add a GPU flag

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
            opts->MAX_ITER = (long int) tempD;
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
    }

    return;
}