/*

Main file will handle program execution.

Last modified 04/01/2025
Andre Adam
*/

#include <lib/TPMS_definitions.hpp>
#include <lib/data_structures.hpp>
#include <lib/usrInput.hpp>
#include <lib/TPMS_helpers.hpp>
#include <math.h>


int main(int argc, char **argv)
{
    // Declare data structure with general user options
    options opts;
    meshInfo mesh;
    saveInfo save;

    // Read user input

    char inputFilename[40];

    sprintf(inputFilename, "input.txt");

    readInputGeneral(inputFilename, &opts);

    if (opts.verbose)
        printOptsGeneral(&opts);

    // Create TPMS

    char *P;

    TPMS_Init(P, &opts, &mesh);


    
    return 0;
}