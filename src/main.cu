#include "helper.cuh"

int main(int argc, char **argv)
{
    // Print efficiency for Linux
    fflush(stdout);

    // Define the data structures
    options opts;
    meshInfo info;
    saveInfo save;
    // Read filename
    char inputFilename[50];

    sprintf(inputFilename, "input.txt");

    readInputGeneral(inputFilename, &opts);
    
    if (opts.verbose)
        printOpts(&opts);

    // Input file read, fill data strcture
    info.numCellsX = opts.nVoxels;
    info.numCellsY = opts.nVoxels;
    info.numCellsZ = opts.nVoxels;
    info.nElements = opts.nVoxels * opts.nVoxels * opts.nVoxels;

    save.nElements = info.nElements;
    save.nVoxel = opts.nVoxels;

    // Open space to hold the structure and initialize it to zeroes
    char *P = (char *)malloc(sizeof(char) * info.nElements);
    memset(P, 0, sizeof(char) * info.nElements);

    // Generate structure according to user input
    int genStructFlag = 0;
    genStructFlag = genStruct(&opts, &info, P);

    if (genStructFlag == 1)
    {
        printf("Error with TPMS type selected. Exiting now.\n");
        return 1;
    }

    // Copy VF

    save.porosity = 1.0 - info.VF;
    
    // Run surface area
    
    SA(P, &info, &save);
    
    // save SA

    // outSA(&opts, &save);

    // tau Sim

    if (opts.Tau == 1)
    {
        TauSim3D(&opts, &info, &save, P);
    }

    // Pore and Particle size distributions

    if(opts.partSD)
        partSD_3D(&opts, &info, &save, P, 1);
    else
        save.part50 = 0;
    
    if (opts.poreSD)
        poreSD_3D(&opts, &info, &save, P, 0);
    else
        save.pore50 = 0;

    // Save the output
    
    outputGeneral(&opts, &save);

    return 0;
}