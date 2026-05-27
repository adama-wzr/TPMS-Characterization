/*

Test function for printing the surface mesh
of a TPMS based on the Marching Cubes algorithm.

Last modified 05/25/2026
Andre Adam.
*/


#include <stdio.h>
#include <lib/TPMS_helpers.hpp>
#include <lib/TPMS_helpers.hpp>
#include <usrInput.hpp>
#include <lib/subDomainFF.hpp>
#include <lib/marchingCubes.hpp>

void printIndex(char max)
{
    printf("----------------------------------\n\n");
    printf("          TPMS Index:\n\n");
    printf("------------------------------------\n\n");
    printf("Index      Name\n");
    printf("_____      ____\n");

    for(int i = 0; i < max; i++)
    {
        printf("%02d         %s\n", i+1, TPMS_Names[i]);
    }

    printf("------------------------------------\n\n");

    return;
}

int main(int argc, char **argv)
{
    // Declare data structure with general user options
    options opts;
    meshInfo mesh;

    // Index max and min (beautifully hardcoded)
    int min, max;
    min = 1;
    max = 27;

    printIndex(max);

    // User entered index
    bool acceptableInput = false;
    int input;
    while (!acceptableInput)
    {
        printf("Enter TPMS Index:\n");
        std::cin >> input;

        if (input >= min && input <= max)
            acceptableInput = true;
    }

    // user entered size

    int nVoxel = 0;
    acceptableInput = false;

    while (!acceptableInput)
    {
        printf("Enter Side Length in Pixels:\n");
        std::cin >> nVoxel;

        if (nVoxel < 25)
        {
            printf("Too small, enter number >= 25.\n");
        }
        else if (nVoxel >= 1000)
        {
            printf("Too large, try a number < 1000.\n");
        }
        else
        {
            acceptableInput = true;
        }
    }

    // user entered bounds

    acceptableInput = false;
    float iso = 0.0;
    while(!acceptableInput)
    {
        printf("Enter isovalue:\n");
        std::cin >> iso;
        if (iso > TPMS_Crit[input - 1])
        {
            printf("Isovalue entered %2.1f is larger than crit value %2.1f. Please try again.\n", iso, TPMS_Crit[input - 1]);
        }
        else
        {
            acceptableInput = true;
        }
    }

    // warn if iso > pinch
    if (iso > TPMS_Pinch[input - 1])
    {
        printf("********************\n\n");
        printf("*     WARNING!     *\n");
        printf("Isovalue entered %2.1f is greater than the pinch value %2.1f.\n", iso, TPMS_Pinch[input - 1]);
        printf("********************\n\n");
    }
    
    // set structs
    opts.nVoxels = nVoxel;
    opts.TPMS_Type = input;
    opts.isoValues = iso;

    // Generate TPMS

    char *P;

    TPMS_Init(&P, &opts, &mesh);

    // marching cubes

    std::vector<std::vector<Point>> triangles;

    triangles = MarchingCubes(&opts, &mesh);

    // User entered index
    acceptableInput = false;
    while (!acceptableInput)
    {
        printf("Save as csv (0) or ply (1)?\n");
        std::cin >> input;
        if (input == 0 || input == 1)
            acceptableInput = true;
    }

    if(input == 0)
    {
        FILE *OUT = fopen("triangles.csv", "w+");

        fprintf(OUT, "x,y,z\n");

        for(int i = 0; i < (int)triangles.size(); i++)
        {
            for(int j = 0; j < 3; j++){
                fprintf(OUT, "%1.2e,%1.2e,%1.2e\n", triangles[i][j].x, triangles[i][j].y, triangles[i][j].z);
            }
        }

        fclose(OUT);
    }
    else
    {
        // Save triangles to ply file
        std::string filename = "triangles.ply";
        std::cout << "Number of triangles: " << (int)triangles.size() << "\n";
        write_to_ply(triangles, filename.c_str());
    }

    return 0;
}