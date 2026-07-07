/*

Test function for printing the surface mesh
of a TPMS based on the Marching Cubes algorithm.

Last modified 05/25/2026
Andre Adam.
*/


#include <stdio.h>
#include <lib/TPMS_helpers.hpp>
#include <usrInput.hpp>
#include <lib/subDomainFF.hpp>
#include <lib/marchingCubes.hpp>
#include <lib/surfaceArea.hpp>
#include <string>

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
    min = 0;
    max = 28;

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
        else
            std::cin.clear();
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
            std::cin.clear();
        }
        else if (nVoxel >= 1000)
        {
            std::cin.clear();
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
            std::cin.clear();
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
        printf("Save as csv (0) or ply (1)? Skip (2)...\n");
        std::cin >> input;
        if (input == 0 || input == 1 || input == 2)
            acceptableInput = true;
        else
            std::cin.clear();
    }

    std::string filename;

    if(input == 0)
    {
        printf("Enter file name (without extension):\n");
        std::cin >> filename;

        std::string file_ext = filename + ".csv";

        FILE *OUT = fopen(file_ext.c_str(), "w+");

        fprintf(OUT, "x,y,z\n");

        for(int i = 0; i < (int)triangles.size(); i++)
        {
            for(int j = 0; j < 3; j++){
                fprintf(OUT, "%1.2e,%1.2e,%1.2e\n", triangles[i][j].x, triangles[i][j].y, triangles[i][j].z);
            }
        }

        fclose(OUT);
    }
    else if(input == 1)
    {
        printf("Enter file name (withiut extension):\n");
        // Save triangles to ply file
        std::cin >> filename;
        std::string file_ext = filename + ".ply";
        std::cout << "Number of triangles: " << (int)triangles.size() << "\n";
        write_to_ply(triangles, file_ext.c_str());
    }
    else
    {
        // not saving anything bye
        printf("Not saving, option selected was not recognized\n");
    }

    // Calculate Surface Area?

    acceptableInput = false;

    while (!acceptableInput)
    {
        printf("Surface Area? (0) No (1) Yes\n");
        std::cin >> input;
        if (input == 0 || input == 1)
            acceptableInput = true;
        else
            std::cin.clear();
    }

    if(input == 1)
    {
        saveInfo save;
        SA_Triangles(&opts, triangles, &save);
    }



    return 0;
}
