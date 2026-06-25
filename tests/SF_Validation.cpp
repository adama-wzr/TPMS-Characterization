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
#include <lib/data_structures.hpp> 
#include <usrInput.hpp>
#include <lib/sfSim.hpp>

typedef struct{
    float SF_TH;
    float SF_sim;
    float area_MC;
    float area_Voxel;
} th_sol;

typedef struct{
    float L;
    float D1;       // D1 is "a" for square flow passage
    float D2;       // D2 is "b" for square flow passage
} geometry;


int main()
{
    // declare data structure with general user options
    options opts;
    meshInfo mesh;

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
    mesh.numCellsZ = mesh.numCellsZ;

    mesh.dx = 2.0 * PI / mesh.numCellsX;
    mesh.dy = 2.0 * PI / mesh.numCellsY;
    mesh.dz = 2.0 * PI / mesh.numCellsZ;

    mesh.nElements = mesh.numCellsX * mesh.numCellsY * mesh.numCellsZ;

    char *P = (char *)malloc(sizeof(char) * mesh.nElements);
    memset(P, 0, sizeof(char) * mesh.nElements);

    float *DC = (float *)malloc(sizeof(float) * mesh.nElements);
    memset(DC, 0.0, sizeof(float) * mesh.nElements);

    if(case_num == 1)
    {
        printf("------------------------------------------\n");
        printf("    Creating Square Flow Passage\n");
        printf("------------------------------------------\n");
        
    }
    else if(case_num == 2)
    {
        printf("------------------------------------------\n");
        printf("    Creating Long Cylindrical Layer\n");
        printf("------------------------------------------\n");
    }

    return 0;
}
