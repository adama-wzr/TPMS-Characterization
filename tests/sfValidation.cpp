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
#include <math.h>
#include <lib/data_structures.hpp> 
#include <usrInput.hpp>

/*
 *
 * Special Data Structures
 *
 * */

typedef struct{
    float SF_TH;
    float SF_sim;
    float area_MC;
    float area_Voxel;
} solutions_sf;

typedef struct{
    float L;
    float D1;       // D1 is "a" for square flow passage
    float D2;       // D2 is "b" for square flow passage
} geometry;

typedef float (*th_case_f_ptr)(float, float, float, float);

/*
 *
 *      Implicit Functions:
 *
 * */

float sqPassage_F(float x, float y, float z, float w)
{
    float f = fabs(std::max(fabs(x),fabs(y)) - w/2.0);
    return f;
}

float cylTube_F(float x, float y, float z, float r)
{
    float f = x + y + z + r;
    return f;
}

/*
 *
 *      Lookup Tables
 *
 * */

static const char *TH_Cases[] =
{
    "Square Flow Passage", "Long Cylindrical Layer"
};

th_case_f_ptr TH_CaseF[] = {
    sqPassage_F,
    cylTube_F
};

/*
 *
 * Helper Functions (save/opts/print/generate)
 *
 * */

void printIndex(char max)
{
    printf("----------------------------------\n\n");
    printf("    Theoreical Case Index:\n\n");
    printf("------------------------------------\n\n");
    printf("Index      Name\n");
    printf("_____      ____\n");

    for(int i = 0; i < max; i++)
    {
        printf("%02d         %s\n", i+1, TH_Cases[i]);
    }

    printf("------------------------------------\n\n");

    return;
}

void generateFlowPassage(options *opts, meshInfo *mesh, geometry *geo, char *P)
{
    /*
     *  Function generateFlowPassage:
     *  
     *  Inputs:
     *      - pointer to options opts;
     *      - pointer to meshInfo mesh;
     *      - pointer to geometry geo;
     *      - pointer to char array holding structure;
     *  Outputs;
     *      - none
     *
     *  Function will create the flow passage according to the implicit equation
     *  and the user-entered information.
     *
     * */

    // define some stuff to make life easier
    float dx = mesh->dx;
    float dy = mesh->dy;
    float dz = mesh->dz;
    
    int nCols = mesh->numCellsX;
    int nRows = mesh->numCellsY;
    int nSlices = mesh->numCellsZ;

    float x, y, z;
    int col, row, slice;

    for(int i = 0; i < mesh->nElements; i++)
    {
        // get col, row, slice
        slice = i / (nCols * nRows);
        row = (i - slice * nCols * nRows) / nCols;
        col = i - slice * nRows * nCols - row * nCols;

        x = -PI + (float) col * dx + dx / 2.0;
        y = -PI + (float) row * dy + dy / 2.0;
        z = -PI + (float) slice * dz + dz / 2.0;

        float f = TH_CaseF[opts->TPMS_Type - 1](x,y,z, (geo->D2 + geo->D1)/2.0);

        if(f < opts->isoValues && f > -opts->isoValues)
        {
            P[i] = 1;
        }
    }

    return;
}

void saveTH_Geometry(meshInfo *mesh, char *P, std::string filename)
{
    FILE *SAVE = fopen(filename.c_str(), "w+");

    fprintf(SAVE, "x,y,z\n");

    for(int i = 0; i < mesh->nElements; i++)
    {
        if(P[i] == 1)
        {
            int slice = i / (mesh->numCellsX * mesh->numCellsY);
            int row = (i - slice * mesh->numCellsX * mesh->numCellsY)/mesh->numCellsX;
            int col = i - slice * mesh->numCellsX * mesh->numCellsY - row * mesh->numCellsX;
            fprintf(SAVE, "%d,%d,%d\n",col, row, slice);
        }
    }

    fclose(SAVE);

    return;
}

int main()
{
    // declare data structure with general user options
    options opts;
    meshInfo mesh;

    geometry geo;
    solutions_sf sol;

    // Initialize some stuff
    geo.D1 = 0.0;
    geo.D2 = 0.0;

    printIndex(2);      // hardcoded number of theoretical cases

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

    // set TPMS_Type = case_num to facilitate things
    opts.TPMS_Type = case_num;

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
    mesh.numCellsZ = mesh.numCellsX;

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
        
        acceptableInput = false;

        while(!acceptableInput)
        {
            printf("Enter inner-side length (\"a\"):\n");
            std::cin >> geo.D1;
            if(geo.D1 < 0.2 || geo.D1 > 5.75)
            {
                printf("Error, try a value between 0.2 and 5.75\n");
                std::cin.clear();
            }
            else
            {
                acceptableInput = true;
            }
        }

        acceptableInput = false;
        while(!acceptableInput)
        {
            printf("Enter outer-side length (\"b\"):\n");
            std::cin >> geo.D2;
            if(geo.D2 < geo.D1 || geo.D2 > 5.75)
            {
                printf("Error, try a value between %1.2f and 5.75\n", geo.D1);
                std::cin.clear();
            }
            else
            {
                acceptableInput = true;
            }
        }

        // set isovalue
        opts.isoValues = (geo.D2 - geo.D1)/2;

        // generate structure

        generateFlowPassage(&opts, &mesh, &geo, P);

    }
    else if(case_num == 2)
    {
        printf("------------------------------------------\n");
        printf("    Creating Long Cylindrical Layer\n");
        printf("------------------------------------------\n");
    }

    // Save structure?
    acceptableInput = false;

    bool saveStruct = false;
    
    printf("Save Geometry? (0) No (1) Yes\n");
    std::cin >> saveStruct;

    std::cin.clear();

    // If user entered nonsense, then don't save

    std::string geometry_filename;

    if(saveStruct)
    {
        printf("Enter file name (without extension):\n");
        std::cin >> geometry_filename;
        std::string file_ext = geometry_filename + ".csv";
        std::cin.clear();

        saveTH_Geometry(&mesh, P, file_ext);
    }
    else
    {
        printf("Geometry won't be saved, moving on...\n");
    }

    return 0;
}
