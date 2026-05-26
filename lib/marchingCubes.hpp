#include <data_structures.hpp>
#include <constants.hpp>
#include <TPMS_definitions.hpp>

Point LinearInterpolate(Point one, Point two)
{
    /*
        Function to reduce the size for the interpolation step.

        - This function should be temporary.
    */

    Point interpolated_result;

    /*
    
        Code to interpolate.

    */



    return interpolated_result;
}

void MarchingCubes(char *P, options *opts, meshInfo *mesh, MarchingCubesTriangles *MCT)
{
    /*
        Function Marching Cubes:
        Inputs:
            - pointer to TPMS
            - pointer to usr options
            - pointer to mesh struct
            - pointer to MarchingCubesTriangles struct
        Outputs:
            - None
        
        Function will get TPMS structure as input, calculate the surface mesh using
        the marching cubes method. Triangle locations are stored in the MCT
        struct.
    */

    long int nTriangles = 0;

    // allocate MCT struct

    MCT = (MarchingCubesTriangles *)malloc(sizeof(MarchingCubesTriangles)*mesh->nElements);

    int row, col, slice;

    // iterate over the whole image
    for(int i = 0; i < mesh->nElements; i++)
    {
        // Break down i into the integer indexes
        int slice = i / (mesh->numCellsX * mesh->numCellsX);
        int row = (i - slice * mesh->numCellsX * mesh->numCellsX) / mesh->numCellsX;
        int col = i - slice * mesh->numCellsX * mesh->numCellsX - row * mesh->numCellsX;

        int cubeIndex = 0;

        // Find corners inside of the TPMS structure, update binary array

        for(int j = 0; j < 8; j++)
        {
            /*
                Index:
                0 = {x      , y     , z + dz}
                1 = {x + dx , y     , z + dz}
                2 = {x + dx , y     , z     }
                3 = {x      , y     , z     }
                4 = {x      , y + dy, z + dz}
                5 = {x + dx , y + dy, z + dz}
                6 = {x + dx , y + dy, z     }
                7 = {x      , y + dy, z     }
            */

            // 8 corners in a cube
            if(j == 0)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  1;
            }
            else if(j == 1)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;
                
                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  2;
            }
            else if(j == 2)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  4;
            }
            else if(j == 3)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  8;
            }
            else if(j == 4)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  16;
            }
            else if(j == 5)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  32;
            }
            else if(j == 6)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues)) 
                    cubeIndex |=  64;
            }
            else if(j == 7)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues)) 
                    cubeIndex |=  128;
            }
        } // endfor cubeIdx

        // if not interface skip
        if(cubeIndex == 0 || cubeIndex == 255)
            continue;

        // Lookup edge (indexes)

        int edgeKey = edgeTable[cubeIndex];

        // find edge intercepts
        int intercept[16][3];
        int idx = 0;
        
        while(edgeKey)
        {
            if(edgeKey&1)
            {
                int v1 = edgeToVertices[idx].first;
                int v2 = edgeToVertices[idx].second;

                // interpolate

            }
        }

        for(int j = 0; j < 16; j++)
        {
            intercept[j][0] = -1;   // x
            intercept[j][1] = -1;   // y
            intercept[j][2] = -1;   // z
        }

        // Populate edge intercepts w/ correct index





    } // endfor

    return;
}