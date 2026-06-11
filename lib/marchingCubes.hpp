
#ifndef _MC
#define _MC

#include <data_structures.hpp>
#include <constants.hpp>
#include <TPMS_definitions.hpp>

Point temp_optMC(options *opts, meshInfo *mesh, Point p1, Point p2)
{
    /*
    
        This is just a placeholder function to test the code out.

        if dir = 0, x
        if dir = 1, y
        if dir = 2, z
    
    */

    Point interceptPoint;

    // initialize things

    interceptPoint.x = -1;
    interceptPoint.y = -1;
    interceptPoint.z = -1;

    float eps = 1.000;
    float ib_delta = 0.0;
    
    int steps = 100;

    float step_size = 0;

    int lb, ub;

    // determine the direction

    float dir_search = -1;

    if(fabs(p1.x - p2.x) > mesh->dx/2)
    {
        dir_search = 0;
    }
    if(fabs(p1.y - p2.y) > mesh->dy/2)
    {
        if (dir_search != -1)
            printf("Inaccurate vertex information?\n");
        dir_search = 1;
    }
    if(fabs(p1.z - p2.z) > mesh->dz/2)
    {
        if (dir_search != -1)
            printf("Inaccurate vertex information?\n");
        dir_search = 2;
    }

    if(dir_search == 0)
    {
        // other intercepts and step size
        interceptPoint.y = p1.y;
        interceptPoint.z = p1.z;

        step_size = fabs(p1.x - p2.x)/steps;

        // find upper and lower bounds

        if(p1.x < p2.x)
        {
            ub = p2.x;
            lb = p1.x;
        }
        else if(p1.x > p2.x)
        {
            ub = p1.x;
            lb = p2.x;
        }
        else
        {
            printf("Big mistake, deltaX not calculated properly\n");
            return interceptPoint;
        }

        // find provisory eps and intercept

        eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p1.x, p1.y, p1.z));

        if(eps > fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p2.x, p2.y, p2.z)))
        {
            eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p2.x, p2.y, p2.z));
            interceptPoint.x = p2.x;
        }
        else
        {
            interceptPoint.x = p1.x;
        }

        // iterate over other values

        for(float i = lb; i < ub; i+=step_size)
        {

            float new_eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](i, p1.y, p1.z));
            
            if(new_eps < eps)
            {
                eps = new_eps;
                interceptPoint.x = i;
            }
        }
    }
    else if(dir_search == 1)
    {
        // set other intercepts and step size
        interceptPoint.x = p1.x;
        interceptPoint.z = p1.z;
        
        step_size = fabs(p1.y - p2.y)/steps;

        // upper and lower bounds

        if(p1.y < p2.y)
        {
            ub = p2.y;
            lb = p1.y;
        }
        else if(p1.y > p2.y)
        {
            ub = p1.y;
            lb = p2.y;
        }
        else
        {
            printf("Big mistake, deltaX not calculated properly\n");
            return interceptPoint;
        }

        // find provisory eps and intercept

        eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p1.x, p1.y, p1.z));

        if(eps > fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p2.x, p2.y, p2.z)))
        {
            eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p2.x, p2.y, p2.z));
            interceptPoint.y = p2.y;
        }
        else
        {
            interceptPoint.y = p1.y;
        }

        // iterate

        for(float i = lb; i < ub; i+=step_size)
        {

            float new_eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p1.x, i, p1.z));

            if(new_eps < eps)
            {
                eps = new_eps;
                interceptPoint.y = i;
            }
        }
    }
    else if(dir_search == 2)
    {
        // set other intercepts
        interceptPoint.y = p1.y;
        interceptPoint.x = p1.x;

        step_size = fabs(p1.z - p2.z)/steps;

        // upper and lower bounds

        if(p1.z < p2.z)
        {
            ub = p2.z;
            lb = p1.z;
        }
        else if(p1.z > p2.z)
        {
            ub = p1.z;
            lb = p2.z;
        }
        else
        {
            printf("Big mistake, deltaX not calculated properly\n");
            return interceptPoint;
        }

        // find provisory eps and intercept

        eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p1.x, p1.y, p1.z));

        if(eps > fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p2.x, p2.y, p2.z)))
        {
            eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p2.x, p2.y, p2.z));
            interceptPoint.z = p2.z;
        }
        else
        {
            interceptPoint.z = p1.z;
        }

        // iterate

        for(float i = lb; i < ub; i+=step_size)
        {

            float new_eps = fabs(opts->isoValues - TPMS_F[opts->TPMS_Type - 1](p1.x, p1.y, i));
            
            if(new_eps < eps)
            {
                eps = new_eps;
                interceptPoint.z = i;
            }
        }
    }

    // check results

    // if(dir_search == 0)
    // {
    //     if(p1.x > interceptPoint.x && p2.x > interceptPoint.x)
    //     {
    //         printf("Both p1.x and p2.x are greater than intercept!\n");
    //     } else if(p1.x < interceptPoint.x && p2.x < interceptPoint.x)
    //     {
    //         printf("Both p1.x and p2.x are smaller than intercept!\n");
    //     }
    // }else if(dir_search == 1)
    // {
    //     if(p1.y > interceptPoint.y && p2.y > interceptPoint.y)
    //     {
    //         printf("Both p1.y and p2.y are greater than intercept!\n");
    //     } else if(p1.y < interceptPoint.y && p2.y < interceptPoint.y)
    //     {
    //         printf("Both p1.y and p2.y are smaller than intercept!\n");
    //     }
    // }else if(dir_search == 2)
    // {
    //     if(p1.z > interceptPoint.z && p2.z > interceptPoint.z)
    //     {
    //         printf("Both p1.z and p2.z are greater than intercept!\n");
    //     } else if(p1.z < interceptPoint.z && p2.z < interceptPoint.z)
    //     {
    //         printf("Both p1.z and p2.z are smaller than intercept!\n");
    //     }
    // }
    
    return interceptPoint;
}

VertexContainer hash_vertices_to_indices(std::vector<std::vector<Point>> &triangles)
{
    VertexContainer container;
    int cnt = 0;
    for (auto &triangle: triangles)
    {
        std::vector<int> indices;
        for (auto &vertex: triangle)
        {
            if (container.vertexMap.count(vertex) == 0)
            {
                container.vertexMap[vertex] = cnt;
                cnt++;
            }
            indices.push_back(container.vertexMap[vertex]);
        }
        container.triangles.push_back(indices);
    }

    return container;
}

void write_to_ply(std::vector<std::vector<Point>> &triangles, const char* path)
{
    VertexContainer container = hash_vertices_to_indices(triangles);
    
    std::ofstream outputFile;
    outputFile.open(path);

    outputFile << "ply\n";
    outputFile << "format ascii 1.0\n";
    outputFile << "element vertex " << container.vertexMap.size() << "\n";
    outputFile << "property float32 x\n"; 
    outputFile << "property float32 y\n";
    outputFile << "property float32 z\n";
    outputFile << "element face " << container.triangles.size() << "\n";
    outputFile << "property list uint8 int32 vertex_indices\n";
    outputFile << "end_header\n";

    std::vector<Point> vertices (container.vertexMap.size());
    for (auto &vertex: container.vertexMap)
        vertices[vertex.second] = vertex.first;

    for (auto &vertex: vertices)
        outputFile << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    for (auto &triangle: container.triangles)
    {
        outputFile << 3 << " ";
        for (int index: triangle)
            outputFile << index << " ";
        outputFile << "\n";
    }
}

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

std::vector<std::vector<Point>> MarchingCubes(options *opts, meshInfo *mesh)
{
    /*
        Function Marching Cubes:
        Inputs:
            - pointer to TPMS
            - pointer to usr options
            - pointer to mesh struct
        Outputs:
            - None
        
        Function will get TPMS structure as input, calculate the surface mesh using
        the marching cubes method. Triangle locations are stored in the MCT
        struct.
    */

    // Get triangles

    std::vector<std::vector<Point>> triangles;

    int row, col, slice;

    // iterate over the whole image
    for(int i = 0; i < mesh->nElements; i++)
    {
        // Break down i into the integer indexes
        int slice = i / (mesh->numCellsX * mesh->numCellsX);
        int row = (i - slice * mesh->numCellsX * mesh->numCellsX) / mesh->numCellsX;
        int col = i - slice * mesh->numCellsX * mesh->numCellsX - row * mesh->numCellsX;

        int cubeIndex = 0;

        GridCell cell;

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

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  1;
            }
            else if(j == 1)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  2;
            }
            else if(j == 2)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  4;
            }
            else if(j == 3)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  8;
            }
            else if(j == 4)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  16;
            }
            else if(j == 5)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz + mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues))
                    cubeIndex |=  32;
            }
            else if(j == 6)
            {
                float x = -PI + (float)col * mesh->dx + mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
                // update bits on cubeIndex accordingly
                if(fabs(f) < fabs(opts->isoValues)) 
                    cubeIndex |=  64;
            }
            else if(j == 7)
            {
                float x = -PI + (float)col * mesh->dx;
                float y = -PI + (float)row * mesh->dy + mesh->dy;
                float z = -PI + (float)slice * mesh->dz;

                cell.vertex[j].x = x;
                cell.vertex[j].y = y;
                cell.vertex[j].z = z;

                float f = TPMS_F[opts->TPMS_Type - 1](x, y, z);

                cell.value[j] = f;
                
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
        int idx = 0;

        std::vector<Point> intercepts(12);
        
        while(edgeKey)
        {
            if(edgeKey&1)
            {
                int v1 = edgeToVertices[idx].first;
                int v2 = edgeToVertices[idx].second;

                Point p1 = cell.vertex[v1];
                Point p2 = cell.vertex[v2];

                // interpolate

                /*

                    In the future, we can use the opt algorithm here.
                
                */

                Point interceptPoint;

                /*
                
                    New Approach: Same as Immersed Boundary method

                */

                interceptPoint = temp_optMC(opts, mesh, p1, p2);

                intercepts[idx] = interceptPoint;

                /* 
                
                    Old Approach: Linear Interpolation
                
                */
                
                // float intIso = (opts->isoValues - fabs(cell.value[v1]))/fabs((cell.value[v2] - cell.value[v1]));

                // interceptPoint.x = intIso * (p2.x - p1.x) + p1.x;
                // interceptPoint.y = intIso * (p2.y - p1.y) + p1.y;
                // interceptPoint.z = intIso * (p2.z - p1.z) + p1.z;

                // intercepts[idx] = interceptPoint;
                // if(fabs(interceptPoint.x) > PI || fabs(interceptPoint.y) > PI || fabs(interceptPoint.z) > PI)
                // {
                //     printf("Intercept: %1.3e,%1.3e,%1.3e\n", interceptPoint.x, interceptPoint.y, interceptPoint.z);
                //     printf("Iso: %1.3e,v1 %1.3e, v2 %1.3e\n", opts->isoValues, cell.value[v1], cell.value[v2]);
                // }
            }
            idx++;
            edgeKey >>= 1;
        }

        // save triangles

        std::vector<std::vector<Point>> cellTriangles;
        
        for(int j = 0; (int)triangleTable[cubeIndex][j] != -1; j += 3)
        {
            std::vector<Point> tri (3);
            for(int k = 0; k < 3; k++)
            {
                tri[k] = intercepts[triangleTable[cubeIndex][j+k]];
            }
            cellTriangles.push_back(tri);
        }

        for(int j = 0; j < (int)cellTriangles.size(); j++)
        {
            triangles.push_back(cellTriangles[j]);
        }

    } // endfor    

    return triangles;
}


#endif