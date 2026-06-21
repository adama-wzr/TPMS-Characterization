
#ifndef _MC
#define _MC

#include <data_structures.hpp>
#include <constants.hpp>
#include <TPMS_definitions.hpp>

float zeroFunction(options *opts, Point p)
{
    /*
        Function zeroFunction:
        Inputs:
            - pointer to opts
            - point (x,y,z)
        Output:
            - returns r = C - |S_a(x,y,z)|, where 
            C is the isovalue and S_a is the TPMS
            function evaluated at px, py, and pz.
    */
    float result = opts->isoValues - fabs(TPMS_F[opts->TPMS_Type - 1](p.x, p.y, p.z));

    return result;
}

Point bisectionMethod(options *opts, meshInfo *mesh, Point p1, Point p2, float v1, float v2)
{
    /*
        Function bisectionMethod:
        Inputs:
            - pointer to opts
            - pointer to mesh
            - point p1
            - point p2
            - isovalue at p1
            - isovalue at p2
        Outputs:
            - TPMS intercept point
        
        Function implements the bisection method to find the exact TPMS
        intercept location based on a marching cubes algorithm.
    */

    // setup

    float tol = 1e-4;
    int maxIter = 1e4;
    int iteration = 0;
    float eps;

    Point interceptPoint;

    Point a, b, c;

    // Step 1: check that the TPMS crosses in this interval

    float p1_value = zeroFunction(opts, p1);
    float p2_value = zeroFunction(opts, p2);

    // check if either p1 or p2 are below tol

    if(p1_value <= tol)
        return p1;
    else if(p2_value <= tol)
        return p2;

    // if p1.v or p2.v is within eps, return

    if(p1_value < 0 && p2_value > 0)
    {
        a = p1;
        b = p2;
    }
    else if(p1_value > 0 && p2_value < 0)
    {
        a = p2;
        b = p1;
    }
    else
    {
        printf("Error! Intercept is out of range!\n");
        // Simple interpolation
        float intIso = (opts->isoValues - fabs(v1))/fabs((v2 - v1));

        interceptPoint.x = intIso * (p2.x - p1.x) + p1.x;
        interceptPoint.y = intIso * (p2.y - p1.y) + p1.y;
        interceptPoint.z = intIso * (p2.z - p1.z) + p1.z;
        return interceptPoint;
    }

    // find c

    c.x = (p1.x + p2.x)/2;
    c.y = (p1.y + p2.y)/2;
    c.z = (p1.z + p2.z)/2;

    eps = fabs(zeroFunction(opts, c));

    if(eps <= tol)
        return c;
    
    // otherwise, iterate
    while (eps > tol && iteration < maxIter)
    {
        // assign new bounds
        if(eps > 0)
            b = c;
        else if(eps < 0)
            a = c;
        
        // find new c

        c.x = (a.x + b.x)/2;
        c.y = (a.y + b.y)/2;
        c.z = (a.z + b.z)/2;

        eps = fabs(zeroFunction(opts, c));

        iteration++;
    }

    // assign intercept point
    interceptPoint = c;

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

                interceptPoint = bisectionMethod(opts, mesh, p1, p2, cell.value[v1], cell.value[v2]);

                intercepts[idx] = interceptPoint;
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
