/*
    SubDomain Partitioning w/ Flood-Fill approach.


    Andre Adam
    07/12/2025
*/

#ifndef _SUBFF
#define _SUBFF

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <data_structures.hpp>

void subDomainFC(meshInfo *mesh, char *subDomain)
{
    /*
        Function subDomainFC:
        Inputs:
            - pointer to mesh struct
            - pointer to subDomain array
        Outputs:
            - none
        
        For each subdomain, the function will determine if
        the sub-domain is fully connected from boundary to
        boundary, assuming periodic boundary conditions.

        This means the domain allows for flow between boundaries,
        and then properties of interest can be assigned.
    */

    // useful variables

    int nCols = mesh->numCellsX;
    int nRows = mesh->numCellsY;
    int nSlices = mesh->numCellsZ;

    bool *checked = (bool *)malloc(sizeof(bool) * mesh->nElements);

    // search each sub-domain fully
    for(int nSub = 1; nSub <= mesh->nChannels; nSub++)
    {
        // set FC to 0 by default

        mesh->sdInfo[nSub-1].FC = 0;

        // initialize checked to 0
        memset(checked, 0, mesh->nElements);
        
        // open list

        std::set<coord> cList;

        // local VF

        long int count = 0;

        // start from the left side

        for(int row = 0; row < nCols; row++)
        {
            for(int slice = 0; slice < nSlices; slice ++)
            {
                if(subDomain[slice *nRows * nCols + row * nCols] == nSub)
                {
                    cList.insert(std::tuple(0, row, slice));
                    count++;
                }
                checked[slice *nRows * nCols + row * nCols] = 1;
            }
        }

        // if there were none at the boundary, FC = false

        if(count == 0)
        {
            mesh->sdInfo[nSub-1].FC = 0;
        }

        // if there were thing at the boundary, scan for FC

        while (!cList.empty())
        {
            // pop first item on the list
            coord pop = *cList.begin();

            // remove from open list
            cList.erase(cList.begin());

            // get coordinates from the list
            int col = std::get<0>(pop);
            int row = std::get<1>(pop);
            int slice = std::get<2>(pop);

            if(col == nCols - 1)
            {
                mesh->sdInfo[nSub - 1].FC = 1;
            }

            /*
                We need to check North, South, East, West, Back, and Front for more fluid:

                North = col + 0, row - 1, slice + 0
                South = col + 0, row + 1, slice + 0
                East  = col + 1, row + 0, slice + 0
                West  = col - 1, row + 0, slice + 0
                Front = col + 0, row + 0, slice - 1
                Back  = col + 0, row + 0, slice + 1

                Note that diagonals are not considered a connection.
            */
           
            int tempRow, tempCol, tempSlice;
            long int tempIndex;

            // North

            tempCol = col;
            tempSlice = slice;

            if (row != 0)
            {
                tempRow = row - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            } else
            {
                // Periodic
                tempRow = mesh->numCellsY - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            }

            if (subDomain[tempIndex] == nSub && checked[tempIndex] == 0)
            {
                count++;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }
            
            checked[tempIndex] = 1;

            // South

            tempCol = col;
            tempSlice = slice;

            if (row != mesh->numCellsY - 1)
            {
                tempRow = row + 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            } else
            {
                // Periodic
                tempRow = 0;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            }

            if (subDomain[tempIndex] == nSub && checked[tempIndex] == 0)
            {
                count++;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }

            checked[tempIndex] = 1;

            // Front

            tempCol = col;
            tempRow = row;

            if (slice != 0)
            {
                tempSlice = slice - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            } else
            {
                // Periodic
                tempSlice = mesh->numCellsZ - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            }

            if (subDomain[tempIndex] == nSub && checked[tempIndex] == 0)
            {
                count++;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }

            checked[tempIndex] = 1;

            // Back

            if (slice != mesh->numCellsZ - 1)
            {
                tempSlice = slice + 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            }
            else
            {
                // Periodic
                tempSlice = 0;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
            }

            if (subDomain[tempIndex] == nSub && checked[tempIndex] == 0)
            {
                count++;
                cList.insert(std::tuple(tempCol, tempRow, tempSlice));
            }

            checked[tempIndex] = 1;
            

            // West

            tempRow = row;
            tempSlice = slice;

            if (col != 0)
            {
                tempCol = col - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomain[tempIndex] == nSub && checked[tempIndex] == 0)
                {
                    count++;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
                checked[tempIndex] = 1;
            }

            // East

            tempRow = row;
            tempSlice = slice;

            if (col != mesh->numCellsX - 1)
            {
                tempCol = col + 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomain[tempIndex] == nSub && checked[tempIndex] == 0)
                {
                    count++;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
                checked[tempIndex] = 1;
            }
        } // end while

        // calculate VF for this method

        mesh->sdInfo[nSub-1].VF = (float)count / mesh->nElements;
        printf("Sub-Domain %d, VF = %1.3e\n", nSub, mesh->sdInfo[nSub-1].VF);
    }

    return;
}

void subDomainFF(meshInfo *mesh, char *P, char *subDomains)
{
    /*
        Funciton subDomainFF:
        Inputs:
            - pointer to mesh struc
            - pointer to binary structure char
            - pointer to subDomains array
        Outputs:
            - none.
        
        Function will get the number of independent channels. While
        this should be two for all TPMS structures, there might be some
        incongruent types that, as we increase the porosity, have pinch
        points that further create more disconnected channels.

        This function assumes periodic boundary conditions (otherwise
        TPMS don't make sense).

        Function starts scanning at the top, left, front corner, 
        with index (0, 0, 0). It will progressively scan all of the
        pixels on the entire domain that are connected. Then it will
        proceed to look for the next unlabelled pixel with the number
        of channels incremented.
    */

    // make sure all entries in subDomain are -1 if fluid, 0 if solid.

    for(int i = 0;  i < mesh->nElements; i++)
    {
        if(P[i] == 1)
        {
            subDomains[i] = 0;
        }
        else
        {
            subDomains[i] = -1;
        }
    }

    // setup search

    int nCols = mesh->numCellsX;
    int nRows = mesh->numCellsY;

    bool scan = true;

    int lastIdxChecked = 0;

    int nChannels = 0;

    int row, col, slice;

    std::set<coord> cList;

    while (scan)
    {
        // find any fluids that haven't been assigned yet
        for(int i = lastIdxChecked; i < mesh->nElements; i++)
        {
            if (subDomains[i] == -1)
            {
                lastIdxChecked = i;

                // open lists and assign starting point
                
                slice = i / (nRows * nCols);
                row = (i - slice * nRows * nCols)/nCols;
                col = i - slice * nRows * nCols - row * nCols;

                cList.insert(std::tuple(col, row, slice));

                // increment nChannels and assign
                nChannels++;
                subDomains[i] = nChannels;

                break;
            }
        }

        if(cList.empty())
            scan = false;

        while(!cList.empty())
        {
            // pop first item on the list
            coord pop = *cList.begin();

            // remove from open list
            cList.erase(cList.begin());

            // get coordinates from the list
            col = std::get<0>(pop);
            row = std::get<1>(pop);
            slice = std::get<2>(pop);

            /*
                We need to check North, South, East, West, Back, and Front for more fluid:

                North = col + 0, row - 1, slice + 0
                South = col + 0, row + 1, slice + 0
                East  = col + 1, row + 0, slice + 0
                West  = col - 1, row + 0, slice + 0
                Front = col + 0, row + 0, slice - 1
                Back  = col + 0, row + 0, slice + 1

                Note that diagonals are not considered a connection.
                Periodic N/S, F/B
            */

            int tempRow, tempCol, tempSlice;
            long int tempIndex;

            // North

            tempCol = col;
            tempSlice = slice;

            if (row != 0)
            {
                tempRow = row - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            } else
            {
                // Periodic
                tempRow = mesh->numCellsY - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }

            // South

            if (row != mesh->numCellsY - 1)
            {
                tempRow = row + 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            } else
            {
                // Periodic
                tempRow = 0;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }

            // Front

            tempCol = col;
            tempRow = row;

            if (slice != 0)
            {
                tempSlice = slice - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            } else
            {
                // Periodic
                tempSlice = mesh->numCellsZ - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }

            // Back

            if (slice != mesh->numCellsZ - 1)
            {
                tempSlice = slice + 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }
            else
            {
                // Periodic
                tempSlice = 0;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }

            // West

            tempRow = row;
            tempSlice = slice;

            if (col != 0)
            {
                tempCol = col - 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }

            // East

            if (col != mesh->numCellsX - 1)
            {
                tempCol = col + 1;
                tempIndex = tempSlice * mesh->numCellsX * mesh->numCellsY + tempRow * mesh->numCellsX + tempCol;
                if (subDomains[tempIndex] == -1)
                {
                    subDomains[tempIndex] = nChannels;
                    cList.insert(std::tuple(tempCol, tempRow, tempSlice));
                }
            }
            // repeat until cList is empty
        } // end while
    }// end while
    
    mesh->nChannels = nChannels;
    return;
}

#endif