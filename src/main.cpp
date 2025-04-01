#include <lib/surfaceArea.hpp>
#include <lib/TPMS_definitions.hpp>
#include <lib/data_structures.hpp>
#include <math.h>


int main(int argc, char **argv)
{
    meshInfo mesh;
    saveInfo save;

    int sideLength = 100;
    int nVoxels = pow(sideLength, 3);
    char *P = (char *)malloc(sizeof(char) * nVoxels);
    printf("Here!");

    C_Y(P, 0.3, sideLength, &mesh);

    SA(P, &mesh, &save);

    printf("SSA = %1.3e\n", save.SA);

    return 0;
}