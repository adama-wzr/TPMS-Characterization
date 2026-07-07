#ifndef _UTIL
#define _UTIL

#include <data_structures.hpp>
#include <math.h>

float VectorMagnitude(Point v)
{
    float mag;

    mag = sqrt(pow(v.x,2) + pow(v.y,2) + pow(v.z,2));

    return mag;
}

Point CrossProduct(Point v1, Point v2)
{

    Point CrossProduct;

    CrossProduct.x = v1.y * v2.z - v1.z * v2.y;
    CrossProduct.y = v1.z * v2.x - v1.x * v2.z;
    CrossProduct.z = v1.x * v2.y - v1.y * v2.x;

    return CrossProduct;
} 

#endif
