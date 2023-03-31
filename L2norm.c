#include "image_quilting.h"
#include <math.h>
#include <assert.h>


double l2norm(slice_t s1, slice_t s2)
{
    double error = 0;
    assert(s1.width == s2.width);
    assert(s1.height == s2.height);
    assert(s1.channels == s2.channels);

    for (int i = 0; i < s1.height; i++)
    {
        for (int j = 0; j < s1.channels * s1.width; j++)
        {
            double s1_data = s1.data[i * s1.jumpsize + j];
            double s2_data = s2.data[i * s2.jumpsize + j];

            error += pow(s1_data - s2_data, 2);
        }
    }
    return sqrt(error);
}