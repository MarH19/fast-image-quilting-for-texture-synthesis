#include "image_quilting.h"
#include <math.h>
#include <assert.h>

error_t l2norm(slice_t s1, slice_t s2)
{
    error_t error = 0;
    assert(s1.width == s2.width);
    assert(s1.height == s2.height);
    assert(s1.channels == s2.channels);

    for (int i = 0; i < s1.height; i++)
    {
        for (int j = 0; j < s1.channels * s1.width; j++)
        {
            pixel_t s1_data = s1.data[i * s1.jumpsize + j];
            pixel_t s2_data = s2.data[i * s2.jumpsize + j];
            error_t diff = (error_t) s1_data - (error_t) s2_data;
            error += diff * diff;
        }
    }
    return error;
}