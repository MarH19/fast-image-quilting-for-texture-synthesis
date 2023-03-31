#include "image_quilting.h"
#include <math.h>

double error;

double l2norm(slice_t inp_slice, slice_t out_slice)
{
    for (int i = 0; i < inp_slice.height; i++)
    {
        for (int j = 0; j < (inp_slice.channels * inp_slice.width); j++)
        {

            double inp_data = inp_slice.data[i * inp_slice.jumpsize + j];
            double out_data = out_slice.data[i * inp_slice.jumpsize + j];

            error += pow(inp_data - out_data, 2);
        }
    }

    return sqrt(error);
}