#include "image_quilting.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// slice_1 and slice_2 of same size are merged, c = 0 for vertical case, c = 1 for horizontal case
void calcErrors(slice_t slice_1, slice_t slice_2, pixel_t *errors);

pixel_t* dpcut(slice_t slice_1, slice_t slice_2, int c)
{
    // initialize errors array
    pixel_t *errors = malloc(sizeof(pixel_t) * slice_1.width * slice_1.height);
    calcErrors(slice_1, slice_2, errors);
    
    if (c == 1)
    {
        // transpose errors
    }

    // initialize dp array

    // fill dp table

    // do backtracking while filling slice_1 -- aliasing a problem?

    if (c == 1)
    {
        // traspose back?? how can i do this while i write?
    }
    return errors;
}

// errors(i,j) = sum of squared differences of the 3 rgb values
void calcErrors(slice_t slice_1, slice_t slice_2, pixel_t *errors)
{
    
    for (int i = 0; i < slice_1.height; i++)
    {
        for (int j = 0; j < slice_1.width; j++)
        {
            // ssd for one pixel
            pixel_t error = (pow(slice_1.data[i*slice_1.jumpsize + j] - slice_2.data[i*slice_2.jumpsize + j], 2)) 
                            + (pow(slice_1.data[i*slice_1.jumpsize + j + 1] - slice_2.data[i*slice_2.jumpsize + j + 1], 2)) 
                            + (pow(slice_1.data[i*slice_1.jumpsize + j + 2] - slice_2.data[i*slice_2.jumpsize + j + 2], 2));
            errors[i*slice_1.width + j] = error;
        }
    }
}
