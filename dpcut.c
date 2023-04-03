#include "image_quilting.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// slice_1 and slice_2 of same size are merged, c = 0 for vertical case, c = 1 for horizontal case
void calcErrors(slice_t slice_1, slice_t slice_2, pixel_t *errors);

// this implementation currently only works for vertical overlapping blocks
void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int c)
{
    // initialize errors array
    // does /3 work
    int width = slice_1.width;
    int height = slice_1.height;
    pixel_t *errors = malloc(sizeof(pixel_t) * width * height);

    calcErrors(slice_1, slice_2, errors);

    if (c == 1)
    {
        // transpose errors

        // probably best/easier to just code two cases
    }

    // initialize dp array
    pixel_t *dp = malloc(sizeof(pixel_t) * width * height);

    // initialize the first row of the dp array with the first row of errors
    for (int i = 0; i < width; i++)
    {
        dp[i] = errors[i];
    }

    // watch out this changes when transpose
    pixel_t max_val = pow(255, 2) * 3 * height + 1;

    // fill dp table
    for (int i = 1; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            // create array of size 3 for the 3 possible paths
            // since we take the min i initialize with the maximum possible value
            pixel_t paths[3] = {max_val, max_val, max_val};
            paths[1] = dp[(i - 1) * width + j];

            if (j != 0)
            {
                paths[0] = dp[(i - 1) * width + j - 1];
            }

            if (j != width - 1)
            {
                paths[2] = dp[(i - 1) * width + j + 1];
            }

            // calculate min element of paths
            // precendence: middle, left, right
            pixel_t min_path = paths[0] < paths[1] ? (paths[0] < paths[2] ? paths[0] : paths[2]) : (paths[2] < paths[1] ? paths[2] : paths[1]);
            dp[i * width + j] = errors[i * width + j] + min_path;
        }
    }
    //return errors;
    //return dp;
    // do backtracking while filling out
    // aliasing shouldn't be a problem but i need to write update the out slice/image

    // find the min element in the last row of the dp table
    int start;
    pixel_t min = max_val;

    for (int i = 0; i < width; i++)
    {
        pixel_t temp_min = dp[(height - 1) * width + i] < max_val ? dp[(height - 1) * width + i] : max_val;
        if (temp_min < min)
        {
            start = i;
            min = temp_min;
        }
    }
    // set last row of out image
    for (int i = 0; i < slice_1.width * 3; i++)
    {
        out.data[(out.height - 1) * out.jumpsize + i] = i/3 < start ? slice_1.data[(height - 1) * slice_1.jumpsize + i] : slice_2.data[(slice_2.height - 1) * slice_2.jumpsize + i];
    }


    // for all rows above do
    // find new start
    // set current row of out image
    for (int i = height - 2; i >= 0; i--)
    {
        // is the j loop needed?
        // find minimum in 3 possible paths
        pixel_t paths[3] = {max_val, max_val, max_val};
        paths[1] = dp[i * width + start];

        if (start != 0)
        {
            paths[0] = dp[i * width + start - 1];
        }

        if (start != width - 1)
        {
            paths[2] = dp[i * width + start + 1];
        }

        pixel_t min_path = paths[0] < paths[1] ? (paths[0] < paths[2] ? paths[0] : paths[2]) : (paths[2] < paths[1] ? paths[2] : paths[1]);

        if (start < width - 1)
        {
            if (dp[width * i + start + 1] == min_path)
            {
                start++;
            }
        }
        if (start > 0)
        {
            if (dp[width * i + start + -1] == min_path)
            {
                start--;
            }
        }
        printf("%d\n", start);
        // fill the row of the output
        for (int j = 0; j < out.width * 3; j++)
        {
            out.data[i * out.jumpsize + j] = j/3 < start ? slice_1.data[i * slice_1.jumpsize + j] : slice_2.data[i * slice_2.jumpsize + j];
        }
    }
    // free errors
    free(errors);
    // free dp
    free(dp);

}

// errors(i,j) = sum of squared differences of the 3 rgb values
void calcErrors(slice_t slice_1, slice_t slice_2, pixel_t *errors)
{

    for (int i = 0; i < slice_1.height; i++)
    {
        for (int j = 0; j < slice_1.width * 3; j += 3)
        {
            // ssd for one pixel
            pixel_t error = (pow(slice_1.data[i * slice_1.jumpsize + j] - slice_2.data[i * slice_2.jumpsize + j], 2)) + (pow(slice_1.data[i * slice_1.jumpsize + j + 1] - slice_2.data[i * slice_2.jumpsize + j + 1], 2)) + (pow(slice_1.data[i * slice_1.jumpsize + j + 2] - slice_2.data[i * slice_2.jumpsize + j + 2], 2));
            errors[i * slice_1.width + j/3] = error;
        }
    }
}
