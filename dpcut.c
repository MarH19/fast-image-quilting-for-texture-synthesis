#include "image_quilting.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void calcErrors(slice_t slice_1, slice_t slice_2, pixel_t *errors);
pixel_t *transpose(pixel_t *mat, int width, int height);

// slice_1 and slice_2 of same size are merged into out, c = 0 for vertical case, c = 1 for horizontal case
// flop count: flops(transpose) + flops(calcerrors) + (s_width * s_height) * (3 * min + add) + s_width * (min + LT) + (s_height - 1) * (3 * min + 2 * EQ)
// this is approximately: 0 + (s_height * s_width) * 3 * 2 + (s_width * s_height) * (3 + 1) + s_width * 2 + (s_height-1) * (3 + 2) 
// = (s_height * s_width) * 10 + s_width * 2 + (s_height -1) * 5 flops  
void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int c)
{

    int width;
    int height;

    if (c == 1)
    {
        width = slice_1.height;
        height = slice_1.width;
    }
    else
    {
        width = slice_1.width;
        height = slice_1.height;
    }

    // initialize errors array
    pixel_t *errors = malloc(sizeof(pixel_t) * width * height);

    calcErrors(slice_1, slice_2, errors);

    if (c == 1)
    {
        // transpose errors (now in row format)
        errors = transpose(errors, slice_1.width, slice_1.height);
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

    //  do backtracking while filling out slice

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
    if (c == 1){
        for (int i = 0; i < out.height; i++)
        {
            for (int k = 0; k < out.channels; k++)
            {
                out.data[i * out.jumpsize + (out.width - 1) * 3 + k] = i < start ? slice_1.data[i * slice_1.jumpsize + (slice_1.width - 1) * 3 + k] : slice_2.data[(i * slice_2.jumpsize + slice_2.width - 1) * 3 + k];
            }
        }
    } else {
        for (int i = 0; i < slice_1.width * 3; i++)
        {
            out.data[(out.height - 1) * out.jumpsize + i] = i / 3 < start ? slice_1.data[(height - 1) * slice_1.jumpsize + i] : slice_2.data[(slice_2.height - 1) * slice_2.jumpsize + i];
        }
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
        //printf("%d\n", start);
        // fill the row of the output
        if (c == 1){
            for (int j = 0; j < out.height; j++)
            {
                for (int k = 0; k < out.channels; k++)
                {
                    out.data[j * out.jumpsize + i * 3 + k] = j < start ? slice_1.data[j * slice_1.jumpsize + i * 3 + k] : slice_2.data[j * slice_2.jumpsize + i * 3 + k];
                }
            }
        } else {
            for (int j = 0; j < out.width * 3; j++)
            {
                out.data[i * out.jumpsize + j] = j / 3 < start ? slice_1.data[i * slice_1.jumpsize + j] : slice_2.data[i * slice_2.jumpsize + j];
            }
        }
    }
    // free errors
    free(errors);
    // free dp
    free(dp);
}

// errors(i,j) = sum of squared differences of the 3 rgb values
// flop count: (s_height * s_width) * 3 * (pow + add)
void calcErrors(slice_t slice_1, slice_t slice_2, pixel_t *errors)
{

    for (int i = 0; i < slice_1.height; i++)
    {
        for (int j = 0; j < slice_1.width * 3; j += 3)
        {
            // ssd for one pixel
            pixel_t error = (pow(slice_1.data[i * slice_1.jumpsize + j] - slice_2.data[i * slice_2.jumpsize + j], 2)) + (pow(slice_1.data[i * slice_1.jumpsize + j + 1] - slice_2.data[i * slice_2.jumpsize + j + 1], 2)) + (pow(slice_1.data[i * slice_1.jumpsize + j + 2] - slice_2.data[i * slice_2.jumpsize + j + 2], 2));
            errors[i * slice_1.width + j / 3] = error;
        }
    }
}

// transpose a matrix
// flop count: 0
pixel_t *transpose(pixel_t *mat, int width, int height)
{
    pixel_t *mat_t = malloc(sizeof(pixel_t) * width * height);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            mat_t[j * height + i] = mat[i * width + j];
        }
    }
    free(mat);
    mat = NULL;
    return mat_t;
}
