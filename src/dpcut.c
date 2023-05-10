#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "image_quilting.h"

static void calc_errors(slice_t slice_1, slice_t slice_2, pixel_t *errors);
static pixel_t *transpose(pixel_t *mat, int width, int height);

// slice_1 and slice_2 of same size are merged into out, left2right = 0 for vertical case, left2right = 1 for left2right case
// flop count: flops(transpose) + flops(calc_errors) + (s_width * s_height) * (3 * min + add) + s_width * (min + LT) + (s_height - 1) * (3 * min + 2 * EQ)
// this is approximately: 0 + (s_height * s_width) * 3 * 2 + (s_width * s_height) * (3 + 1) + s_width * 2 + (s_height-1) * (3 + 2)
// = (s_height * s_width) * 10 + s_width * 2 + (s_height -1) * 5 flops
void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int left2right)
{
    int width;
    int height;

    if (left2right)
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
    pixel_t *dp = malloc(sizeof(pixel_t) * width * height);

    calc_errors(slice_1, slice_2, dp);

    // transpose errors (now in row format)
    if (left2right)
    {
        pixel_t *temp = transpose(dp, slice_1.width, slice_1.height);
        free(dp);
        dp = NULL;
        dp = temp;
    }

    // fill dp table
    for (int i = 1; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            // per default check middle
            pixel_t prev = dp[(i - 1) * width + j];
            // check left
            if (j > 0)
                prev = MIN(prev, dp[(i - 1) * width + j - 1]);
            // check right
            if (j < width - 1)
                prev = MIN(prev, dp[(i - 1) * width + j + 1]);
            dp[i * width + j] += prev;
        }
    }

    //  do backtracking while filling out slice
    // find the min element in the last row of the dp table
    int previdx = 0;
    pixel_t prevval = INT_MAX;
    for (int i = 0; i < width; i++)
    {
        pixel_t temp_min = dp[(height - 1) * width + i];
        if (temp_min < prevval)
        {
            previdx = i;
            prevval = temp_min;
        }
    }

    // for all rows find nextidx, nextval
    // set current row of out image
    for (int i = height - 1; i >= 0; i--)
    {
        int nextidx = previdx;
        pixel_t nextval = dp[i * width + nextidx];
        if (previdx > 0 && dp[i * width + previdx - 1] < nextval)
        {
            nextval = dp[i * width + previdx - 1];
            nextidx = previdx - 1;
        }
        if (previdx < width - 1 && dp[i * width + previdx + 1] < nextval)
            nextidx = previdx + 1;
        //  fill the row of the output
        if (left2right)
        {
            for (int j = 0; j < out.height; j++)
            {
                for (int k = 0; k < out.channels; k++)
                {
                    int idxo = j * out.jumpsize + i * out.channels + k;
                    int idx1 = j * slice_1.jumpsize + i * out.channels + k;
                    int idx2 = j * slice_2.jumpsize + i * out.channels + k;
                    out.data[idxo] = j < nextidx ? slice_1.data[idx1] : slice_2.data[idx2];
                }
            }
        }
        else
        {
            for (int j = 0; j < out.width * out.channels; j++)
            {
                int idxo = i * out.jumpsize + j;
                int idx1 = i * slice_1.jumpsize + j;
                int idx2 = i * slice_2.jumpsize + j;
                out.data[idxo] = (j / out.channels) < nextidx ? slice_1.data[idx1] : slice_2.data[idx2];
            }
        }
        previdx = nextidx;
    }
    free(dp);
}

// errors(i,j) = sum of squared differences of the 3 rgb values
// flop count: (s_height * s_width) * 3 * (pow + add)
static void calc_errors(slice_t s1, slice_t s2, pixel_t *errors)
{
    for (int i = 0; i < s1.height; i++)
    {
        for (int j = 0; j < s1.width; j++)
        {
            // ssd for one pixel
            pixel_t error = 0.0;
            for (int c = 0; c < s1.channels; c++)
            {
                int idx1 = i * s1.jumpsize + j * s1.channels + c;
                int idx2 = i * s2.jumpsize + j * s2.channels + c;
                error += pow(s1.data[idx1] - s2.data[idx2], 2);
            }
            errors[i * s1.width + j] = error;
        }
    }
}

// transpose a matrix
// flop count: 0
/* transpose: transposes matrix into new matrix */
static pixel_t *transpose(pixel_t *mat, int width, int height)
{
    pixel_t *mat_t = malloc(sizeof(pixel_t) * width * height);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            mat_t[j * height + i] = mat[i * width + j];
        }
    }
    return mat_t;
}
