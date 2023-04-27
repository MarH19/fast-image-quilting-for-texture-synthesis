#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "image_quilting.h"
/* What do we need to do?
 *  1. load the image
 *  2. calculate all the parameters
 *  3. malloc the output image
 *  4. call function for random patch
 *  5. copy patches to something
 *  6. find next random patch
 *  7.
 */

/*
Calculates the error between the slice of an input block (i.e., for each block of the input image) and the slice of
a fixed output block
*/
void calc_errors(image_t in, int blocksize, slice_t in_slice, slice_t out_slice, pixel_t *errors, int add)
{
    int error_jumpsize = in.width - blocksize + 1;

    if (add)
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {

                in_slice.data = in.data + in_slice.jumpsize * i + j * in_slice.channels;
                errors[i * error_jumpsize + j] += l2norm(in_slice, out_slice); // if add > 0, the operation is plus (+)
            }
        }
    }
    else
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {

                in_slice.data = in.data + in_slice.jumpsize * i + j * in_slice.channels;
                pixel_t l2 = l2norm(in_slice, out_slice);
                if (j == 1)
                {
                    // printf("l2: %lf\n", l2);
                }
                errors[i * error_jumpsize + j] -= l2norm(in_slice, out_slice); // if add == 0, the operation is minus (-)
            }
        }
    }
}

/*
The function find finds the coordinates of a candidate block that is within a certain range (specified by tolerance)
from the best fitting block
 */
static int counter = 0;
coord find(pixel_t *errors, int height, int width, pixel_t tolerance)
{

    pixel_t min_error = INFINITY;

    // search for the minimum error in the errors array and store it in a variable.
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (errors[i * width + j] < min_error)
            {

                min_error = errors[i * width + j];
            }
        }
    }
    //printf("err: %lf\n", min_error);

    // Count how many canditates exist in order to know the size of the array of candidates

    pixel_t tol_range = (1.0 + tolerance) * min_error;
    int nr_candidates = 0;
    // printf("toll: %lf, height: %d, width: %d\n", tol_range, height, width);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (errors[i * width + j] <= tol_range)
            {
                nr_candidates++;
            }
        }
    }

    // Create the array with candidates and populate it with
    coord candidates[nr_candidates];
    int idx = 0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {

            if (errors[i * width + j] <= tol_range)
            {
                coord candid = {i, j};
                candidates[idx] = candid;
                idx++;
            }
        }
    }
    // Choose randomly a candidate
    // printf("cand: %d\n", nr_candidates);
    int random_idx = rand() % nr_candidates;
    // printf("%d \n", counter++);
    coord random_candidate = candidates[random_idx];

    // return coordinates of the random candidate
    return random_candidate;
}

// Command to compile the code:   gcc dpcut.c imageio.c L2norm.c  -o imageio -lm
image_t image_quilting(image_t in, int blocksize, int num_blocks, int overlap, pixel_t tolerance)
{

    int out_size = num_blocks * blocksize - (num_blocks - 1) * overlap;
    int n = out_size * out_size * 3;
    double *out_image = (double *)calloc(n, sizeof(double));
    assert(out_image);
    image_t out = {out_image, out_size, out_size, 3};
    int errorlen = (in.height - blocksize + 1) * (in.width - blocksize + 1) * sizeof(pixel_t);
    pixel_t *errors = (pixel_t *)malloc(errorlen);
    assert(errors);

    for (int row = 0; row < num_blocks; row++)
    {
        for (int col = 0; col < num_blocks; col++)
        {
            printf("row: %d, col: %d\n", row, col);
            int si = row * (blocksize - overlap);
            int sj = col * (blocksize - overlap);
            memset(errors, 0, errorlen);

            // Very first case, so pick one at random
            if (row == 0 && col == 0)
            {
                int row_idx = 0; //rand() % (in.height - blocksize + 1);
                int col_idx = 0; //rand() % (in.width - blocksize + 1);

                slice_t out_block = slice_image(out, si, sj, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, row_idx, col_idx, row_idx + blocksize, col_idx + blocksize);

                slice_cpy(in_block, out_block);
                continue;
            }
            if (row != 0) // safe to check top
            {
                slice_t out_slice = slice_image(out, si, sj, si + overlap, sj + blocksize);
                slice_t in_slice = slice_image(in, 0, 0, overlap, blocksize);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }

            if (col != 0) // safe to check left
            {
                slice_t out_slice = slice_image(out, si, sj, si + blocksize, sj + overlap);
                slice_t in_slice = slice_image(in, 0, 0, blocksize, overlap);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }
            if (col != 0 && row != 0) // overlap so remove common
            {
                slice_t out_slice = slice_image(out, si, sj, si + overlap, sj + overlap);
                slice_t in_slice = slice_image(in, 0, 0, overlap, overlap);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 0);
            }
            // Search for a random candidate block among the best matching blocks (determined by the tolerance)
            coord random_candidate = find(errors, in.height - blocksize + 1, in.width - blocksize + 1, tolerance);
            int rand_row = random_candidate.row;
            int rand_col = random_candidate.col;
            

            if (row == 0)
            {
                slice_t out_block = slice_image(out, si, sj+overlap, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, rand_row, rand_col+overlap, rand_row + blocksize, rand_col + blocksize);
                slice_cpy(in_block, out_block);
            }
            else if (col == 0)
            {
                slice_t out_block = slice_image(out, si+overlap, sj, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, rand_row+overlap, rand_col, rand_row + blocksize, rand_col + blocksize);
                slice_cpy(in_block, out_block);
            }
            else
            {
                slice_t out_block = slice_image(out, si+overlap, sj+overlap, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, rand_row+overlap, rand_col+overlap, rand_row + blocksize, rand_col + blocksize);
                slice_cpy(in_block, out_block);
            }
            // slice_t out_block = slice_image(out, si, sj, si + blocksize, sj + blocksize);
            // slice_t in_block = slice_image(in, rand_row, rand_col, rand_row + blocksize, rand_col + blocksize);
            // slice_cpy(in_block, out_block);

            if (row != 0)
            {
                slice_t out_slice = slice_image(out, si, sj, si + overlap, sj + blocksize);
                slice_t in_slice = slice_image(in, rand_row, rand_col, rand_row + overlap, rand_col + blocksize);
                dpcut(out_slice, in_slice, out_slice, 1);
            }
            if (col != 0)
            {
                slice_t out_slice = slice_image(out, si, sj, si + blocksize, sj + overlap);
                slice_t in_slice = slice_image(in, rand_row, rand_col, rand_row + blocksize, rand_col + overlap);
                dpcut(out_slice, in_slice, out_slice, 0);
            }
            // if (row > 1 && col >= 1){
            //     exit(-1);
            // }
        }
    }

    free(errors);

    return out;
}
// gcc dpcut.c imageio.c L2norm.c  -o imageio -lm
