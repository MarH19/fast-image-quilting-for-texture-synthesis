#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_quilting.h"

/*
Calculates the error between the slice of an input block (i.e., for each block of the input image) and the slice of
a fixed output block
*/
void calc_errors_(image_t in_, image_t out_, int row_out, int col_out, int blocksize, int overlap, pixel_t *errors)
{
    int error_height = in_.height - blocksize + 1;
    int error_width = in_.width - blocksize + 1;
    pixel_t *in = in_.data;
    pixel_t *out = out_.data;
    /* calculate jumpsize aka stride */
    int inj = in_.channels * in_.width;
    int outj = out_.channels * out_.width;
    pixel_t error1, error2, error3;
    // one diff per subblock and color
    pixel_t diff11, diff12, diff13;
    pixel_t diff21, diff22, diff23;
    pixel_t diff31, diff32, diff33;

    /* iterate over all possible blocks */

    if (col_out == 0)
    {
        for (int i = 0; i < error_height; i++)
        {
            for (int j = 0; j < error_width; j++)
            {
                error1 = error2 = error3 = 0;

                /* we have 3 different subblocks
                 *  1 2 2
                 *  3 . .
                 *  3 . .
                 */
                /* I use a variant where I can calculate all three with two loops */
                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < (blocksize - overlap); m++)
                    {
                        if (m < overlap) // subblock 1
                        {
                            diff11 = in[(i + k) * inj + (j + m) * 3 + 0] - out[(row_out + k) * outj + (col_out + m) * 3 + 0];
                            diff12 = in[(i + k) * inj + (j + m) * 3 + 1] - out[(row_out + k) * outj + (col_out + m) * 3 + 1];
                            diff13 = in[(i + k) * inj + (j + m) * 3 + 2] - out[(row_out + k) * outj + (col_out + m) * 3 + 2];
                            error1 = error1 + diff11 * diff11 + diff12 * diff12 + diff13 * diff13;
                        }
                        diff31 = in[(i + m + overlap) * inj + (j + k) * 3 + 0] - out[(row_out + m + overlap) * outj + (col_out + k) * 3 + 0];
                        diff32 = in[(i + m + overlap) * inj + (j + k) * 3 + 1] - out[(row_out + m + overlap) * outj + (col_out + k) * 3 + 1];
                        diff33 = in[(i + m + overlap) * inj + (j + k) * 3 + 2] - out[(row_out + m + overlap) * outj + (col_out + k) * 3 + 2];
                        error3 = error3 + diff31 * diff31 + diff32 * diff32 + diff33 * diff33;
                    }
                }
                errors[i * error_width + j] = error1 + error2 + error3;
            }
        }
    }
    else if (row_out == 0)
    {
        for (int i = 0; i < error_height; i++)
        {
            for (int j = 0; j < error_width; j++)
            {
                error1 = error2 = error3 = 0;

                /* we have 3 different subblocks
                 *  1 2 2
                 *  3 . .
                 *  3 . .
                 */
                /* I use a variant where I can calculate all three with two loops */
                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < (blocksize - overlap); m++)
                    {
                        if (m < overlap) // subblock 1
                        {
                            diff11 = in[(i + k) * inj + (j + m) * 3 + 0] - out[(row_out + k) * outj + (col_out + m) * 3 + 0];
                            diff12 = in[(i + k) * inj + (j + m) * 3 + 1] - out[(row_out + k) * outj + (col_out + m) * 3 + 1];
                            diff13 = in[(i + k) * inj + (j + m) * 3 + 2] - out[(row_out + k) * outj + (col_out + m) * 3 + 2];
                            error1 = error1 + diff11 * diff11 + diff12 * diff12 + diff13 * diff13;
                        }
                        diff21 = in[(i + k) * inj + (j + m + overlap) * 3 + 0] - out[(row_out + k) * outj + (col_out + m + overlap) * 3 + 0];
                        diff22 = in[(i + k) * inj + (j + m + overlap) * 3 + 1] - out[(row_out + k) * outj + (col_out + m + overlap) * 3 + 1];
                        diff23 = in[(i + k) * inj + (j + m + overlap) * 3 + 2] - out[(row_out + k) * outj + (col_out + m + overlap) * 3 + 2];
                        error2 = error2 + diff21 * diff21 + diff22 * diff22 + diff23 * diff23;
                    }
                }
                errors[i * error_width + j] = error1 + error2 + error3;
            }
        }
    }
    else
    {
        for (int i = 0; i < error_height; i++)
        {
            for (int j = 0; j < error_width; j++)
            {
                error1 = error2 = error3 = 0;

                /* we have 3 different subblocks
                 *  1 2 2
                 *  3 . .
                 *  3 . .
                 */
                /* I use a variant where I can calculate all three with two loops */
                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < (blocksize - overlap); m++)
                    {
                        if (m < overlap) // subblock 1
                        {
                            diff11 = in[(i + k) * inj + (j + m) * 3 + 0] - out[(row_out + k) * outj + (col_out + m) * 3 + 0];
                            diff12 = in[(i + k) * inj + (j + m) * 3 + 1] - out[(row_out + k) * outj + (col_out + m) * 3 + 1];
                            diff13 = in[(i + k) * inj + (j + m) * 3 + 2] - out[(row_out + k) * outj + (col_out + m) * 3 + 2];
                            error1 = error1 + diff11 * diff11 + diff12 * diff12 + diff13 * diff13;
                        }
                        diff21 = in[(i + k) * inj + (j + m + overlap) * 3 + 0] - out[(row_out + k) * outj + (col_out + m + overlap) * 3 + 0];
                        diff22 = in[(i + k) * inj + (j + m + overlap) * 3 + 1] - out[(row_out + k) * outj + (col_out + m + overlap) * 3 + 1];
                        diff23 = in[(i + k) * inj + (j + m + overlap) * 3 + 2] - out[(row_out + k) * outj + (col_out + m + overlap) * 3 + 2];
                        error2 = error2 + diff21 * diff21 + diff22 * diff22 + diff23 * diff23;
                        diff31 = in[(i + m + overlap) * inj + (j + k) * 3 + 0] - out[(row_out + m + overlap) * outj + (col_out + k) * 3 + 0];
                        diff32 = in[(i + m + overlap) * inj + (j + k) * 3 + 1] - out[(row_out + m + overlap) * outj + (col_out + k) * 3 + 1];
                        diff33 = in[(i + m + overlap) * inj + (j + k) * 3 + 2] - out[(row_out + m + overlap) * outj + (col_out + k) * 3 + 2];
                        error3 = error3 + diff31 * diff31 + diff32 * diff32 + diff33 * diff33;
                    }
                }
                errors[i * error_width + j] = error1 + error2 + error3;
            }
        }
    }
}

void calc_errors(image_t in, int blocksize, slice_t in_slice, slice_t out_slice, pixel_t *errors, int add)
{
    int error_jumpsize = in.width - blocksize + 1;

    if (add) // if add, the operation is plus (+)
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {
                in_slice.data = in.data + in_slice.jumpsize * i + j * in_slice.channels;
                errors[i * error_jumpsize + j] += l2norm(in_slice, out_slice);
            }
        }
    }
    else // the  operation is minus (-)
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {
                in_slice.data = in.data + in_slice.jumpsize * i + j * in_slice.channels;
                errors[i * error_jumpsize + j] -= l2norm(in_slice, out_slice);
            }
        }
    }
}

pixel_t l2norm(slice_t s1, slice_t s2)
{
    pixel_t error = 0;
    for (int i = 0; i < s1.height; i++)
    {
        for (int j = 0; j < s1.channels * s1.width; j++)
        {
            pixel_t s1_data = s1.data[i * s1.jumpsize + j];
            pixel_t s2_data = s2.data[i * s2.jumpsize + j];
            int diff = s1_data - s2_data;
            error += diff * diff;
        }
    }
    return error;
}

/*
The function find finds the coordinates of a candidate block that is within a certain range (specified by tolerance)
from the best fitting block
 */
coord find(pixel_t *errors, int height, int width, pixel_t tol_nom, pixel_t tol_den)
{
    pixel_t min_error = PIXEL_T_MAX;

    // search for the minimum error in the errors array and store it in a variable.
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            if (errors[i * width + j] < min_error)
                min_error = errors[i * width + j];

    // Count how many canditates exist in order to know the size of the array of candidates
    pixel_t tol_range = min_error + (min_error * tol_nom) / tol_den;
    int nr_candidates = 0;
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            if (errors[i * width + j] <= tol_range)
                nr_candidates++;

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
    int random_idx = rand() % nr_candidates;
    coord random_candidate = candidates[random_idx];

    // return coordinates of the random candidate
    return random_candidate;
}

// Command to compile the code:   gcc dpcut.c imageio.c L2norm.c  -o imageio -lm
// (n_blocks)(n_blocks-1) * 2 * (flops(calcerrors) + flops(dpcut)) + (n_blocks-1)^2 * flops(calcerrors) + (n_blocks^2-1) * find
// = (n_blocks)(n_blocks-1) * 2 * ((in_height - blocksize + 1)(in_width - blocksize + 1) * (1 + overlap*blocksize*3*3 + 1) + (overlap * blocksize) * 10 + overlap * 2 + (blocksize -1) * 5) + (n_blocks-1)^2 * (overlap*overlap*3*3+1) + (n_blocks^2-1) * ((in_height - blocksize + 1)(in_width - blocksize + 1)*3 +2)
// nb = n_blocks, ih = in_height, iw = in_width, ov = overlap, bs = blocksize
// = nb*(nb-1) * 2 * ((ih - bs + 1) * (iw -bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs -1) * 5) + (nb-1)^2 * (ov * ov * 3 * 3 + 1) + (nb * nb -1) * ((ih -bs + 1)*(iw - bs + 1) * 3 + 2)

image_t image_quilting(image_t in, int blocksize, int num_blocks, int overlap, pixel_t tol_nom, pixel_t tol_den)
{
    int out_size = num_blocks * blocksize - (num_blocks - 1) * overlap;
    int n = out_size * out_size * 3;
    pixel_t *out_image = (pixel_t *)calloc(n, sizeof(pixel_t));
    assert(out_image);
    image_t out = {out_image, out_size, out_size, 3};
    int errorlen = (in.height - blocksize + 1) * (in.width - blocksize + 1) * sizeof(pixel_t);
    pixel_t *errors = (pixel_t *)malloc(errorlen);
    assert(errors);

    for (int row = 0; row < num_blocks; row++)
    {
        for (int col = 0; col < num_blocks; col++)
        {
            int si = row * (blocksize - overlap);
            int sj = col * (blocksize - overlap);
            memset(errors, 0, errorlen);

            // Very first case, so pick one at random
            if (row == 0 && col == 0)
            {
                int row_idx = 0; // rand() % (in.height - blocksize + 1);
                int col_idx = 0; // rand() % (in.width - blocksize + 1);

                slice_t out_block = slice_image(out, si, sj, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, row_idx, col_idx, row_idx + blocksize, col_idx + blocksize);

                slice_cpy(in_block, out_block);
                continue;
            }
            calc_errors_(in, out, si, sj, blocksize, overlap, errors);
            // Search for a random candidate block among the best matching blocks (determined by the tolerance)
            coord random_candidate = find(errors, in.height - blocksize + 1, in.width - blocksize + 1, tol_nom, tol_den);
            int rand_row = random_candidate.row;
            int rand_col = random_candidate.col;
            if (row == 0)
            {
                slice_t out_block = slice_image(out, si, sj + overlap, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, rand_row, rand_col + overlap, rand_row + blocksize, rand_col + blocksize);
                slice_cpy(in_block, out_block);
            }
            else if (col == 0)
            {
                slice_t out_block = slice_image(out, si + overlap, sj, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, rand_row + overlap, rand_col, rand_row + blocksize, rand_col + blocksize);
                slice_cpy(in_block, out_block);
            }
            else
            {
                slice_t out_block = slice_image(out, si + overlap, sj + overlap, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, rand_row + overlap, rand_col + overlap, rand_row + blocksize, rand_col + blocksize);
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
        }
    }
    free(errors);
    return out;
}
// gcc dpcut.c imageio.c L2norm.c  -o imageio -lm
