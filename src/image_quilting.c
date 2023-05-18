#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_quilting.h"




/* Here the integral image of the input image is calculated, where the left and upper border consists of zeros to facilitate the calculation of the integral */
void generate_integral_image(image_t in, pixel_t *s2)
{

    int s2_jumpsize = in.width + 1;
    int in_jumpsize = in.width * in.channels;

    for (int i = 1; i < in.height + 1; i++)
    {
        for (int j = 1; j < in.width + 1; j++)
        {
            //printf("i: %d, j: %d\n", i, j);
            s2[i * s2_jumpsize + j] = s2[i * s2_jumpsize + j - 1] + s2[(i - 1) * s2_jumpsize + j] - s2[(i - 1) * s2_jumpsize + j - 1] + in.data[(i - 1) * in_jumpsize + in.channels * j - 3] * in.data[(i - 1) * in_jumpsize + in.channels * j - 3] +
                                      in.data[(i - 1) * in_jumpsize + in.channels * j-2] * in.data[(i - 1) * in_jumpsize + in.channels * j-2] +
                                      in.data[(i - 1) * in_jumpsize + in.channels * j-1] * in.data[(i - 1) * in_jumpsize + in.channels * j-1];

            /*
            printf("out: %d, left: %d, top: %d, angle: %d, p1: %d, p2: %d, p3: %d\n",s2[i * s2_jumpsize + j], s2[i * s2_jumpsize + j - 1], s2[(i - 1) * s2_jumpsize + j], s2[(i - 1) * s2_jumpsize + j - 1],
            in.data[(i - 1) * in_jumpsize + in.channels * j - 3], in.data[(i - 1) * in_jumpsize + in.channels * j-2], in.data[(i - 1) * in_jumpsize + in.channels * j-1]);
            if(j==3)exit(1);
            */
        }
    }
}

/*
Calculates the error between the slice of an input block (i.e., for each block of the input image) and the slice of
a fixed output block
*/
void calc_errors(image_t in, int blocksize, slice_t in_slice, slice_t out_slice, pixel_t *errors, int add)
{
    int error_jumpsize = in.width - blocksize + 1;
    pixel_t *in_slice_start = in_slice.data;

    if (add) // if add, the operation is plus (+)
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {
                in_slice.data = in_slice_start + in_slice.jumpsize * i + j * in_slice.channels;
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
    pixel_t tol_range = min_error + (min_error / tol_den) * tol_nom;
    int nr_candidates = 0;
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
    int random_idx = rand() % nr_candidates;
    coord random_candidate = candidates[random_idx];

    // return coordinates of the random candidate
    return random_candidate;
}

/* Sum of in_overlap(x,y) * out_overlap(x,y) */
pixel_t mul_in_out(slice_t s1, slice_t s2)
{
    pixel_t error = 0;
    assert(s1.width == s2.width);
    assert(s1.height == s2.height);
    assert(s1.channels == s2.channels);

    for (int i = 0; i < s1.height; i++)
    {
        for (int j = 0; j < s1.channels * s1.width; j++)
        {
            pixel_t s1_data = s1.data[i * s1.jumpsize + j];
            pixel_t s2_data = s2.data[i * s2.jumpsize + j];
            error += s1_data * s2_data;
        }
    }
    return error;
}

/* This function is used to calculate (a part of) the overlap error using the precomputed integral image */
void integral_sum(image_t in, int blocksize, int overlap, slice_t in_slice, int in_integral_start, slice_t out_slice, int out_integral_start, pixel_t *errors, pixel_t *integral)
{
    // Compute the integral for the overlap region of the output block
    int integral_jumpsize = in.width + 1;
    pixel_t out_integral = integral[out_integral_start + out_slice.height * integral_jumpsize + out_slice.width] 
                            - integral[out_integral_start + out_slice.width] 
                            - integral[out_integral_start + out_slice.height * integral_jumpsize] 
                            + integral[out_integral_start];

    
    int error_jumpsize = in.width - blocksize + 1;
    pixel_t *in_slice_start = in_slice.data;

    for (int i = 0; i < in.height - blocksize + 1; i++)
    {
        for (int j = 0; j < in.width - blocksize + 1; j++)
        {   
            // Compute integral for the overlap region of the input block(s)
            int new_start = in_integral_start + i*integral_jumpsize + j;
            pixel_t in_integral = integral[new_start + in_slice.height * integral_jumpsize + in_slice.width] 
                            - integral[new_start + in_slice.width] 
                            - integral[new_start + in_slice.height * integral_jumpsize] 
                            + integral[new_start];


            in_slice.data = in_slice_start + in_slice.jumpsize * i + j * in_slice.channels;
            // Calculate (a part of) the overlap error for an input block
            errors[i * error_jumpsize + j] += (in_integral - 2 * mul_in_out(in_slice, out_slice) + out_integral);
        }
    }
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

    // Initialize the integral image with size (in.width + 1) * (in.height + 1) in order to have later zeros on the left and upper border
    pixel_t *integral = (pixel_t *)malloc((in.height + 1) * (in.width + 1) * sizeof(pixel_t));
    assert(integral);
    memset(integral, 0, (in.height + 1) * (in.width + 1) * sizeof(pixel_t));

    // Precompute the integral image
    generate_integral_image(in, integral);

    /* These arrays are needed to determine where in the input image the block above the one to be 
       inserted in the output image comes from in order to calculate the integral */
    int *row_above_list = (int *)malloc(num_blocks * sizeof(int));
    assert(row_above_list);
    int *col_above_list = (int *)malloc(num_blocks * sizeof(int));
    assert(col_above_list);

    memset(row_above_list, 0, num_blocks * sizeof(int));
    memset(col_above_list, 0, num_blocks * sizeof(int));

    /* Here we store the indices of the block that lies above the block we want to insert, 
    so that we can find it in the input image and calculate the integral */
    int row_above = 0;
    int col_above = 0;

    /* Here we store the indices of the block that is to the left of the block we want to insert, 
    so that we can find it in the input image and calculate the integral */
    int row_left = 0;
    int col_left = 0;

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

                // update the indices with the indices of the inserted block
                row_left = row_idx;
                col_left = col_idx;
                row_above_list[col] = row_idx;
                col_above_list[col] = col_idx;

                slice_t out_block = slice_image(out, si, sj, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, row_idx, col_idx, row_idx + blocksize, col_idx + blocksize);

                slice_cpy(in_block, out_block);
                continue;
            }
            /* In this case we consider the top row of the output image and the overlap 
            error can be fully calculated using the formula that contains the integral image part */
            if (row == 0 && col != 0)
            {
                slice_t out_slice = slice_image(in, row_left, col_left + (blocksize - overlap), row_left + blocksize, col_left + blocksize);
                slice_t in_slice = slice_image(in, 0, 0, blocksize, overlap);
                int integral_width = in.width + 1;
                int in_integral_start = 0;
                int out_integral_start = row_left * integral_width + (col_left + (blocksize - overlap));
                integral_sum(in, blocksize, overlap, in_slice, in_integral_start, out_slice, out_integral_start, errors, integral);
            }
            /* In this case we consider the left column of the output image starting from the second row */
            if (row != 0 && col == 0)
            {
                // A left part of the upper overlap can be calculated using integral images
                row_above = row_above_list[col];
                col_above = col_above_list[col];
                slice_t out_slice = slice_image(in, row_above + (blocksize - overlap), col_above, row_above + blocksize, col_above + (blocksize - overlap));
                slice_t in_slice = slice_image(in, 0, 0, overlap, blocksize - overlap);
                int integral_width = in.width + 1;
                int in_integral_start = 0;
                int out_integral_start = row_above * integral_width + col_above + integral_width * (blocksize - overlap);
                integral_sum(in, blocksize, overlap, in_slice, in_integral_start, out_slice, out_integral_start, errors, integral);

                // The top right square (of size overlap * overlap) of the upper overlap can be calculated using calc_errors
                out_slice = slice_image(out, si, sj + (blocksize - overlap) ,si + overlap, sj + blocksize);
                in_slice = slice_image(in, 0, blocksize - overlap, overlap, blocksize);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }
            /* In this case we consider all blocks that are between the top row, the left column and the right column */
            if (row != 0 && col != 0 && col != num_blocks - 1)
            {   
                // The middle part of the upper overlap can be calculated using integral images
                row_above = row_above_list[col];
                col_above = col_above_list[col];
                slice_t out_slice = slice_image(in, row_above + (blocksize - overlap), col_above + overlap, row_above + blocksize, col_above + (blocksize - overlap));
                slice_t in_slice = slice_image(in, 0, overlap, overlap, blocksize - overlap);
                int integral_width = in.width + 1;
                int in_integral_start = overlap;
                int out_integral_start = row_above * integral_width + col_above + integral_width * (blocksize - overlap) + overlap;
                integral_sum(in, blocksize, overlap, in_slice, in_integral_start, out_slice, out_integral_start, errors, integral);

                // The lower part of the left overlap can be calculated using integral images
                out_slice = slice_image(in, row_left + overlap, col_left + (blocksize - overlap), row_left + blocksize, col_left + blocksize);
                in_slice = slice_image(in, overlap, 0, blocksize, overlap);
                integral_width = in.width + 1;
                in_integral_start = integral_width * overlap;
                out_integral_start = row_left * integral_width + col_left + integral_width * overlap + (blocksize - overlap);
                integral_sum(in, blocksize, overlap, in_slice, in_integral_start, out_slice, out_integral_start, errors, integral);

                // The right square part of the upper overlap can be calculated using calc_errors
                out_slice = slice_image(out, si, sj + (blocksize - overlap) ,si + overlap, sj + blocksize);
                in_slice = slice_image(in, 0, blocksize - overlap, overlap, blocksize);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);

                // The left square part of the upper overlap can be calculated using calc_errors
                out_slice = slice_image(out, si, sj ,si + overlap, sj + overlap);
                in_slice = slice_image(in, 0, 0, overlap, overlap);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }
            /* In this case we consider all blocks that are on the right column */
            if (row != 0 && col == num_blocks - 1)
            {
                // The right part of the upper overlap can be calculated using integral images
                row_above = row_above_list[col];
                col_above = col_above_list[col];
                slice_t out_slice = slice_image(in, row_above + (blocksize - overlap), col_above + overlap, row_above + blocksize, col_above + blocksize);
                slice_t in_slice = slice_image(in, 0, overlap, overlap, blocksize);
                int integral_width = in.width + 1;
                int in_integral_start = overlap;
                int out_integral_start = row_above * integral_width + col_above + integral_width * (blocksize - overlap) + overlap;
                integral_sum(in, blocksize, overlap, in_slice, in_integral_start, out_slice, out_integral_start, errors, integral);

                // The lower part of the overlap on the left can be calculated using integral images
                out_slice = slice_image(in, row_left + overlap, col_left + (blocksize - overlap), row_left + blocksize, col_left + blocksize);
                in_slice = slice_image(in, overlap, 0, blocksize, overlap);
                integral_width = in.width + 1;
                in_integral_start = integral_width * overlap;
                out_integral_start = row_left * integral_width + col_left + integral_width * overlap + (blocksize - overlap);
                integral_sum(in, blocksize, overlap, in_slice, in_integral_start, out_slice, out_integral_start, errors, integral);

                // The left square of the upper overlap can be calculated using calc_errors
                out_slice = slice_image(out, si, sj ,si + overlap, sj + overlap);
                in_slice = slice_image(in, 0, 0, overlap, overlap);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }

            // Search for a random candidate block among the best matching blocks (determined by the tolerance)
            coord random_candidate = find(errors, in.height - blocksize + 1, in.width - blocksize + 1, tol_nom, tol_den);
            int rand_row = random_candidate.row;
            int rand_col = random_candidate.col;

            // update the indices with the indices of the inserted block
            row_left = rand_row;
            col_left = rand_col;
            row_above_list[col] = rand_row;
            col_above_list[col] = rand_col;

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
    free(integral);
    free(row_above_list);
    free(col_above_list);
    return out;
}
// gcc dpcut.c imageio.c L2norm.c  -o imageio -lm
