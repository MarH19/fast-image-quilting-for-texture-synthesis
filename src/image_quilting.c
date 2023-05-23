#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_quilting.h"

#define INTEGRAL(start, height, width, jumpsize) (integral[start + height * jumpsize + width] - integral[start + width] - integral[start + height * jumpsize] + integral[start])
#define VEC32ADD(v) (v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7])

/*
orow: output row
ocol: output col
lrow: left predecessor row
lcol: left predecessor col
arow: above predecessor row
acol: above predecessor col
*/
void fill_error_matrix(image_t in_, image_t out_, int orow, int ocol, error_t *errors, error_t *integral, int lrow, int lcol, int arow, int acol, int blocksize, int overlap, int num_blocks)
{
    int error_width = in_.width - blocksize + 1;
    int error_height = in_.height - blocksize + 1;
    int integral_width = in_.width + 1;
    pixel_t *in = in_.data;
    pixel_t *out = out_.data;
    int ijump = in_.channels * in_.width;
    int ojump = out_.channels * out_.width;
    error_t diff0r, diff0g, diff0b, error0r, error0g, error0b;
    error_t diff1r, diff1g, diff1b, error1r, error1g, error1b;
    error_t diff2r, diff2g, diff2b, error2r, error2g, error2b;
    error_t diff3r, diff3g, diff3b, error3r, error3g, error3b;

    if (orow == 0)
    {
        /* we got 1 subblock to calculate
         * 3 . .
         * 3 . .
         * 3 . .
         */
        int block3_in_integral_base = 0;

        int height3 = blocksize;
        int width3 = overlap;
        // precalculation for block 3
        int block3_out_integral_start = lrow * integral_width + lcol + (blocksize - overlap);
        error_t block3_out_integral = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);
        __m256i zero = _mm256_setzero_si256();

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width; icol++)
            {
                error3r = error3g = error3b = 0;
                __m256i error00 = zero;
                error_t error00_arr[64 / sizeof(error_t)];

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                error_t block3_in_integral = INTEGRAL(block3_in_integral_start, height3, width3, integral_width);

                for (int k = 0; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        __m256i min00 = _mm256_loadu_si256((__m256i_u *) &in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i mout = _mm256_loadu_si256((__m256i_u *) &in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i min00_x_out = _mm256_madd_epi16(min00, mout);
                        error00 = _mm256_add_epi32(min00_x_out, error00);
                    }
                    for (; m < overlap * 3 - 7; m += 8)  // should run only once everything else is a bug
                    {
                        __m256i min00 = _mm256_loadu_si256((__m256i_u *) &in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i mout = _mm256_loadu_si256((__m256i_u *) &in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i min00_x_out = _mm256_madd_epi16(min00, mout);
                        error00 = _mm256_add_epi32(_mm256_blend_epi32(min00_x_out, zero, 0b11110000), error00);
                    }
                }
                _mm256_storeu_si256((__m256i_u *) error00_arr, error00);
                errors[irow * error_width + icol] = block3_out_integral - 2 * VEC32ADD(error00_arr) + block3_in_integral;
            }
        }
    }
    else if (ocol == 0)
    {
        /* we got 4 subblocks to calculate
         * 1 1 2
         * . . .
         * . . .
         */
        int block1_in_integral_base = 0;

        int height1 = overlap;
        int width1 = blocksize - overlap;

        // precalculation for block 1
        int block1_out_integral_start = (arow + blocksize - overlap) * integral_width + acol;
        error_t block1_out_integral = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width; icol++)
            {
                error1r = error2r = error1g = error2g = error1b = error2b = 0;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                error_t block1_in_integral = INTEGRAL(block1_in_integral_start, height1, width1, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < overlap; m++)
                    {
                        // block2 l2norm
                        diff2r = (error_t)in[(irow + k) * ijump + (icol + m + blocksize - overlap) * 3 + 0] - (error_t)out[(orow + k) * ojump + (ocol + m + blocksize - overlap) * 3 + 0];
                        diff2g = (error_t)in[(irow + k) * ijump + (icol + m + blocksize - overlap) * 3 + 1] - (error_t)out[(orow + k) * ojump + (ocol + m + blocksize - overlap) * 3 + 1];
                        diff2b = (error_t)in[(irow + k) * ijump + (icol + m + blocksize - overlap) * 3 + 2] - (error_t)out[(orow + k) * ojump + (ocol + m + blocksize - overlap) * 3 + 2];
                        error2r += diff2r * diff2r;
                        error2g += diff2g * diff2g;
                        error2b += diff2b * diff2b;
                    }
                    for (int m = 0; m < (blocksize - overlap); m++)
                    {
                        // block 1 -> mul_sum
                        diff1r = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 0] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m) * 3 + 0];
                        diff1g = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 1] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m) * 3 + 1];
                        diff1b = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 2] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m) * 3 + 2];
                        error1r += diff1r;
                        error1g += diff1g;
                        error1b += diff1b;
                    }
                }
                errors[irow * error_width + icol] = error2r + error2g + error2b + block1_out_integral - 2 * (error1r + error1g + error1b) + block1_in_integral;
            }
        }
    }
    else if (ocol / (blocksize - overlap) != num_blocks - 1)
    {
        /* we got 4 subblocks to calculate
         * 0 1 2
         * 3 . .
         * 3 . .
         */
        int block1_in_integral_base = overlap;
        int block3_in_integral_base = overlap * integral_width;

        int height1 = overlap;
        int width1 = blocksize - 2 * overlap;

        int height3 = blocksize - overlap;
        int width3 = overlap;
        // precalculation for block 1
        int block1_out_integral_start = (arow + blocksize - overlap) * integral_width + acol + overlap;
        error_t block1_out_integral = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        // precalculation for block 3
        int block3_out_integral_start = (lrow + overlap) * integral_width + lcol + (blocksize - overlap);
        error_t block3_out_integral = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width; icol++)
            {
                error0r = error1r = error2r = error3r = error0g = error1g = error2g = error3g = error0b = error1b = error2b = error3b = 0;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                error_t block1_in_integral = INTEGRAL(block1_in_integral_start, height1, width1, integral_width);

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                error_t block3_in_integral = INTEGRAL(block3_in_integral_start, height3, width3, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < overlap; m++)
                    {
                        // block0 l2norm
                        diff0r = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 0] - (error_t)out[(orow + k) * ojump + (ocol + m) * 3 + 0];
                        diff0g = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 1] - (error_t)out[(orow + k) * ojump + (ocol + m) * 3 + 1];
                        diff0b = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 2] - (error_t)out[(orow + k) * ojump + (ocol + m) * 3 + 2];
                        error0r += diff0r * diff0r;
                        error0g += diff0g * diff0g;
                        error0b += diff0b * diff0b;

                        // block2 l2norm
                        diff2r = (error_t)in[(irow + k) * ijump + (icol + m + blocksize - overlap) * 3 + 0] - (error_t)out[(orow + k) * ojump + (ocol + m + blocksize - overlap) * 3 + 0];
                        diff2g = (error_t)in[(irow + k) * ijump + (icol + m + blocksize - overlap) * 3 + 1] - (error_t)out[(orow + k) * ojump + (ocol + m + blocksize - overlap) * 3 + 1];
                        diff2b = (error_t)in[(irow + k) * ijump + (icol + m + blocksize - overlap) * 3 + 2] - (error_t)out[(orow + k) * ojump + (ocol + m + blocksize - overlap) * 3 + 2];
                        error2r += diff2r * diff2r;
                        error2g += diff2g * diff2g;
                        error2b += diff2b * diff2b;
                    }
                    for (int m = 0; m < (blocksize - 2 * overlap); m++)
                    {
                        // block 1 -> mul_sum
                        diff1r = (error_t)in[(irow + k) * ijump + (icol + m + overlap) * 3 + 0] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m + overlap) * 3 + 0];
                        diff1g = (error_t)in[(irow + k) * ijump + (icol + m + overlap) * 3 + 1] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m + overlap) * 3 + 1];
                        diff1b = (error_t)in[(irow + k) * ijump + (icol + m + overlap) * 3 + 2] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m + overlap) * 3 + 2];
                        error1r += diff1r;
                        error1g += diff1g;
                        error1b += diff1b;
                    }
                }
                for (int k = 0; k < blocksize - overlap; k++)
                {
                    for (int m = 0; m < overlap; m++)
                    {
                        // block3 mul_sum part 1
                        diff3r = (error_t)in[(irow + k + overlap) * ijump + (icol + m) * 3 + 0] * (error_t)in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap + m) * 3 + 0];
                        diff3g = (error_t)in[(irow + k + overlap) * ijump + (icol + m) * 3 + 1] * (error_t)in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap + m) * 3 + 1];
                        diff3b = (error_t)in[(irow + k + overlap) * ijump + (icol + m) * 3 + 2] * (error_t)in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap + m) * 3 + 2];
                        error3r += diff3r;
                        error3g += diff3g;
                        error3b += diff3b;
                    }
                }
                errors[irow * error_width + icol] = error0r + error0g + error0b + error2r + error2g + error2b + block1_out_integral - 2 * (error1r + error1g + error1b) + block1_in_integral + block3_out_integral - 2 * (error3r + error3g + error3b) + block3_in_integral;
            }
        }
    }
    else
    {
        /* we got 4 subblocks to calculate
         * 0 1 1
         * 3 . .
         * 3 . .
         */
        int block1_in_integral_base = overlap;
        int block3_in_integral_base = overlap * integral_width;

        int height1 = overlap;
        int width1 = blocksize - overlap;

        int height3 = blocksize - overlap;
        int width3 = overlap;
        // precalculation for block 1
        int block1_out_integral_start = (arow + blocksize - overlap) * integral_width + acol + overlap;
        error_t block1_out_integral = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        // precalculation for block 3
        int block3_out_integral_start = (lrow + overlap) * integral_width + lcol + (blocksize - overlap);
        error_t block3_out_integral = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width; icol++)
            {
                error0r = error1r = error3r = error0g = error1g = error3g = error0b = error1b = error3b = 0;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                error_t block1_in_integral = INTEGRAL(block1_in_integral_start, height1, width1, integral_width);

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                error_t block3_in_integral = INTEGRAL(block3_in_integral_start, height3, width3, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < overlap; m++)
                    {
                        // block0 l2norm
                        diff0r = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 0] - (error_t)out[(orow + k) * ojump + (ocol + m) * 3 + 0];
                        diff0g = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 1] - (error_t)out[(orow + k) * ojump + (ocol + m) * 3 + 1];
                        diff0b = (error_t)in[(irow + k) * ijump + (icol + m) * 3 + 2] - (error_t)out[(orow + k) * ojump + (ocol + m) * 3 + 2];
                        error0r += diff0r * diff0r;
                        error0g += diff0g * diff0g;
                        error0b += diff0b * diff0b;
                    }
                    for (int m = 0; m < (blocksize - overlap); m++)
                    {
                        // block 1 -> mul_sum
                        diff1r = (error_t)in[(irow + k) * ijump + (icol + m + overlap) * 3 + 0] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m + overlap) * 3 + 0];
                        diff1g = (error_t)in[(irow + k) * ijump + (icol + m + overlap) * 3 + 1] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m + overlap) * 3 + 1];
                        diff1b = (error_t)in[(irow + k) * ijump + (icol + m + overlap) * 3 + 2] * (error_t)in[(arow + (blocksize - overlap) + k) * ijump + (acol + m + overlap) * 3 + 2];
                        error1r += diff1r;
                        error1g += diff1g;
                        error1b += diff1b;
                    }
                }
                for (int k = 0; k < blocksize - overlap; k++)
                {
                    for (int m = 0; m < overlap; m++)
                    {
                        // block3 mul_sum part 1
                        diff3r = (error_t)in[(irow + k + overlap) * ijump + (icol + m) * 3 + 0] * (error_t)in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap + m) * 3 + 0];
                        diff3g = (error_t)in[(irow + k + overlap) * ijump + (icol + m) * 3 + 1] * (error_t)in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap + m) * 3 + 1];
                        diff3b = (error_t)in[(irow + k + overlap) * ijump + (icol + m) * 3 + 2] * (error_t)in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap + m) * 3 + 2];
                        error3r += diff3r;
                        error3g += diff3g;
                        error3b += diff3b;
                    }
                }
                errors[irow * error_width + icol] = error0r + error0g + error0b + block1_out_integral - 2 * (error1r + error1g + error1b) + block1_in_integral + block3_out_integral - 2 * (error3r + error3g + error3b) + block3_in_integral;
            }
        }
    }
}

/* Here the integral image of the input image is calculated, where the left and upper border consists of zeros to facilitate the calculation of the integral */
void generate_integral_image(image_t in, error_t *out)
{
    int outjump = in.width + 1;
    assert(in.channels == 3);
    int injump = in.width * 3;
    for (int i = 0; i < in.height + 1; i++)
        out[i * outjump] = 0;
    for (int i = 0; i < outjump; i++)
        out[i] = 0;
    for (int i = 1; i < in.height + 1; i++)
    {
        for (int j = 1; j < in.width + 1; j++)
        {
            error_t temp = out[i * outjump + j - 1] + out[(i - 1) * outjump + j] - out[(i - 1) * outjump + j - 1];
            // offset -1 in row and col as we start with (1, 1)
            temp += (error_t)in.data[(i - 1) * injump + 3 * j - 3] * (error_t)in.data[(i - 1) * injump + 3 * j - 3];
            temp += (error_t)in.data[(i - 1) * injump + 3 * j - 2] * (error_t)in.data[(i - 1) * injump + 3 * j - 2];
            temp += (error_t)in.data[(i - 1) * injump + 3 * j - 1] * (error_t)in.data[(i - 1) * injump + 3 * j - 1];
            out[i * outjump + j] = temp;
        }
    }
}

/*
The function find finds the coordinates of a candidate block that is within a certain range (specified by tolerance)
from the best fitting block
 */
coord find(error_t *errors, int height, int width, pixel_t tol_nom, pixel_t tol_den)
{
    error_t min_error = ERROR_T_MAX;

    // search for the minimum error in the errors array and store it in a variable.
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            if (errors[i * width + j] < min_error)
                min_error = errors[i * width + j];

    // Count how many canditates exist in order to know the size of the array of candidates
    error_t tol_range = min_error + (min_error * tol_nom) / tol_den;
    // Create the array with maximum amount of possible candidates
    coord candidates[height * width];
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
    int random_idx = rand() % idx;
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
    int errorlen = (in.height - blocksize + 1) * (in.width - blocksize + 1);
    error_t errors[errorlen];

    // Integral image pads with 0 by 1 row and 1 col to simplify calculations
    error_t integral[(in.height + 1) * (in.width + 1)];

    generate_integral_image(in, integral);

    /* These arrays are needed to determine where in the input image the block above the one to be
       inserted in the output image comes from in order to calculate the integral */
    int row_above_list[num_blocks];
    int col_above_list[num_blocks];

    /* Here we store the indices of the block that lies above the block we want to insert,
    so that we can find it in the input image and calculate the integral */
    int row_above;
    int col_above;

    /* Here we store the indices of the block that is to the left of the block we want to insert,
    so that we can find it in the input image and calculate the integral */
    int row_left;
    int col_left;

    for (int row = 0; row < num_blocks; row++)
    {
        for (int col = 0; col < num_blocks; col++)
        {

            int si = row * (blocksize - overlap);
            int sj = col * (blocksize - overlap);

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
            row_above = row_above_list[col];
            col_above = col_above_list[col];
            fill_error_matrix(in, out, si, sj, errors, integral, row_left, col_left, row_above, col_above, blocksize, overlap, num_blocks);

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
    return out;
}
