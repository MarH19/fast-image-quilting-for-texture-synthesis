#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_quilting.h"

#define INTEGRAL(start, height, width, jumpsize) (integral[start + height * jumpsize + width] - integral[start + width] - integral[start + height * jumpsize] + integral[start])
#define VEC32ADD(v) (v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7])
/* reference: https://stackoverflow.com/questions/60108658/fastest-method-to-calculate-sum-of-all-packed-32-bit-integers-using-avx512-or-av */
// static inline
static inline uint32_t hsum_epi32_avx(__m128i x)
{
    __m128i hi64  = _mm_unpackhi_epi64(x, x); // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
    __m128i sum64 = _mm_add_epi32(hi64, x);
    __m128i hi32  = _mm_shuffle_epi32(sum64, _MM_SHUFFLE(2, 3, 0, 1)); // Swap the low two elements
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32);       // movd
}

// only needs AVX2
static inline uint32_t hsum_8x32(__m256i v)
{
    __m128i sum128 = _mm_add_epi32( 
                 _mm256_castsi256_si128(v),
                 _mm256_extracti128_si256(v, 1));
    return hsum_epi32_avx(sum128);
}

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
    __m256i zero = _mm256_setzero_si256();

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

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width; icol++)
            {
                __m256i error00 = zero;

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                error_t block3_in_integral = INTEGRAL(block3_in_integral_start, height3, width3, integral_width);

                for (int k = 0; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        error00 = _mm256_add_epi32(mull00, error00);
                    }
                    if (overlap % 16 != 0)
                    {
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        error00 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00);
                    }
                }
                errors[irow * error_width + icol] = block3_out_integral - 2 * hsum_8x32(error00) + block3_in_integral;
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
                __m256i error00, error00_block1, error00_block2, error00_block02, error00_block13;
                error00 = error00_block1 = error00_block2 = error00_block02 = error00_block13 = zero;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                error_t block1_in_integral = INTEGRAL(block1_in_integral_start, height1, width1, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < (blocksize - overlap) * 3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block1 mulsum part3
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b11110000));

                        // block2 l2norm part1
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_blend_epi32(sqrdiff00, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < blocksize * 3 - 15; m += 16)
                    {
                        // block2 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_madd_epi16(diff00, diff00));
                    }
                }

                error00_block02 = error00_block2;
                error00_block13 = _mm256_slli_epi32(error00_block1, 1);
                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                errors[irow * error_width + icol] = block1_out_integral + block1_in_integral + hsum_8x32(error00);
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
                __m256i error00, error00_block0, error00_block1, error00_block2, error00_block3, error00_block02, error00_block13;
                error00 = error00_block0 = error00_block1 = error00_block2 = error00_block3 = error00_block02 = error00_block13 = zero;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                error_t block1_in_integral = INTEGRAL(block1_in_integral_start, height1, width1, integral_width);

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                error_t block3_in_integral = INTEGRAL(block3_in_integral_start, height3, width3, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block0 l2norm part1
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_madd_epi16(diff00, diff00));
                    }
                    if (overlap % 16 != 0)
                    {
                        // block0 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_blend_epi32(sqrdiff00, zero, 0b11110000));

                        // block1 mulsum part1
                        // above_in = above-block input number 00
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < (blocksize - overlap) * 3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block1 mulsum part3
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b11110000));

                        // block2 l2norm part1
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_blend_epi32(sqrdiff00, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < blocksize * 3 - 15; m += 16)
                    {
                        // block2 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_madd_epi16(diff00, diff00));
                    }
                }
                // we start with row overlap to ensure equivalence in handling block 3 (to other cases)
                for (int k = overlap; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block3 mulsum part1
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        error00_block3 = _mm256_add_epi32(mull00, error00_block3);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block3 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        error00_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00_block3);
                    }
                }

                error00_block02 = _mm256_add_epi32(error00_block0, error00_block2);
                error00_block13 = _mm256_slli_epi32(_mm256_add_epi32(error00_block1, error00_block3), 1);
                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                errors[irow * error_width + icol] = block1_out_integral + block1_in_integral + hsum_8x32(error00) + block3_out_integral + block3_in_integral;
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
                __m256i error00, error00_block0, error00_block1, error00_block3, error00_block02, error00_block13;
                error00 = error00_block0 = error00_block1 = error00_block3 = error00_block02 = error00_block13 = zero;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                error_t block1_in_integral = INTEGRAL(block1_in_integral_start, height1, width1, integral_width);

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                error_t block3_in_integral = INTEGRAL(block3_in_integral_start, height3, width3, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block0 l2norm part1
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_madd_epi16(diff00, diff00));
                    }
                    if (overlap % 16 != 0)
                    {
                        // block0 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_blend_epi32(sqrdiff00, zero, 0b11110000));

                        // block1 mulsum part1
                        // above_in = above-block input number 00
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < (blocksize)*3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol)*3 + m]);
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                    }
                }
                // we start with row overlap to ensure equivalence in handling block 3 (to other cases)
                for (int k = overlap; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block3 mulsum part1
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        error00_block3 = _mm256_add_epi32(mull00, error00_block3);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block3 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        error00_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00_block3);
                    }
                }

                error00_block02 = error00_block0;
                error00_block13 = _mm256_slli_epi32(_mm256_add_epi32(error00_block1, error00_block3), 1);
                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                errors[irow * error_width + icol] = block1_out_integral + block1_in_integral + hsum_8x32(error00) + block3_out_integral + block3_in_integral;
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
