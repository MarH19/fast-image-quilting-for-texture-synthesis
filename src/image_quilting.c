#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_quilting.h"

#define INTEGRAL(start, height, width, jumpsize) (integral[start + height * jumpsize + width] - integral[start + width] - integral[start + height * jumpsize] + integral[start])
#define VEC32ADD(v) (v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7])

// reference https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
static inline uint32_t hsum_epi32_avx(__m128i x)
{
    __m128i hi64 = _mm_unpackhi_epi64(x, x); // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
    __m128i sum64 = _mm_add_epi32(hi64, x);
    __m128i hi32 = _mm_shuffle_epi32(sum64, _MM_SHUFFLE(2, 3, 0, 1)); // Swap the low two elements
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32); // movd
}

// only needs AVX2
static inline uint32_t hsum_8x32(__m256i v)
{
    __m128i sum128 = _mm_add_epi32(
        _mm256_castsi256_si128(v),
        _mm256_extracti128_si256(v, 1));
    return hsum_epi32_avx(sum128);
}

static inline __m256i loadu_8x16(pixel_t *v)
{
    return _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i_u *)v));
}

// reference https://stackoverflow.com/questions/51274287/computing-8-horizontal-sums-of-eight-avx-single-precision-floating-point-vectors
static inline __m256i hsum_8x8x32(__m256i v0, __m256i v1, __m256i v2, __m256i v3, __m256i v4, __m256i v5, __m256i v6, __m256i v7)
{
    const __m256i s01 = _mm256_hadd_epi32(v0, v1);
    const __m256i s23 = _mm256_hadd_epi32(v2, v3);
    const __m256i s45 = _mm256_hadd_epi32(v4, v5);
    const __m256i s67 = _mm256_hadd_epi32(v6, v7);
    const __m256i s0123 = _mm256_hadd_epi32(s01, s23);
    const __m256i s4567 = _mm256_hadd_epi32(s45, s67);

    // inter-lane shuffle
    v0 = _mm256_blend_epi32(s0123, s4567, 0xF0);
    v1 = _mm256_permute2f128_si256(s0123, s4567, 0x21);

    return _mm256_add_epi32(v0, v1);
}

// #define INTEGRAL(start, height, width, jumpsize) (integral[start + height * jumpsize + width] - integral[start + width] - integral[start + height * jumpsize] + integral[start])
static inline __m256i integralx8(int start, int height, int width, int jumpsize, error_t *integral)
{
    /* 0 . . 1
       . . . .
       2 . . 3
    */
    const __m256i v0 = _mm256_loadu_si256((__m256i_u *)&integral[start]);
    const __m256i v1 = _mm256_loadu_si256((__m256i_u *)&integral[start + width]);
    const __m256i v2 = _mm256_loadu_si256((__m256i_u *)&integral[start + height * jumpsize]);
    const __m256i v3 = _mm256_loadu_si256((__m256i_u *)&integral[start + height * jumpsize + width]);
    return _mm256_sub_epi32(_mm256_add_epi32(v3, v0), _mm256_add_epi32(v2, v1));
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
        error_t block3_out_integral_scalar = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);
        __m256i block3_out_integral = _mm256_set1_epi32(block3_out_integral_scalar);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width - 7; icol += 8)
            {
                __m256i error00 = zero;
                __m256i error01 = zero;
                __m256i error02 = zero;
                __m256i error03 = zero;
                __m256i error04 = zero;
                __m256i error05 = zero;
                __m256i error06 = zero;
                __m256i error07 = zero;

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                __m256i block3_in_integral = integralx8(block3_in_integral_start, height3, width3, integral_width, integral);
                __m256i block_integral = _mm256_add_epi32(block3_in_integral, block3_out_integral);

                for (int k = 0; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);

                        error00 = _mm256_add_epi32(mull00, error00);
                        error01 = _mm256_add_epi32(mull01, error01);
                        error02 = _mm256_add_epi32(mull02, error02);
                        error03 = _mm256_add_epi32(mull03, error03);
                        error04 = _mm256_add_epi32(mull04, error04);
                        error05 = _mm256_add_epi32(mull05, error05);
                        error06 = _mm256_add_epi32(mull06, error06);
                        error07 = _mm256_add_epi32(mull07, error07);
                    }
                    if (overlap % 16 != 0)
                    {
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);

                        error00 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00);
                        error01 = _mm256_add_epi32(_mm256_blend_epi32(mull01, zero, 0b11110000), error01);
                        error02 = _mm256_add_epi32(_mm256_blend_epi32(mull02, zero, 0b11110000), error02);
                        error03 = _mm256_add_epi32(_mm256_blend_epi32(mull03, zero, 0b11110000), error03);
                        error04 = _mm256_add_epi32(_mm256_blend_epi32(mull04, zero, 0b11110000), error04);
                        error05 = _mm256_add_epi32(_mm256_blend_epi32(mull05, zero, 0b11110000), error05);
                        error06 = _mm256_add_epi32(_mm256_blend_epi32(mull06, zero, 0b11110000), error06);
                        error07 = _mm256_add_epi32(_mm256_blend_epi32(mull07, zero, 0b11110000), error07);
                    }
                }
                //__m256i error = _mm256_set_epi32(hsum_8x32(error07), hsum_8x32(error06), hsum_8x32(error05), hsum_8x32(error04), hsum_8x32(error03), hsum_8x32(error02), hsum_8x32(error01), hsum_8x32(error00));
                __m256i error = hsum_8x8x32(error00, error01, error02, error03, error04, error05, error06, error07);
                __m256i twice = _mm256_slli_epi32(error, 1);
                __m256i error_with_integral = _mm256_sub_epi32(block_integral, twice);
                _mm256_storeu_si256((__m256i_u *)&errors[irow * error_width + icol], error_with_integral);
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
        error_t block1_out_integral_scalar = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        __m256i block1_out_integral = _mm256_set1_epi32(block1_out_integral_scalar);
        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width - 7; icol += 8)
            {
                __m256i error00, error00_block1, error00_block2, error00_block02, error00_block13;
                __m256i error01, error01_block1, error01_block2, error01_block02, error01_block13;
                __m256i error02, error02_block1, error02_block2, error02_block02, error02_block13;
                __m256i error03, error03_block1, error03_block2, error03_block02, error03_block13;
                __m256i error04, error04_block1, error04_block2, error04_block02, error04_block13;
                __m256i error05, error05_block1, error05_block2, error05_block02, error05_block13;
                __m256i error06, error06_block1, error06_block2, error06_block02, error06_block13;
                __m256i error07, error07_block1, error07_block2, error07_block02, error07_block13;

                error00 = error00_block1 = error00_block2 = error00_block02 = error00_block13 = zero;
                error01 = error01_block1 = error01_block2 = error01_block02 = error01_block13 = zero;
                error02 = error02_block1 = error02_block2 = error02_block02 = error02_block13 = zero;
                error03 = error03_block1 = error03_block2 = error03_block02 = error03_block13 = zero;
                error04 = error04_block1 = error04_block2 = error04_block02 = error04_block13 = zero;
                error05 = error05_block1 = error05_block2 = error05_block02 = error05_block13 = zero;
                error06 = error06_block1 = error06_block2 = error06_block02 = error06_block13 = zero;
                error07 = error07_block1 = error07_block2 = error07_block02 = error07_block13 = zero;

                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                __m256i block1_in_integral = integralx8(block1_in_integral_start, height1, width1, integral_width, integral);
                __m256i block_integral = _mm256_add_epi32(block1_in_integral, block1_out_integral);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < (blocksize - overlap) * 3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                        error01_block1 = _mm256_add_epi32(error01_block1, mul01);
                        error02_block1 = _mm256_add_epi32(error02_block1, mul02);
                        error03_block1 = _mm256_add_epi32(error03_block1, mul03);
                        error04_block1 = _mm256_add_epi32(error04_block1, mul04);
                        error05_block1 = _mm256_add_epi32(error05_block1, mul05);
                        error06_block1 = _mm256_add_epi32(error06_block1, mul06);
                        error07_block1 = _mm256_add_epi32(error07_block1, mul07);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block1 mulsum part3
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b11110000));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b11110000));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b11110000));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b11110000));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b11110000));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b11110000));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b11110000));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b11110000));

                        // block2 l2norm part1
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);

                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_blend_epi32(sqrdiff00, zero, 0b00001111));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_blend_epi32(sqrdiff01, zero, 0b00001111));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_blend_epi32(sqrdiff02, zero, 0b00001111));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_blend_epi32(sqrdiff03, zero, 0b00001111));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_blend_epi32(sqrdiff04, zero, 0b00001111));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_blend_epi32(sqrdiff05, zero, 0b00001111));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_blend_epi32(sqrdiff06, zero, 0b00001111));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_blend_epi32(sqrdiff07, zero, 0b00001111));

                        m += 16;
                    }
                    for (; m < blocksize * 3 - 15; m += 16)
                    {
                        // block2 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_madd_epi16(diff00, diff00));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_madd_epi16(diff01, diff01));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_madd_epi16(diff02, diff02));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_madd_epi16(diff03, diff03));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_madd_epi16(diff04, diff04));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_madd_epi16(diff05, diff05));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_madd_epi16(diff06, diff06));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_madd_epi16(diff07, diff07));
                    }
                }

                error00_block02 = error00_block2;
                error01_block02 = error01_block2;
                error02_block02 = error02_block2;
                error03_block02 = error03_block2;
                error04_block02 = error04_block2;
                error05_block02 = error05_block2;
                error06_block02 = error06_block2;
                error07_block02 = error07_block2;

                error00_block13 = _mm256_slli_epi32(error00_block1, 1);
                error01_block13 = _mm256_slli_epi32(error01_block1, 1);
                error02_block13 = _mm256_slli_epi32(error02_block1, 1);
                error03_block13 = _mm256_slli_epi32(error03_block1, 1);
                error04_block13 = _mm256_slli_epi32(error04_block1, 1);
                error05_block13 = _mm256_slli_epi32(error05_block1, 1);
                error06_block13 = _mm256_slli_epi32(error06_block1, 1);
                error07_block13 = _mm256_slli_epi32(error07_block1, 1);

                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                error01 = _mm256_sub_epi32(error01_block02, error01_block13);
                error02 = _mm256_sub_epi32(error02_block02, error02_block13);
                error03 = _mm256_sub_epi32(error03_block02, error03_block13);
                error04 = _mm256_sub_epi32(error04_block02, error04_block13);
                error05 = _mm256_sub_epi32(error05_block02, error05_block13);
                error06 = _mm256_sub_epi32(error06_block02, error06_block13);
                error07 = _mm256_sub_epi32(error07_block02, error07_block13);

                //__m256i error = _mm256_set_epi32(hsum_8x32(error07), hsum_8x32(error06), hsum_8x32(error05), hsum_8x32(error04), hsum_8x32(error03), hsum_8x32(error02), hsum_8x32(error01), hsum_8x32(error00));
                __m256i error = hsum_8x8x32(error00, error01, error02, error03, error04, error05, error06, error07);
                __m256i error_with_integral = _mm256_add_epi32(error, block_integral);
                _mm256_storeu_si256((__m256i_u *)&errors[irow * error_width + icol], error_with_integral);
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
        error_t block1_out_integral_scalar = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        __m256i block1_out_integral = _mm256_set1_epi32(block1_out_integral_scalar);
        // precalculation for block 3
        int block3_out_integral_start = (lrow + overlap) * integral_width + lcol + (blocksize - overlap);
        error_t block3_out_integral_scalar = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);
        __m256i block3_out_integral = _mm256_set1_epi32(block3_out_integral_scalar);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width - 7; icol += 8)
            {
                __m256i error00, error00_block0, error00_block1, error00_block2, error00_block3, error00_block02, error00_block13;
                __m256i error01, error01_block0, error01_block1, error01_block2, error01_block3, error01_block02, error01_block13;
                __m256i error02, error02_block0, error02_block1, error02_block2, error02_block3, error02_block02, error02_block13;
                __m256i error03, error03_block0, error03_block1, error03_block2, error03_block3, error03_block02, error03_block13;
                __m256i error04, error04_block0, error04_block1, error04_block2, error04_block3, error04_block02, error04_block13;
                __m256i error05, error05_block0, error05_block1, error05_block2, error05_block3, error05_block02, error05_block13;
                __m256i error06, error06_block0, error06_block1, error06_block2, error06_block3, error06_block02, error06_block13;
                __m256i error07, error07_block0, error07_block1, error07_block2, error07_block3, error07_block02, error07_block13;

                error00 = error00_block0 = error00_block1 = error00_block2 = error00_block3 = error00_block02 = error00_block13 = zero;
                error01 = error01_block0 = error01_block1 = error01_block2 = error01_block3 = error01_block02 = error01_block13 = zero;
                error02 = error02_block0 = error02_block1 = error02_block2 = error02_block3 = error02_block02 = error02_block13 = zero;
                error03 = error03_block0 = error03_block1 = error03_block2 = error03_block3 = error03_block02 = error03_block13 = zero;
                error04 = error04_block0 = error04_block1 = error04_block2 = error04_block3 = error04_block02 = error04_block13 = zero;
                error05 = error05_block0 = error05_block1 = error05_block2 = error05_block3 = error05_block02 = error05_block13 = zero;
                error06 = error06_block0 = error06_block1 = error06_block2 = error06_block3 = error06_block02 = error06_block13 = zero;
                error07 = error07_block0 = error07_block1 = error07_block2 = error07_block3 = error07_block02 = error07_block13 = zero;

                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                __m256i block1_in_integral = integralx8(block1_in_integral_start, height1, width1, integral_width, integral);

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                __m256i block3_in_integral = integralx8(block3_in_integral_start, height3, width3, integral_width, integral);

                __m256i block1_integral = _mm256_add_epi32(block1_in_integral, block1_out_integral);
                __m256i block3_integral = _mm256_add_epi32(block3_in_integral, block3_out_integral);
                __m256i block_integral = _mm256_add_epi32(block1_integral, block3_integral);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block0 l2norm part1
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_madd_epi16(diff00, diff00));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_madd_epi16(diff01, diff01));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_madd_epi16(diff02, diff02));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_madd_epi16(diff03, diff03));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_madd_epi16(diff04, diff04));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_madd_epi16(diff05, diff05));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_madd_epi16(diff06, diff06));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_madd_epi16(diff07, diff07));
                    }
                    if (overlap % 16 != 0)
                    {
                        // block0 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);

                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_blend_epi32(sqrdiff00, zero, 0b11110000));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_blend_epi32(sqrdiff01, zero, 0b11110000));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_blend_epi32(sqrdiff02, zero, 0b11110000));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_blend_epi32(sqrdiff03, zero, 0b11110000));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_blend_epi32(sqrdiff04, zero, 0b11110000));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_blend_epi32(sqrdiff05, zero, 0b11110000));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_blend_epi32(sqrdiff06, zero, 0b11110000));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_blend_epi32(sqrdiff07, zero, 0b11110000));

                        // block1 mulsum part1
                        // above_in = above-block input number 00
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b00001111));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b00001111));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b00001111));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b00001111));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b00001111));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b00001111));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b00001111));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b00001111));

                        m += 16;
                    }
                    for (; m < (blocksize - overlap) * 3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                        error01_block1 = _mm256_add_epi32(error01_block1, mul01);
                        error02_block1 = _mm256_add_epi32(error02_block1, mul02);
                        error03_block1 = _mm256_add_epi32(error03_block1, mul03);
                        error04_block1 = _mm256_add_epi32(error04_block1, mul04);
                        error05_block1 = _mm256_add_epi32(error05_block1, mul05);
                        error06_block1 = _mm256_add_epi32(error06_block1, mul06);
                        error07_block1 = _mm256_add_epi32(error07_block1, mul07);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block1 mulsum part3
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b11110000));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b11110000));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b11110000));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b11110000));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b11110000));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b11110000));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b11110000));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b11110000));

                        // block2 l2norm part1
                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);

                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_blend_epi32(sqrdiff00, zero, 0b00001111));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_blend_epi32(sqrdiff01, zero, 0b00001111));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_blend_epi32(sqrdiff02, zero, 0b00001111));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_blend_epi32(sqrdiff03, zero, 0b00001111));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_blend_epi32(sqrdiff04, zero, 0b00001111));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_blend_epi32(sqrdiff05, zero, 0b00001111));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_blend_epi32(sqrdiff06, zero, 0b00001111));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_blend_epi32(sqrdiff07, zero, 0b00001111));

                        m += 16;
                    }
                    for (; m < blocksize * 3 - 15; m += 16)
                    {
                        // block2 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_madd_epi16(diff00, diff00));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_madd_epi16(diff01, diff01));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_madd_epi16(diff02, diff02));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_madd_epi16(diff03, diff03));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_madd_epi16(diff04, diff04));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_madd_epi16(diff05, diff05));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_madd_epi16(diff06, diff06));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_madd_epi16(diff07, diff07));
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
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);

                        error00_block3 = _mm256_add_epi32(mull00, error00_block3);
                        error01_block3 = _mm256_add_epi32(mull01, error01_block3);
                        error02_block3 = _mm256_add_epi32(mull02, error02_block3);
                        error03_block3 = _mm256_add_epi32(mull03, error03_block3);
                        error04_block3 = _mm256_add_epi32(mull04, error04_block3);
                        error05_block3 = _mm256_add_epi32(mull05, error05_block3);
                        error06_block3 = _mm256_add_epi32(mull06, error06_block3);
                        error07_block3 = _mm256_add_epi32(mull07, error07_block3);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block3 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);

                        error00_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00_block3);
                        error01_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull01, zero, 0b11110000), error01_block3);
                        error02_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull02, zero, 0b11110000), error02_block3);
                        error03_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull03, zero, 0b11110000), error03_block3);
                        error04_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull04, zero, 0b11110000), error04_block3);
                        error05_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull05, zero, 0b11110000), error05_block3);
                        error06_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull06, zero, 0b11110000), error06_block3);
                        error07_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull07, zero, 0b11110000), error07_block3);
                    }
                }

                error00_block02 = _mm256_add_epi32(error00_block0, error00_block2);
                error01_block02 = _mm256_add_epi32(error01_block0, error01_block2);
                error02_block02 = _mm256_add_epi32(error02_block0, error02_block2);
                error03_block02 = _mm256_add_epi32(error03_block0, error03_block2);
                error04_block02 = _mm256_add_epi32(error04_block0, error04_block2);
                error05_block02 = _mm256_add_epi32(error05_block0, error05_block2);
                error06_block02 = _mm256_add_epi32(error06_block0, error06_block2);
                error07_block02 = _mm256_add_epi32(error07_block0, error07_block2);

                error00_block13 = _mm256_slli_epi32(_mm256_add_epi32(error00_block1, error00_block3), 1);
                error01_block13 = _mm256_slli_epi32(_mm256_add_epi32(error01_block1, error01_block3), 1);
                error02_block13 = _mm256_slli_epi32(_mm256_add_epi32(error02_block1, error02_block3), 1);
                error03_block13 = _mm256_slli_epi32(_mm256_add_epi32(error03_block1, error03_block3), 1);
                error04_block13 = _mm256_slli_epi32(_mm256_add_epi32(error04_block1, error04_block3), 1);
                error05_block13 = _mm256_slli_epi32(_mm256_add_epi32(error05_block1, error05_block3), 1);
                error06_block13 = _mm256_slli_epi32(_mm256_add_epi32(error06_block1, error06_block3), 1);
                error07_block13 = _mm256_slli_epi32(_mm256_add_epi32(error07_block1, error07_block3), 1);

                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                error01 = _mm256_sub_epi32(error01_block02, error01_block13);
                error02 = _mm256_sub_epi32(error02_block02, error02_block13);
                error03 = _mm256_sub_epi32(error03_block02, error03_block13);
                error04 = _mm256_sub_epi32(error04_block02, error04_block13);
                error05 = _mm256_sub_epi32(error05_block02, error05_block13);
                error06 = _mm256_sub_epi32(error06_block02, error06_block13);
                error07 = _mm256_sub_epi32(error07_block02, error07_block13);

                //__m256i error = _mm256_set_epi32(hsum_8x32(error07), hsum_8x32(error06), hsum_8x32(error05), hsum_8x32(error04), hsum_8x32(error03), hsum_8x32(error02), hsum_8x32(error01), hsum_8x32(error00));
                __m256i error = hsum_8x8x32(error00, error01, error02, error03, error04, error05, error06, error07);
                __m256i error_with_integral = _mm256_add_epi32(error, block_integral);
                _mm256_storeu_si256((__m256i_u *)&errors[irow * error_width + icol], error_with_integral);
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
        error_t block1_out_integral_scalar = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        __m256i block1_out_integral = _mm256_set1_epi32(block1_out_integral_scalar);
        // precalculation for block 3
        int block3_out_integral_start = (lrow + overlap) * integral_width + lcol + (blocksize - overlap);
        error_t block3_out_integral_scalar = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);
        __m256i block3_out_integral = _mm256_set1_epi32(block3_out_integral_scalar);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width - 7; icol += 8)
            {
                __m256i error00, error00_block0, error00_block1, error00_block3, error00_block02, error00_block13;
                __m256i error01, error01_block0, error01_block1, error01_block3, error01_block02, error01_block13;
                __m256i error02, error02_block0, error02_block1, error02_block3, error02_block02, error02_block13;
                __m256i error03, error03_block0, error03_block1, error03_block3, error03_block02, error03_block13;
                __m256i error04, error04_block0, error04_block1, error04_block3, error04_block02, error04_block13;
                __m256i error05, error05_block0, error05_block1, error05_block3, error05_block02, error05_block13;
                __m256i error06, error06_block0, error06_block1, error06_block3, error06_block02, error06_block13;
                __m256i error07, error07_block0, error07_block1, error07_block3, error07_block02, error07_block13;

                error00 = error00_block0 = error00_block1 = error00_block3 = error00_block02 = error00_block13 = zero;
                error01 = error01_block0 = error01_block1 = error01_block3 = error01_block02 = error01_block13 = zero;
                error02 = error02_block0 = error02_block1 = error02_block3 = error02_block02 = error02_block13 = zero;
                error03 = error03_block0 = error03_block1 = error03_block3 = error03_block02 = error03_block13 = zero;
                error04 = error04_block0 = error04_block1 = error04_block3 = error04_block02 = error04_block13 = zero;
                error05 = error05_block0 = error05_block1 = error05_block3 = error05_block02 = error05_block13 = zero;
                error06 = error06_block0 = error06_block1 = error06_block3 = error06_block02 = error06_block13 = zero;
                error07 = error07_block0 = error07_block1 = error07_block3 = error07_block02 = error07_block13 = zero;

                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start = block1_in_integral_base + irow * integral_width + icol;
                __m256i block1_in_integral = integralx8(block1_in_integral_start, height1, width1, integral_width, integral);

                int block3_in_integral_start = block3_in_integral_base + irow * integral_width + icol;
                __m256i block3_in_integral = integralx8(block3_in_integral_start, height3, width3, integral_width, integral);

                __m256i block1_integral = _mm256_add_epi32(block1_in_integral, block1_out_integral);
                __m256i block3_integral = _mm256_add_epi32(block3_in_integral, block3_out_integral);
                __m256i block_integral = _mm256_add_epi32(block1_integral, block3_integral);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block0 l2norm part1
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_madd_epi16(diff00, diff00));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_madd_epi16(diff01, diff01));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_madd_epi16(diff02, diff02));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_madd_epi16(diff03, diff03));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_madd_epi16(diff04, diff04));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_madd_epi16(diff05, diff05));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_madd_epi16(diff06, diff06));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_madd_epi16(diff07, diff07));
                    }
                    if (overlap % 16 != 0)
                    {
                        // block0 l2norm part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i curr_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);

                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);

                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_blend_epi32(sqrdiff00, zero, 0b11110000));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_blend_epi32(sqrdiff01, zero, 0b11110000));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_blend_epi32(sqrdiff02, zero, 0b11110000));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_blend_epi32(sqrdiff03, zero, 0b11110000));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_blend_epi32(sqrdiff04, zero, 0b11110000));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_blend_epi32(sqrdiff05, zero, 0b11110000));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_blend_epi32(sqrdiff06, zero, 0b11110000));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_blend_epi32(sqrdiff07, zero, 0b11110000));

                        // block1 mulsum part1
                        // above_in = above-block input number 00
                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b00001111));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b00001111));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b00001111));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b00001111));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b00001111));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b00001111));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b00001111));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b00001111));

                        m += 16;
                    }
                    for (; m < (blocksize)*3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i above_in = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);

                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                        error01_block1 = _mm256_add_epi32(error01_block1, mul01);
                        error02_block1 = _mm256_add_epi32(error02_block1, mul02);
                        error03_block1 = _mm256_add_epi32(error03_block1, mul03);
                        error04_block1 = _mm256_add_epi32(error04_block1, mul04);
                        error05_block1 = _mm256_add_epi32(error05_block1, mul05);
                        error06_block1 = _mm256_add_epi32(error06_block1, mul06);
                        error07_block1 = _mm256_add_epi32(error07_block1, mul07);
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
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);

                        error00_block3 = _mm256_add_epi32(mull00, error00_block3);
                        error01_block3 = _mm256_add_epi32(mull01, error01_block3);
                        error02_block3 = _mm256_add_epi32(mull02, error02_block3);
                        error03_block3 = _mm256_add_epi32(mull03, error03_block3);
                        error04_block3 = _mm256_add_epi32(mull04, error04_block3);
                        error05_block3 = _mm256_add_epi32(mull05, error05_block3);
                        error06_block3 = _mm256_add_epi32(mull06, error06_block3);
                        error07_block3 = _mm256_add_epi32(mull07, error07_block3);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block3 mulsum part2
                        __m256i curr_in_00 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);

                        __m256i left_in = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);

                        error00_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00_block3);
                        error01_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull01, zero, 0b11110000), error01_block3);
                        error02_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull02, zero, 0b11110000), error02_block3);
                        error03_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull03, zero, 0b11110000), error03_block3);
                        error04_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull04, zero, 0b11110000), error04_block3);
                        error05_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull05, zero, 0b11110000), error05_block3);
                        error06_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull06, zero, 0b11110000), error06_block3);
                        error07_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull07, zero, 0b11110000), error07_block3);
                    }
                }

                error00_block02 = error00_block0;
                error01_block02 = error01_block0;
                error02_block02 = error02_block0;
                error03_block02 = error03_block0;
                error04_block02 = error04_block0;
                error05_block02 = error05_block0;
                error06_block02 = error06_block0;
                error07_block02 = error07_block0;

                error00_block13 = _mm256_slli_epi32(_mm256_add_epi32(error00_block1, error00_block3), 1);
                error01_block13 = _mm256_slli_epi32(_mm256_add_epi32(error01_block1, error01_block3), 1);
                error02_block13 = _mm256_slli_epi32(_mm256_add_epi32(error02_block1, error02_block3), 1);
                error03_block13 = _mm256_slli_epi32(_mm256_add_epi32(error03_block1, error03_block3), 1);
                error04_block13 = _mm256_slli_epi32(_mm256_add_epi32(error04_block1, error04_block3), 1);
                error05_block13 = _mm256_slli_epi32(_mm256_add_epi32(error05_block1, error05_block3), 1);
                error06_block13 = _mm256_slli_epi32(_mm256_add_epi32(error06_block1, error06_block3), 1);
                error07_block13 = _mm256_slli_epi32(_mm256_add_epi32(error07_block1, error07_block3), 1);

                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                error01 = _mm256_sub_epi32(error01_block02, error01_block13);
                error02 = _mm256_sub_epi32(error02_block02, error02_block13);
                error03 = _mm256_sub_epi32(error03_block02, error03_block13);
                error04 = _mm256_sub_epi32(error04_block02, error04_block13);
                error05 = _mm256_sub_epi32(error05_block02, error05_block13);
                error06 = _mm256_sub_epi32(error06_block02, error06_block13);
                error07 = _mm256_sub_epi32(error07_block02, error07_block13);
                //__m256i error = _mm256_set_epi32(hsum_8x32(error07), hsum_8x32(error06), hsum_8x32(error05), hsum_8x32(error04), hsum_8x32(error03), hsum_8x32(error02), hsum_8x32(error01), hsum_8x32(error00));
                __m256i error = hsum_8x8x32(error00, error01, error02, error03, error04, error05, error06, error07);
                __m256i error_with_integral = _mm256_add_epi32(error, block_integral);
                _mm256_storeu_si256((__m256i_u *)&errors[irow * error_width + icol], error_with_integral);
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
