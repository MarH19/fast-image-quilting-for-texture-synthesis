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

static inline __m256i loadu_8x16(pixel_t* v)
{
    return _mm256_cvtepu8_epi16 (_mm_loadu_si128((__m128i_u *) v));
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
            for (int icol = 0; icol < error_width - 15; icol += 16)
            {
                __m256i error00 = zero;
                __m256i error01 = zero;
                __m256i error02 = zero;
                __m256i error03 = zero;
                __m256i error04 = zero;
                __m256i error05 = zero;
                __m256i error06 = zero;
                __m256i error07 = zero;
                __m256i error08 = zero;
                __m256i error09 = zero;
                __m256i error10 = zero;
                __m256i error11 = zero;
                __m256i error12 = zero;
                __m256i error13 = zero;
                __m256i error14 = zero;
                __m256i error15 = zero;

                int block3_in_integral_start00 = block3_in_integral_base + irow * integral_width + icol + 0;
                int block3_in_integral_start01 = block3_in_integral_base + irow * integral_width + icol + 1;
                int block3_in_integral_start02 = block3_in_integral_base + irow * integral_width + icol + 2;
                int block3_in_integral_start03 = block3_in_integral_base + irow * integral_width + icol + 3;
                int block3_in_integral_start04 = block3_in_integral_base + irow * integral_width + icol + 4;
                int block3_in_integral_start05 = block3_in_integral_base + irow * integral_width + icol + 5;
                int block3_in_integral_start06 = block3_in_integral_base + irow * integral_width + icol + 6;
                int block3_in_integral_start07 = block3_in_integral_base + irow * integral_width + icol + 7;
                int block3_in_integral_start08 = block3_in_integral_base + irow * integral_width + icol + 8;
                int block3_in_integral_start09 = block3_in_integral_base + irow * integral_width + icol + 9;
                int block3_in_integral_start10 = block3_in_integral_base + irow * integral_width + icol + 10;
                int block3_in_integral_start11 = block3_in_integral_base + irow * integral_width + icol + 11;
                int block3_in_integral_start12 = block3_in_integral_base + irow * integral_width + icol + 12;
                int block3_in_integral_start13 = block3_in_integral_base + irow * integral_width + icol + 13;
                int block3_in_integral_start14 = block3_in_integral_base + irow * integral_width + icol + 14;
                int block3_in_integral_start15 = block3_in_integral_base + irow * integral_width + icol + 15;
                error_t block3_in_integral00 = INTEGRAL(block3_in_integral_start00, height3, width3, integral_width);
                error_t block3_in_integral01 = INTEGRAL(block3_in_integral_start01, height3, width3, integral_width);
                error_t block3_in_integral02 = INTEGRAL(block3_in_integral_start02, height3, width3, integral_width);
                error_t block3_in_integral03 = INTEGRAL(block3_in_integral_start03, height3, width3, integral_width);
                error_t block3_in_integral04 = INTEGRAL(block3_in_integral_start04, height3, width3, integral_width);
                error_t block3_in_integral05 = INTEGRAL(block3_in_integral_start05, height3, width3, integral_width);
                error_t block3_in_integral06 = INTEGRAL(block3_in_integral_start06, height3, width3, integral_width);
                error_t block3_in_integral07 = INTEGRAL(block3_in_integral_start07, height3, width3, integral_width);
                error_t block3_in_integral08 = INTEGRAL(block3_in_integral_start08, height3, width3, integral_width);
                error_t block3_in_integral09 = INTEGRAL(block3_in_integral_start09, height3, width3, integral_width);
                error_t block3_in_integral10 = INTEGRAL(block3_in_integral_start10, height3, width3, integral_width);
                error_t block3_in_integral11 = INTEGRAL(block3_in_integral_start11, height3, width3, integral_width);
                error_t block3_in_integral12 = INTEGRAL(block3_in_integral_start12, height3, width3, integral_width);
                error_t block3_in_integral13 = INTEGRAL(block3_in_integral_start13, height3, width3, integral_width);
                error_t block3_in_integral14 = INTEGRAL(block3_in_integral_start14, height3, width3, integral_width);
                error_t block3_in_integral15 = INTEGRAL(block3_in_integral_start15, height3, width3, integral_width);

                for (int k = 0; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i left_in = loadu_8x16(&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);
                        __m256i mull08 = _mm256_madd_epi16(curr_in_08, left_in);
                        __m256i mull09 = _mm256_madd_epi16(curr_in_09, left_in);
                        __m256i mull10 = _mm256_madd_epi16(curr_in_10, left_in);
                        __m256i mull11 = _mm256_madd_epi16(curr_in_11, left_in);
                        __m256i mull12 = _mm256_madd_epi16(curr_in_12, left_in);
                        __m256i mull13 = _mm256_madd_epi16(curr_in_13, left_in);
                        __m256i mull14 = _mm256_madd_epi16(curr_in_14, left_in);
                        __m256i mull15 = _mm256_madd_epi16(curr_in_15, left_in);
                        error00 = _mm256_add_epi32(mull00, error00);
                        error01 = _mm256_add_epi32(mull01, error01);
                        error02 = _mm256_add_epi32(mull02, error02);
                        error03 = _mm256_add_epi32(mull03, error03);
                        error04 = _mm256_add_epi32(mull04, error04);
                        error05 = _mm256_add_epi32(mull05, error05);
                        error06 = _mm256_add_epi32(mull06, error06);
                        error07 = _mm256_add_epi32(mull07, error07);
                        error08 = _mm256_add_epi32(mull08, error08);
                        error09 = _mm256_add_epi32(mull09, error09);
                        error10 = _mm256_add_epi32(mull10, error10);
                        error11 = _mm256_add_epi32(mull11, error11);
                        error12 = _mm256_add_epi32(mull12, error12);
                        error13 = _mm256_add_epi32(mull13, error13);
                        error14 = _mm256_add_epi32(mull14, error14);
                        error15 = _mm256_add_epi32(mull15, error15);
                    }
                    if (overlap % 16 != 0)
                    {
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i left_in = loadu_8x16(&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);
                        __m256i mull08 = _mm256_madd_epi16(curr_in_08, left_in);
                        __m256i mull09 = _mm256_madd_epi16(curr_in_09, left_in);
                        __m256i mull10 = _mm256_madd_epi16(curr_in_10, left_in);
                        __m256i mull11 = _mm256_madd_epi16(curr_in_11, left_in);
                        __m256i mull12 = _mm256_madd_epi16(curr_in_12, left_in);
                        __m256i mull13 = _mm256_madd_epi16(curr_in_13, left_in);
                        __m256i mull14 = _mm256_madd_epi16(curr_in_14, left_in);
                        __m256i mull15 = _mm256_madd_epi16(curr_in_15, left_in);
                        error00 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00);
                        error01 = _mm256_add_epi32(_mm256_blend_epi32(mull01, zero, 0b11110000), error01);
                        error02 = _mm256_add_epi32(_mm256_blend_epi32(mull02, zero, 0b11110000), error02);
                        error03 = _mm256_add_epi32(_mm256_blend_epi32(mull03, zero, 0b11110000), error03);
                        error04 = _mm256_add_epi32(_mm256_blend_epi32(mull04, zero, 0b11110000), error04);
                        error05 = _mm256_add_epi32(_mm256_blend_epi32(mull05, zero, 0b11110000), error05);
                        error06 = _mm256_add_epi32(_mm256_blend_epi32(mull06, zero, 0b11110000), error06);
                        error07 = _mm256_add_epi32(_mm256_blend_epi32(mull07, zero, 0b11110000), error07);
                        error08 = _mm256_add_epi32(_mm256_blend_epi32(mull08, zero, 0b11110000), error08);
                        error09 = _mm256_add_epi32(_mm256_blend_epi32(mull09, zero, 0b11110000), error09);
                        error10 = _mm256_add_epi32(_mm256_blend_epi32(mull10, zero, 0b11110000), error10);
                        error11 = _mm256_add_epi32(_mm256_blend_epi32(mull11, zero, 0b11110000), error11);
                        error12 = _mm256_add_epi32(_mm256_blend_epi32(mull12, zero, 0b11110000), error12);
                        error13 = _mm256_add_epi32(_mm256_blend_epi32(mull13, zero, 0b11110000), error13);
                        error14 = _mm256_add_epi32(_mm256_blend_epi32(mull14, zero, 0b11110000), error14);
                        error15 = _mm256_add_epi32(_mm256_blend_epi32(mull15, zero, 0b11110000), error15);
                    }
                }
                errors[irow * error_width + icol + 0] = block3_out_integral - 2 * hsum_8x32(error00) + block3_in_integral00;
                errors[irow * error_width + icol + 1] = block3_out_integral - 2 * hsum_8x32(error01) + block3_in_integral01;
                errors[irow * error_width + icol + 2] = block3_out_integral - 2 * hsum_8x32(error02) + block3_in_integral02;
                errors[irow * error_width + icol + 3] = block3_out_integral - 2 * hsum_8x32(error03) + block3_in_integral03;
                errors[irow * error_width + icol + 4] = block3_out_integral - 2 * hsum_8x32(error04) + block3_in_integral04;
                errors[irow * error_width + icol + 5] = block3_out_integral - 2 * hsum_8x32(error05) + block3_in_integral05;
                errors[irow * error_width + icol + 6] = block3_out_integral - 2 * hsum_8x32(error06) + block3_in_integral06;
                errors[irow * error_width + icol + 7] = block3_out_integral - 2 * hsum_8x32(error07) + block3_in_integral07;
                errors[irow * error_width + icol + 8] = block3_out_integral - 2 * hsum_8x32(error08) + block3_in_integral08;
                errors[irow * error_width + icol + 9] = block3_out_integral - 2 * hsum_8x32(error09) + block3_in_integral09;
                errors[irow * error_width + icol + 10] = block3_out_integral - 2 * hsum_8x32(error10) + block3_in_integral10;
                errors[irow * error_width + icol + 11] = block3_out_integral - 2 * hsum_8x32(error11) + block3_in_integral11;
                errors[irow * error_width + icol + 12] = block3_out_integral - 2 * hsum_8x32(error12) + block3_in_integral12;
                errors[irow * error_width + icol + 13] = block3_out_integral - 2 * hsum_8x32(error13) + block3_in_integral13;
                errors[irow * error_width + icol + 14] = block3_out_integral - 2 * hsum_8x32(error14) + block3_in_integral14;
                errors[irow * error_width + icol + 15] = block3_out_integral - 2 * hsum_8x32(error15) + block3_in_integral15;
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
            for (int icol = 0; icol < error_width - 15; icol += 16)
            {
                __m256i error00, error00_block1, error00_block2, error00_block02, error00_block13;
                __m256i error01, error01_block1, error01_block2, error01_block02, error01_block13;
                __m256i error02, error02_block1, error02_block2, error02_block02, error02_block13;
                __m256i error03, error03_block1, error03_block2, error03_block02, error03_block13;
                __m256i error04, error04_block1, error04_block2, error04_block02, error04_block13;
                __m256i error05, error05_block1, error05_block2, error05_block02, error05_block13;
                __m256i error06, error06_block1, error06_block2, error06_block02, error06_block13;
                __m256i error07, error07_block1, error07_block2, error07_block02, error07_block13;
                __m256i error08, error08_block1, error08_block2, error08_block02, error08_block13;
                __m256i error09, error09_block1, error09_block2, error09_block02, error09_block13;
                __m256i error10, error10_block1, error10_block2, error10_block02, error10_block13;
                __m256i error11, error11_block1, error11_block2, error11_block02, error11_block13;
                __m256i error12, error12_block1, error12_block2, error12_block02, error12_block13;
                __m256i error13, error13_block1, error13_block2, error13_block02, error13_block13;
                __m256i error14, error14_block1, error14_block2, error14_block02, error14_block13;
                __m256i error15, error15_block1, error15_block2, error15_block02, error15_block13;
                error00 = error00_block1 = error00_block2 = error00_block02 = error00_block13 = zero;
                error01 = error01_block1 = error01_block2 = error01_block02 = error01_block13 = zero;
                error02 = error02_block1 = error02_block2 = error02_block02 = error02_block13 = zero;
                error03 = error03_block1 = error03_block2 = error03_block02 = error03_block13 = zero;
                error04 = error04_block1 = error04_block2 = error04_block02 = error04_block13 = zero;
                error05 = error05_block1 = error05_block2 = error05_block02 = error05_block13 = zero;
                error06 = error06_block1 = error06_block2 = error06_block02 = error06_block13 = zero;
                error07 = error07_block1 = error07_block2 = error07_block02 = error07_block13 = zero;
                error08 = error08_block1 = error08_block2 = error08_block02 = error08_block13 = zero;
                error09 = error09_block1 = error09_block2 = error09_block02 = error09_block13 = zero;
                error10 = error10_block1 = error10_block2 = error10_block02 = error10_block13 = zero;
                error11 = error11_block1 = error11_block2 = error11_block02 = error11_block13 = zero;
                error12 = error12_block1 = error12_block2 = error12_block02 = error12_block13 = zero;
                error13 = error13_block1 = error13_block2 = error13_block02 = error13_block13 = zero;
                error14 = error14_block1 = error14_block2 = error14_block02 = error14_block13 = zero;
                error15 = error15_block1 = error15_block2 = error15_block02 = error15_block13 = zero;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start00 = block1_in_integral_base + irow * integral_width + icol + 0;
                int block1_in_integral_start01 = block1_in_integral_base + irow * integral_width + icol + 1;
                int block1_in_integral_start02 = block1_in_integral_base + irow * integral_width + icol + 2;
                int block1_in_integral_start03 = block1_in_integral_base + irow * integral_width + icol + 3;
                int block1_in_integral_start04 = block1_in_integral_base + irow * integral_width + icol + 4;
                int block1_in_integral_start05 = block1_in_integral_base + irow * integral_width + icol + 5;
                int block1_in_integral_start06 = block1_in_integral_base + irow * integral_width + icol + 6;
                int block1_in_integral_start07 = block1_in_integral_base + irow * integral_width + icol + 7;
                int block1_in_integral_start08 = block1_in_integral_base + irow * integral_width + icol + 8;
                int block1_in_integral_start09 = block1_in_integral_base + irow * integral_width + icol + 9;
                int block1_in_integral_start10 = block1_in_integral_base + irow * integral_width + icol + 10;
                int block1_in_integral_start11 = block1_in_integral_base + irow * integral_width + icol + 11;
                int block1_in_integral_start12 = block1_in_integral_base + irow * integral_width + icol + 12;
                int block1_in_integral_start13 = block1_in_integral_base + irow * integral_width + icol + 13;
                int block1_in_integral_start14 = block1_in_integral_base + irow * integral_width + icol + 14;
                int block1_in_integral_start15 = block1_in_integral_base + irow * integral_width + icol + 15;
                error_t block1_in_integral00 = INTEGRAL(block1_in_integral_start00, height1, width1, integral_width);
                error_t block1_in_integral01 = INTEGRAL(block1_in_integral_start01, height1, width1, integral_width);
                error_t block1_in_integral02 = INTEGRAL(block1_in_integral_start02, height1, width1, integral_width);
                error_t block1_in_integral03 = INTEGRAL(block1_in_integral_start03, height1, width1, integral_width);
                error_t block1_in_integral04 = INTEGRAL(block1_in_integral_start04, height1, width1, integral_width);
                error_t block1_in_integral05 = INTEGRAL(block1_in_integral_start05, height1, width1, integral_width);
                error_t block1_in_integral06 = INTEGRAL(block1_in_integral_start06, height1, width1, integral_width);
                error_t block1_in_integral07 = INTEGRAL(block1_in_integral_start07, height1, width1, integral_width);
                error_t block1_in_integral08 = INTEGRAL(block1_in_integral_start08, height1, width1, integral_width);
                error_t block1_in_integral09 = INTEGRAL(block1_in_integral_start09, height1, width1, integral_width);
                error_t block1_in_integral10 = INTEGRAL(block1_in_integral_start10, height1, width1, integral_width);
                error_t block1_in_integral11 = INTEGRAL(block1_in_integral_start11, height1, width1, integral_width);
                error_t block1_in_integral12 = INTEGRAL(block1_in_integral_start12, height1, width1, integral_width);
                error_t block1_in_integral13 = INTEGRAL(block1_in_integral_start13, height1, width1, integral_width);
                error_t block1_in_integral14 = INTEGRAL(block1_in_integral_start14, height1, width1, integral_width);
                error_t block1_in_integral15 = INTEGRAL(block1_in_integral_start15, height1, width1, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < (blocksize - overlap) * 3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                        error01_block1 = _mm256_add_epi32(error01_block1, mul01);
                        error02_block1 = _mm256_add_epi32(error02_block1, mul02);
                        error03_block1 = _mm256_add_epi32(error03_block1, mul03);
                        error04_block1 = _mm256_add_epi32(error04_block1, mul04);
                        error05_block1 = _mm256_add_epi32(error05_block1, mul05);
                        error06_block1 = _mm256_add_epi32(error06_block1, mul06);
                        error07_block1 = _mm256_add_epi32(error07_block1, mul07);
                        error08_block1 = _mm256_add_epi32(error08_block1, mul08);
                        error09_block1 = _mm256_add_epi32(error09_block1, mul09);
                        error10_block1 = _mm256_add_epi32(error10_block1, mul10);
                        error11_block1 = _mm256_add_epi32(error11_block1, mul11);
                        error12_block1 = _mm256_add_epi32(error12_block1, mul12);
                        error13_block1 = _mm256_add_epi32(error13_block1, mul13);
                        error14_block1 = _mm256_add_epi32(error14_block1, mul14);
                        error15_block1 = _mm256_add_epi32(error15_block1, mul15);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block1 mulsum part3
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b11110000));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b11110000));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b11110000));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b11110000));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b11110000));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b11110000));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b11110000));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b11110000));
                        error08_block1 = _mm256_add_epi32(error08_block1, _mm256_blend_epi32(mul08, zero, 0b11110000));
                        error09_block1 = _mm256_add_epi32(error09_block1, _mm256_blend_epi32(mul09, zero, 0b11110000));
                        error10_block1 = _mm256_add_epi32(error10_block1, _mm256_blend_epi32(mul10, zero, 0b11110000));
                        error11_block1 = _mm256_add_epi32(error11_block1, _mm256_blend_epi32(mul11, zero, 0b11110000));
                        error12_block1 = _mm256_add_epi32(error12_block1, _mm256_blend_epi32(mul12, zero, 0b11110000));
                        error13_block1 = _mm256_add_epi32(error13_block1, _mm256_blend_epi32(mul13, zero, 0b11110000));
                        error14_block1 = _mm256_add_epi32(error14_block1, _mm256_blend_epi32(mul14, zero, 0b11110000));
                        error15_block1 = _mm256_add_epi32(error15_block1, _mm256_blend_epi32(mul15, zero, 0b11110000));

                        // block2 l2norm part1
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);
                        __m256i sqrdiff08 = _mm256_madd_epi16(diff08, diff08);
                        __m256i sqrdiff09 = _mm256_madd_epi16(diff09, diff09);
                        __m256i sqrdiff10 = _mm256_madd_epi16(diff10, diff10);
                        __m256i sqrdiff11 = _mm256_madd_epi16(diff11, diff11);
                        __m256i sqrdiff12 = _mm256_madd_epi16(diff12, diff12);
                        __m256i sqrdiff13 = _mm256_madd_epi16(diff13, diff13);
                        __m256i sqrdiff14 = _mm256_madd_epi16(diff14, diff14);
                        __m256i sqrdiff15 = _mm256_madd_epi16(diff15, diff15);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_blend_epi32(sqrdiff00, zero, 0b00001111));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_blend_epi32(sqrdiff01, zero, 0b00001111));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_blend_epi32(sqrdiff02, zero, 0b00001111));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_blend_epi32(sqrdiff03, zero, 0b00001111));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_blend_epi32(sqrdiff04, zero, 0b00001111));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_blend_epi32(sqrdiff05, zero, 0b00001111));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_blend_epi32(sqrdiff06, zero, 0b00001111));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_blend_epi32(sqrdiff07, zero, 0b00001111));
                        error08_block2 = _mm256_add_epi32(error08_block2, _mm256_blend_epi32(sqrdiff08, zero, 0b00001111));
                        error09_block2 = _mm256_add_epi32(error09_block2, _mm256_blend_epi32(sqrdiff09, zero, 0b00001111));
                        error10_block2 = _mm256_add_epi32(error10_block2, _mm256_blend_epi32(sqrdiff10, zero, 0b00001111));
                        error11_block2 = _mm256_add_epi32(error11_block2, _mm256_blend_epi32(sqrdiff11, zero, 0b00001111));
                        error12_block2 = _mm256_add_epi32(error12_block2, _mm256_blend_epi32(sqrdiff12, zero, 0b00001111));
                        error13_block2 = _mm256_add_epi32(error13_block2, _mm256_blend_epi32(sqrdiff13, zero, 0b00001111));
                        error14_block2 = _mm256_add_epi32(error14_block2, _mm256_blend_epi32(sqrdiff14, zero, 0b00001111));
                        error15_block2 = _mm256_add_epi32(error15_block2, _mm256_blend_epi32(sqrdiff15, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < blocksize * 3 - 15; m += 16)
                    {
                        // block2 l2norm part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_madd_epi16(diff00, diff00));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_madd_epi16(diff01, diff01));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_madd_epi16(diff02, diff02));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_madd_epi16(diff03, diff03));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_madd_epi16(diff04, diff04));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_madd_epi16(diff05, diff05));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_madd_epi16(diff06, diff06));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_madd_epi16(diff07, diff07));
                        error08_block2 = _mm256_add_epi32(error08_block2, _mm256_madd_epi16(diff08, diff08));
                        error09_block2 = _mm256_add_epi32(error09_block2, _mm256_madd_epi16(diff09, diff09));
                        error10_block2 = _mm256_add_epi32(error10_block2, _mm256_madd_epi16(diff10, diff10));
                        error11_block2 = _mm256_add_epi32(error11_block2, _mm256_madd_epi16(diff11, diff11));
                        error12_block2 = _mm256_add_epi32(error12_block2, _mm256_madd_epi16(diff12, diff12));
                        error13_block2 = _mm256_add_epi32(error13_block2, _mm256_madd_epi16(diff13, diff13));
                        error14_block2 = _mm256_add_epi32(error14_block2, _mm256_madd_epi16(diff14, diff14));
                        error15_block2 = _mm256_add_epi32(error15_block2, _mm256_madd_epi16(diff15, diff15));
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
                error08_block02 = error08_block2;
                error09_block02 = error09_block2;
                error10_block02 = error10_block2;
                error11_block02 = error11_block2;
                error12_block02 = error12_block2;
                error13_block02 = error13_block2;
                error14_block02 = error14_block2;
                error15_block02 = error15_block2;
                error00_block13 = _mm256_slli_epi32(error00_block1, 1);
                error01_block13 = _mm256_slli_epi32(error01_block1, 1);
                error02_block13 = _mm256_slli_epi32(error02_block1, 1);
                error03_block13 = _mm256_slli_epi32(error03_block1, 1);
                error04_block13 = _mm256_slli_epi32(error04_block1, 1);
                error05_block13 = _mm256_slli_epi32(error05_block1, 1);
                error06_block13 = _mm256_slli_epi32(error06_block1, 1);
                error07_block13 = _mm256_slli_epi32(error07_block1, 1);
                error08_block13 = _mm256_slli_epi32(error08_block1, 1);
                error09_block13 = _mm256_slli_epi32(error09_block1, 1);
                error10_block13 = _mm256_slli_epi32(error10_block1, 1);
                error11_block13 = _mm256_slli_epi32(error11_block1, 1);
                error12_block13 = _mm256_slli_epi32(error12_block1, 1);
                error13_block13 = _mm256_slli_epi32(error13_block1, 1);
                error14_block13 = _mm256_slli_epi32(error14_block1, 1);
                error15_block13 = _mm256_slli_epi32(error15_block1, 1);
                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                error01 = _mm256_sub_epi32(error01_block02, error01_block13);
                error02 = _mm256_sub_epi32(error02_block02, error02_block13);
                error03 = _mm256_sub_epi32(error03_block02, error03_block13);
                error04 = _mm256_sub_epi32(error04_block02, error04_block13);
                error05 = _mm256_sub_epi32(error05_block02, error05_block13);
                error06 = _mm256_sub_epi32(error06_block02, error06_block13);
                error07 = _mm256_sub_epi32(error07_block02, error07_block13);
                error08 = _mm256_sub_epi32(error08_block02, error08_block13);
                error09 = _mm256_sub_epi32(error09_block02, error09_block13);
                error10 = _mm256_sub_epi32(error10_block02, error10_block13);
                error11 = _mm256_sub_epi32(error11_block02, error11_block13);
                error12 = _mm256_sub_epi32(error12_block02, error12_block13);
                error13 = _mm256_sub_epi32(error13_block02, error13_block13);
                error14 = _mm256_sub_epi32(error14_block02, error14_block13);
                error15 = _mm256_sub_epi32(error15_block02, error15_block13);
                errors[irow * error_width + icol + 0] = block1_out_integral + block1_in_integral00 + hsum_8x32(error00);
                errors[irow * error_width + icol + 1] = block1_out_integral + block1_in_integral01 + hsum_8x32(error01);
                errors[irow * error_width + icol + 2] = block1_out_integral + block1_in_integral02 + hsum_8x32(error02);
                errors[irow * error_width + icol + 3] = block1_out_integral + block1_in_integral03 + hsum_8x32(error03);
                errors[irow * error_width + icol + 4] = block1_out_integral + block1_in_integral04 + hsum_8x32(error04);
                errors[irow * error_width + icol + 5] = block1_out_integral + block1_in_integral05 + hsum_8x32(error05);
                errors[irow * error_width + icol + 6] = block1_out_integral + block1_in_integral06 + hsum_8x32(error06);
                errors[irow * error_width + icol + 7] = block1_out_integral + block1_in_integral07 + hsum_8x32(error07);
                errors[irow * error_width + icol + 8] = block1_out_integral + block1_in_integral08 + hsum_8x32(error08);
                errors[irow * error_width + icol + 9] = block1_out_integral + block1_in_integral09 + hsum_8x32(error09);
                errors[irow * error_width + icol + 10] = block1_out_integral + block1_in_integral10 + hsum_8x32(error10);
                errors[irow * error_width + icol + 11] = block1_out_integral + block1_in_integral11 + hsum_8x32(error11);
                errors[irow * error_width + icol + 12] = block1_out_integral + block1_in_integral12 + hsum_8x32(error12);
                errors[irow * error_width + icol + 13] = block1_out_integral + block1_in_integral13 + hsum_8x32(error13);
                errors[irow * error_width + icol + 14] = block1_out_integral + block1_in_integral14 + hsum_8x32(error14);
                errors[irow * error_width + icol + 15] = block1_out_integral + block1_in_integral15 + hsum_8x32(error15);
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
            for (int icol = 0; icol < error_width - 15; icol += 16)
            {
                __m256i error00, error00_block0, error00_block1, error00_block2, error00_block3, error00_block02, error00_block13;
                __m256i error01, error01_block0, error01_block1, error01_block2, error01_block3, error01_block02, error01_block13;
                __m256i error02, error02_block0, error02_block1, error02_block2, error02_block3, error02_block02, error02_block13;
                __m256i error03, error03_block0, error03_block1, error03_block2, error03_block3, error03_block02, error03_block13;
                __m256i error04, error04_block0, error04_block1, error04_block2, error04_block3, error04_block02, error04_block13;
                __m256i error05, error05_block0, error05_block1, error05_block2, error05_block3, error05_block02, error05_block13;
                __m256i error06, error06_block0, error06_block1, error06_block2, error06_block3, error06_block02, error06_block13;
                __m256i error07, error07_block0, error07_block1, error07_block2, error07_block3, error07_block02, error07_block13;
                __m256i error08, error08_block0, error08_block1, error08_block2, error08_block3, error08_block02, error08_block13;
                __m256i error09, error09_block0, error09_block1, error09_block2, error09_block3, error09_block02, error09_block13;
                __m256i error10, error10_block0, error10_block1, error10_block2, error10_block3, error10_block02, error10_block13;
                __m256i error11, error11_block0, error11_block1, error11_block2, error11_block3, error11_block02, error11_block13;
                __m256i error12, error12_block0, error12_block1, error12_block2, error12_block3, error12_block02, error12_block13;
                __m256i error13, error13_block0, error13_block1, error13_block2, error13_block3, error13_block02, error13_block13;
                __m256i error14, error14_block0, error14_block1, error14_block2, error14_block3, error14_block02, error14_block13;
                __m256i error15, error15_block0, error15_block1, error15_block2, error15_block3, error15_block02, error15_block13;
                error00 = error00_block0 = error00_block1 = error00_block2 = error00_block3 = error00_block02 = error00_block13 = zero;
                error01 = error01_block0 = error01_block1 = error01_block2 = error01_block3 = error01_block02 = error01_block13 = zero;
                error02 = error02_block0 = error02_block1 = error02_block2 = error02_block3 = error02_block02 = error02_block13 = zero;
                error03 = error03_block0 = error03_block1 = error03_block2 = error03_block3 = error03_block02 = error03_block13 = zero;
                error04 = error04_block0 = error04_block1 = error04_block2 = error04_block3 = error04_block02 = error04_block13 = zero;
                error05 = error05_block0 = error05_block1 = error05_block2 = error05_block3 = error05_block02 = error05_block13 = zero;
                error06 = error06_block0 = error06_block1 = error06_block2 = error06_block3 = error06_block02 = error06_block13 = zero;
                error07 = error07_block0 = error07_block1 = error07_block2 = error07_block3 = error07_block02 = error07_block13 = zero;
                error08 = error08_block0 = error08_block1 = error08_block2 = error08_block3 = error08_block02 = error08_block13 = zero;
                error09 = error09_block0 = error09_block1 = error09_block2 = error09_block3 = error09_block02 = error09_block13 = zero;
                error10 = error10_block0 = error10_block1 = error10_block2 = error10_block3 = error10_block02 = error10_block13 = zero;
                error11 = error11_block0 = error11_block1 = error11_block2 = error11_block3 = error11_block02 = error11_block13 = zero;
                error12 = error12_block0 = error12_block1 = error12_block2 = error12_block3 = error12_block02 = error12_block13 = zero;
                error13 = error13_block0 = error13_block1 = error13_block2 = error13_block3 = error13_block02 = error13_block13 = zero;
                error14 = error14_block0 = error14_block1 = error14_block2 = error14_block3 = error14_block02 = error14_block13 = zero;
                error15 = error15_block0 = error15_block1 = error15_block2 = error15_block3 = error15_block02 = error15_block13 = zero;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start00 = block1_in_integral_base + irow * integral_width + icol + 0;
                int block1_in_integral_start01 = block1_in_integral_base + irow * integral_width + icol + 1;
                int block1_in_integral_start02 = block1_in_integral_base + irow * integral_width + icol + 2;
                int block1_in_integral_start03 = block1_in_integral_base + irow * integral_width + icol + 3;
                int block1_in_integral_start04 = block1_in_integral_base + irow * integral_width + icol + 4;
                int block1_in_integral_start05 = block1_in_integral_base + irow * integral_width + icol + 5;
                int block1_in_integral_start06 = block1_in_integral_base + irow * integral_width + icol + 6;
                int block1_in_integral_start07 = block1_in_integral_base + irow * integral_width + icol + 7;
                int block1_in_integral_start08 = block1_in_integral_base + irow * integral_width + icol + 8;
                int block1_in_integral_start09 = block1_in_integral_base + irow * integral_width + icol + 9;
                int block1_in_integral_start10 = block1_in_integral_base + irow * integral_width + icol + 10;
                int block1_in_integral_start11 = block1_in_integral_base + irow * integral_width + icol + 11;
                int block1_in_integral_start12 = block1_in_integral_base + irow * integral_width + icol + 12;
                int block1_in_integral_start13 = block1_in_integral_base + irow * integral_width + icol + 13;
                int block1_in_integral_start14 = block1_in_integral_base + irow * integral_width + icol + 14;
                int block1_in_integral_start15 = block1_in_integral_base + irow * integral_width + icol + 15;
                error_t block1_in_integral00 = INTEGRAL(block1_in_integral_start00, height1, width1, integral_width);
                error_t block1_in_integral01 = INTEGRAL(block1_in_integral_start01, height1, width1, integral_width);
                error_t block1_in_integral02 = INTEGRAL(block1_in_integral_start02, height1, width1, integral_width);
                error_t block1_in_integral03 = INTEGRAL(block1_in_integral_start03, height1, width1, integral_width);
                error_t block1_in_integral04 = INTEGRAL(block1_in_integral_start04, height1, width1, integral_width);
                error_t block1_in_integral05 = INTEGRAL(block1_in_integral_start05, height1, width1, integral_width);
                error_t block1_in_integral06 = INTEGRAL(block1_in_integral_start06, height1, width1, integral_width);
                error_t block1_in_integral07 = INTEGRAL(block1_in_integral_start07, height1, width1, integral_width);
                error_t block1_in_integral08 = INTEGRAL(block1_in_integral_start08, height1, width1, integral_width);
                error_t block1_in_integral09 = INTEGRAL(block1_in_integral_start09, height1, width1, integral_width);
                error_t block1_in_integral10 = INTEGRAL(block1_in_integral_start10, height1, width1, integral_width);
                error_t block1_in_integral11 = INTEGRAL(block1_in_integral_start11, height1, width1, integral_width);
                error_t block1_in_integral12 = INTEGRAL(block1_in_integral_start12, height1, width1, integral_width);
                error_t block1_in_integral13 = INTEGRAL(block1_in_integral_start13, height1, width1, integral_width);
                error_t block1_in_integral14 = INTEGRAL(block1_in_integral_start14, height1, width1, integral_width);
                error_t block1_in_integral15 = INTEGRAL(block1_in_integral_start15, height1, width1, integral_width);

                int block3_in_integral_start00 = block3_in_integral_base + irow * integral_width + icol + 0;
                int block3_in_integral_start01 = block3_in_integral_base + irow * integral_width + icol + 1;
                int block3_in_integral_start02 = block3_in_integral_base + irow * integral_width + icol + 2;
                int block3_in_integral_start03 = block3_in_integral_base + irow * integral_width + icol + 3;
                int block3_in_integral_start04 = block3_in_integral_base + irow * integral_width + icol + 4;
                int block3_in_integral_start05 = block3_in_integral_base + irow * integral_width + icol + 5;
                int block3_in_integral_start06 = block3_in_integral_base + irow * integral_width + icol + 6;
                int block3_in_integral_start07 = block3_in_integral_base + irow * integral_width + icol + 7;
                int block3_in_integral_start08 = block3_in_integral_base + irow * integral_width + icol + 8;
                int block3_in_integral_start09 = block3_in_integral_base + irow * integral_width + icol + 9;
                int block3_in_integral_start10 = block3_in_integral_base + irow * integral_width + icol + 10;
                int block3_in_integral_start11 = block3_in_integral_base + irow * integral_width + icol + 11;
                int block3_in_integral_start12 = block3_in_integral_base + irow * integral_width + icol + 12;
                int block3_in_integral_start13 = block3_in_integral_base + irow * integral_width + icol + 13;
                int block3_in_integral_start14 = block3_in_integral_base + irow * integral_width + icol + 14;
                int block3_in_integral_start15 = block3_in_integral_base + irow * integral_width + icol + 15;
                error_t block3_in_integral00 = INTEGRAL(block3_in_integral_start00, height3, width3, integral_width);
                error_t block3_in_integral01 = INTEGRAL(block3_in_integral_start01, height3, width3, integral_width);
                error_t block3_in_integral02 = INTEGRAL(block3_in_integral_start02, height3, width3, integral_width);
                error_t block3_in_integral03 = INTEGRAL(block3_in_integral_start03, height3, width3, integral_width);
                error_t block3_in_integral04 = INTEGRAL(block3_in_integral_start04, height3, width3, integral_width);
                error_t block3_in_integral05 = INTEGRAL(block3_in_integral_start05, height3, width3, integral_width);
                error_t block3_in_integral06 = INTEGRAL(block3_in_integral_start06, height3, width3, integral_width);
                error_t block3_in_integral07 = INTEGRAL(block3_in_integral_start07, height3, width3, integral_width);
                error_t block3_in_integral08 = INTEGRAL(block3_in_integral_start08, height3, width3, integral_width);
                error_t block3_in_integral09 = INTEGRAL(block3_in_integral_start09, height3, width3, integral_width);
                error_t block3_in_integral10 = INTEGRAL(block3_in_integral_start10, height3, width3, integral_width);
                error_t block3_in_integral11 = INTEGRAL(block3_in_integral_start11, height3, width3, integral_width);
                error_t block3_in_integral12 = INTEGRAL(block3_in_integral_start12, height3, width3, integral_width);
                error_t block3_in_integral13 = INTEGRAL(block3_in_integral_start13, height3, width3, integral_width);
                error_t block3_in_integral14 = INTEGRAL(block3_in_integral_start14, height3, width3, integral_width);
                error_t block3_in_integral15 = INTEGRAL(block3_in_integral_start15, height3, width3, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block0 l2norm part1
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_madd_epi16(diff00, diff00));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_madd_epi16(diff01, diff01));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_madd_epi16(diff02, diff02));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_madd_epi16(diff03, diff03));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_madd_epi16(diff04, diff04));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_madd_epi16(diff05, diff05));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_madd_epi16(diff06, diff06));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_madd_epi16(diff07, diff07));
                        error08_block0 = _mm256_add_epi32(error08_block0, _mm256_madd_epi16(diff08, diff08));
                        error09_block0 = _mm256_add_epi32(error09_block0, _mm256_madd_epi16(diff09, diff09));
                        error10_block0 = _mm256_add_epi32(error10_block0, _mm256_madd_epi16(diff10, diff10));
                        error11_block0 = _mm256_add_epi32(error11_block0, _mm256_madd_epi16(diff11, diff11));
                        error12_block0 = _mm256_add_epi32(error12_block0, _mm256_madd_epi16(diff12, diff12));
                        error13_block0 = _mm256_add_epi32(error13_block0, _mm256_madd_epi16(diff13, diff13));
                        error14_block0 = _mm256_add_epi32(error14_block0, _mm256_madd_epi16(diff14, diff14));
                        error15_block0 = _mm256_add_epi32(error15_block0, _mm256_madd_epi16(diff15, diff15));
                    }
                    if (overlap % 16 != 0)
                    {
                        // block0 l2norm part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);
                        __m256i sqrdiff08 = _mm256_madd_epi16(diff08, diff08);
                        __m256i sqrdiff09 = _mm256_madd_epi16(diff09, diff09);
                        __m256i sqrdiff10 = _mm256_madd_epi16(diff10, diff10);
                        __m256i sqrdiff11 = _mm256_madd_epi16(diff11, diff11);
                        __m256i sqrdiff12 = _mm256_madd_epi16(diff12, diff12);
                        __m256i sqrdiff13 = _mm256_madd_epi16(diff13, diff13);
                        __m256i sqrdiff14 = _mm256_madd_epi16(diff14, diff14);
                        __m256i sqrdiff15 = _mm256_madd_epi16(diff15, diff15);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_blend_epi32(sqrdiff00, zero, 0b11110000));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_blend_epi32(sqrdiff01, zero, 0b11110000));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_blend_epi32(sqrdiff02, zero, 0b11110000));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_blend_epi32(sqrdiff03, zero, 0b11110000));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_blend_epi32(sqrdiff04, zero, 0b11110000));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_blend_epi32(sqrdiff05, zero, 0b11110000));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_blend_epi32(sqrdiff06, zero, 0b11110000));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_blend_epi32(sqrdiff07, zero, 0b11110000));
                        error08_block0 = _mm256_add_epi32(error08_block0, _mm256_blend_epi32(sqrdiff08, zero, 0b11110000));
                        error09_block0 = _mm256_add_epi32(error09_block0, _mm256_blend_epi32(sqrdiff09, zero, 0b11110000));
                        error10_block0 = _mm256_add_epi32(error10_block0, _mm256_blend_epi32(sqrdiff10, zero, 0b11110000));
                        error11_block0 = _mm256_add_epi32(error11_block0, _mm256_blend_epi32(sqrdiff11, zero, 0b11110000));
                        error12_block0 = _mm256_add_epi32(error12_block0, _mm256_blend_epi32(sqrdiff12, zero, 0b11110000));
                        error13_block0 = _mm256_add_epi32(error13_block0, _mm256_blend_epi32(sqrdiff13, zero, 0b11110000));
                        error14_block0 = _mm256_add_epi32(error14_block0, _mm256_blend_epi32(sqrdiff14, zero, 0b11110000));
                        error15_block0 = _mm256_add_epi32(error15_block0, _mm256_blend_epi32(sqrdiff15, zero, 0b11110000));

                        // block1 mulsum part1
                        // above_in = above-block input number 00
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b00001111));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b00001111));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b00001111));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b00001111));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b00001111));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b00001111));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b00001111));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b00001111));
                        error08_block1 = _mm256_add_epi32(error08_block1, _mm256_blend_epi32(mul08, zero, 0b00001111));
                        error09_block1 = _mm256_add_epi32(error09_block1, _mm256_blend_epi32(mul09, zero, 0b00001111));
                        error10_block1 = _mm256_add_epi32(error10_block1, _mm256_blend_epi32(mul10, zero, 0b00001111));
                        error11_block1 = _mm256_add_epi32(error11_block1, _mm256_blend_epi32(mul11, zero, 0b00001111));
                        error12_block1 = _mm256_add_epi32(error12_block1, _mm256_blend_epi32(mul12, zero, 0b00001111));
                        error13_block1 = _mm256_add_epi32(error13_block1, _mm256_blend_epi32(mul13, zero, 0b00001111));
                        error14_block1 = _mm256_add_epi32(error14_block1, _mm256_blend_epi32(mul14, zero, 0b00001111));
                        error15_block1 = _mm256_add_epi32(error15_block1, _mm256_blend_epi32(mul15, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < (blocksize - overlap) * 3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                        error01_block1 = _mm256_add_epi32(error01_block1, mul01);
                        error02_block1 = _mm256_add_epi32(error02_block1, mul02);
                        error03_block1 = _mm256_add_epi32(error03_block1, mul03);
                        error04_block1 = _mm256_add_epi32(error04_block1, mul04);
                        error05_block1 = _mm256_add_epi32(error05_block1, mul05);
                        error06_block1 = _mm256_add_epi32(error06_block1, mul06);
                        error07_block1 = _mm256_add_epi32(error07_block1, mul07);
                        error08_block1 = _mm256_add_epi32(error08_block1, mul08);
                        error09_block1 = _mm256_add_epi32(error09_block1, mul09);
                        error10_block1 = _mm256_add_epi32(error10_block1, mul10);
                        error11_block1 = _mm256_add_epi32(error11_block1, mul11);
                        error12_block1 = _mm256_add_epi32(error12_block1, mul12);
                        error13_block1 = _mm256_add_epi32(error13_block1, mul13);
                        error14_block1 = _mm256_add_epi32(error14_block1, mul14);
                        error15_block1 = _mm256_add_epi32(error15_block1, mul15);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block1 mulsum part3
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b11110000));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b11110000));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b11110000));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b11110000));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b11110000));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b11110000));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b11110000));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b11110000));
                        error08_block1 = _mm256_add_epi32(error08_block1, _mm256_blend_epi32(mul08, zero, 0b11110000));
                        error09_block1 = _mm256_add_epi32(error09_block1, _mm256_blend_epi32(mul09, zero, 0b11110000));
                        error10_block1 = _mm256_add_epi32(error10_block1, _mm256_blend_epi32(mul10, zero, 0b11110000));
                        error11_block1 = _mm256_add_epi32(error11_block1, _mm256_blend_epi32(mul11, zero, 0b11110000));
                        error12_block1 = _mm256_add_epi32(error12_block1, _mm256_blend_epi32(mul12, zero, 0b11110000));
                        error13_block1 = _mm256_add_epi32(error13_block1, _mm256_blend_epi32(mul13, zero, 0b11110000));
                        error14_block1 = _mm256_add_epi32(error14_block1, _mm256_blend_epi32(mul14, zero, 0b11110000));
                        error15_block1 = _mm256_add_epi32(error15_block1, _mm256_blend_epi32(mul15, zero, 0b11110000));

                        // block2 l2norm part1
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);
                        __m256i sqrdiff08 = _mm256_madd_epi16(diff08, diff08);
                        __m256i sqrdiff09 = _mm256_madd_epi16(diff09, diff09);
                        __m256i sqrdiff10 = _mm256_madd_epi16(diff10, diff10);
                        __m256i sqrdiff11 = _mm256_madd_epi16(diff11, diff11);
                        __m256i sqrdiff12 = _mm256_madd_epi16(diff12, diff12);
                        __m256i sqrdiff13 = _mm256_madd_epi16(diff13, diff13);
                        __m256i sqrdiff14 = _mm256_madd_epi16(diff14, diff14);
                        __m256i sqrdiff15 = _mm256_madd_epi16(diff15, diff15);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_blend_epi32(sqrdiff00, zero, 0b00001111));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_blend_epi32(sqrdiff01, zero, 0b00001111));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_blend_epi32(sqrdiff02, zero, 0b00001111));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_blend_epi32(sqrdiff03, zero, 0b00001111));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_blend_epi32(sqrdiff04, zero, 0b00001111));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_blend_epi32(sqrdiff05, zero, 0b00001111));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_blend_epi32(sqrdiff06, zero, 0b00001111));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_blend_epi32(sqrdiff07, zero, 0b00001111));
                        error08_block2 = _mm256_add_epi32(error08_block2, _mm256_blend_epi32(sqrdiff08, zero, 0b00001111));
                        error09_block2 = _mm256_add_epi32(error09_block2, _mm256_blend_epi32(sqrdiff09, zero, 0b00001111));
                        error10_block2 = _mm256_add_epi32(error10_block2, _mm256_blend_epi32(sqrdiff10, zero, 0b00001111));
                        error11_block2 = _mm256_add_epi32(error11_block2, _mm256_blend_epi32(sqrdiff11, zero, 0b00001111));
                        error12_block2 = _mm256_add_epi32(error12_block2, _mm256_blend_epi32(sqrdiff12, zero, 0b00001111));
                        error13_block2 = _mm256_add_epi32(error13_block2, _mm256_blend_epi32(sqrdiff13, zero, 0b00001111));
                        error14_block2 = _mm256_add_epi32(error14_block2, _mm256_blend_epi32(sqrdiff14, zero, 0b00001111));
                        error15_block2 = _mm256_add_epi32(error15_block2, _mm256_blend_epi32(sqrdiff15, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < blocksize * 3 - 15; m += 16)
                    {
                        // block2 l2norm part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        error00_block2 = _mm256_add_epi32(error00_block2, _mm256_madd_epi16(diff00, diff00));
                        error01_block2 = _mm256_add_epi32(error01_block2, _mm256_madd_epi16(diff01, diff01));
                        error02_block2 = _mm256_add_epi32(error02_block2, _mm256_madd_epi16(diff02, diff02));
                        error03_block2 = _mm256_add_epi32(error03_block2, _mm256_madd_epi16(diff03, diff03));
                        error04_block2 = _mm256_add_epi32(error04_block2, _mm256_madd_epi16(diff04, diff04));
                        error05_block2 = _mm256_add_epi32(error05_block2, _mm256_madd_epi16(diff05, diff05));
                        error06_block2 = _mm256_add_epi32(error06_block2, _mm256_madd_epi16(diff06, diff06));
                        error07_block2 = _mm256_add_epi32(error07_block2, _mm256_madd_epi16(diff07, diff07));
                        error08_block2 = _mm256_add_epi32(error08_block2, _mm256_madd_epi16(diff08, diff08));
                        error09_block2 = _mm256_add_epi32(error09_block2, _mm256_madd_epi16(diff09, diff09));
                        error10_block2 = _mm256_add_epi32(error10_block2, _mm256_madd_epi16(diff10, diff10));
                        error11_block2 = _mm256_add_epi32(error11_block2, _mm256_madd_epi16(diff11, diff11));
                        error12_block2 = _mm256_add_epi32(error12_block2, _mm256_madd_epi16(diff12, diff12));
                        error13_block2 = _mm256_add_epi32(error13_block2, _mm256_madd_epi16(diff13, diff13));
                        error14_block2 = _mm256_add_epi32(error14_block2, _mm256_madd_epi16(diff14, diff14));
                        error15_block2 = _mm256_add_epi32(error15_block2, _mm256_madd_epi16(diff15, diff15));
                    }
                }
                // we start with row overlap to ensure equivalence in handling block 3 (to other cases)
                for (int k = overlap; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block3 mulsum part1
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i left_in = loadu_8x16(&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);
                        __m256i mull08 = _mm256_madd_epi16(curr_in_08, left_in);
                        __m256i mull09 = _mm256_madd_epi16(curr_in_09, left_in);
                        __m256i mull10 = _mm256_madd_epi16(curr_in_10, left_in);
                        __m256i mull11 = _mm256_madd_epi16(curr_in_11, left_in);
                        __m256i mull12 = _mm256_madd_epi16(curr_in_12, left_in);
                        __m256i mull13 = _mm256_madd_epi16(curr_in_13, left_in);
                        __m256i mull14 = _mm256_madd_epi16(curr_in_14, left_in);
                        __m256i mull15 = _mm256_madd_epi16(curr_in_15, left_in);
                        error00_block3 = _mm256_add_epi32(mull00, error00_block3);
                        error01_block3 = _mm256_add_epi32(mull01, error01_block3);
                        error02_block3 = _mm256_add_epi32(mull02, error02_block3);
                        error03_block3 = _mm256_add_epi32(mull03, error03_block3);
                        error04_block3 = _mm256_add_epi32(mull04, error04_block3);
                        error05_block3 = _mm256_add_epi32(mull05, error05_block3);
                        error06_block3 = _mm256_add_epi32(mull06, error06_block3);
                        error07_block3 = _mm256_add_epi32(mull07, error07_block3);
                        error08_block3 = _mm256_add_epi32(mull08, error08_block3);
                        error09_block3 = _mm256_add_epi32(mull09, error09_block3);
                        error10_block3 = _mm256_add_epi32(mull10, error10_block3);
                        error11_block3 = _mm256_add_epi32(mull11, error11_block3);
                        error12_block3 = _mm256_add_epi32(mull12, error12_block3);
                        error13_block3 = _mm256_add_epi32(mull13, error13_block3);
                        error14_block3 = _mm256_add_epi32(mull14, error14_block3);
                        error15_block3 = _mm256_add_epi32(mull15, error15_block3);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block3 mulsum part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i left_in = loadu_8x16(&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);
                        __m256i mull08 = _mm256_madd_epi16(curr_in_08, left_in);
                        __m256i mull09 = _mm256_madd_epi16(curr_in_09, left_in);
                        __m256i mull10 = _mm256_madd_epi16(curr_in_10, left_in);
                        __m256i mull11 = _mm256_madd_epi16(curr_in_11, left_in);
                        __m256i mull12 = _mm256_madd_epi16(curr_in_12, left_in);
                        __m256i mull13 = _mm256_madd_epi16(curr_in_13, left_in);
                        __m256i mull14 = _mm256_madd_epi16(curr_in_14, left_in);
                        __m256i mull15 = _mm256_madd_epi16(curr_in_15, left_in);
                        error00_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00_block3);
                        error01_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull01, zero, 0b11110000), error01_block3);
                        error02_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull02, zero, 0b11110000), error02_block3);
                        error03_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull03, zero, 0b11110000), error03_block3);
                        error04_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull04, zero, 0b11110000), error04_block3);
                        error05_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull05, zero, 0b11110000), error05_block3);
                        error06_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull06, zero, 0b11110000), error06_block3);
                        error07_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull07, zero, 0b11110000), error07_block3);
                        error08_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull08, zero, 0b11110000), error08_block3);
                        error09_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull09, zero, 0b11110000), error09_block3);
                        error10_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull10, zero, 0b11110000), error10_block3);
                        error11_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull11, zero, 0b11110000), error11_block3);
                        error12_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull12, zero, 0b11110000), error12_block3);
                        error13_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull13, zero, 0b11110000), error13_block3);
                        error14_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull14, zero, 0b11110000), error14_block3);
                        error15_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull15, zero, 0b11110000), error15_block3);
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
                error08_block02 = _mm256_add_epi32(error08_block0, error08_block2);
                error09_block02 = _mm256_add_epi32(error09_block0, error09_block2);
                error10_block02 = _mm256_add_epi32(error10_block0, error10_block2);
                error11_block02 = _mm256_add_epi32(error11_block0, error11_block2);
                error12_block02 = _mm256_add_epi32(error12_block0, error12_block2);
                error13_block02 = _mm256_add_epi32(error13_block0, error13_block2);
                error14_block02 = _mm256_add_epi32(error14_block0, error14_block2);
                error15_block02 = _mm256_add_epi32(error15_block0, error15_block2);
                error00_block13 = _mm256_slli_epi32(_mm256_add_epi32(error00_block1, error00_block3), 1);
                error01_block13 = _mm256_slli_epi32(_mm256_add_epi32(error01_block1, error01_block3), 1);
                error02_block13 = _mm256_slli_epi32(_mm256_add_epi32(error02_block1, error02_block3), 1);
                error03_block13 = _mm256_slli_epi32(_mm256_add_epi32(error03_block1, error03_block3), 1);
                error04_block13 = _mm256_slli_epi32(_mm256_add_epi32(error04_block1, error04_block3), 1);
                error05_block13 = _mm256_slli_epi32(_mm256_add_epi32(error05_block1, error05_block3), 1);
                error06_block13 = _mm256_slli_epi32(_mm256_add_epi32(error06_block1, error06_block3), 1);
                error07_block13 = _mm256_slli_epi32(_mm256_add_epi32(error07_block1, error07_block3), 1);
                error08_block13 = _mm256_slli_epi32(_mm256_add_epi32(error08_block1, error08_block3), 1);
                error09_block13 = _mm256_slli_epi32(_mm256_add_epi32(error09_block1, error09_block3), 1);
                error10_block13 = _mm256_slli_epi32(_mm256_add_epi32(error10_block1, error10_block3), 1);
                error11_block13 = _mm256_slli_epi32(_mm256_add_epi32(error11_block1, error11_block3), 1);
                error12_block13 = _mm256_slli_epi32(_mm256_add_epi32(error12_block1, error12_block3), 1);
                error13_block13 = _mm256_slli_epi32(_mm256_add_epi32(error13_block1, error13_block3), 1);
                error14_block13 = _mm256_slli_epi32(_mm256_add_epi32(error14_block1, error14_block3), 1);
                error15_block13 = _mm256_slli_epi32(_mm256_add_epi32(error15_block1, error15_block3), 1);
                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                error01 = _mm256_sub_epi32(error01_block02, error01_block13);
                error02 = _mm256_sub_epi32(error02_block02, error02_block13);
                error03 = _mm256_sub_epi32(error03_block02, error03_block13);
                error04 = _mm256_sub_epi32(error04_block02, error04_block13);
                error05 = _mm256_sub_epi32(error05_block02, error05_block13);
                error06 = _mm256_sub_epi32(error06_block02, error06_block13);
                error07 = _mm256_sub_epi32(error07_block02, error07_block13);
                error08 = _mm256_sub_epi32(error08_block02, error08_block13);
                error09 = _mm256_sub_epi32(error09_block02, error09_block13);
                error10 = _mm256_sub_epi32(error10_block02, error10_block13);
                error11 = _mm256_sub_epi32(error11_block02, error11_block13);
                error12 = _mm256_sub_epi32(error12_block02, error12_block13);
                error13 = _mm256_sub_epi32(error13_block02, error13_block13);
                error14 = _mm256_sub_epi32(error14_block02, error14_block13);
                error15 = _mm256_sub_epi32(error15_block02, error15_block13);
                errors[irow * error_width + icol + 0] = block1_out_integral + block1_in_integral00 + hsum_8x32(error00) + block3_out_integral + block3_in_integral00;
                errors[irow * error_width + icol + 1] = block1_out_integral + block1_in_integral01 + hsum_8x32(error01) + block3_out_integral + block3_in_integral01;
                errors[irow * error_width + icol + 2] = block1_out_integral + block1_in_integral02 + hsum_8x32(error02) + block3_out_integral + block3_in_integral02;
                errors[irow * error_width + icol + 3] = block1_out_integral + block1_in_integral03 + hsum_8x32(error03) + block3_out_integral + block3_in_integral03;
                errors[irow * error_width + icol + 4] = block1_out_integral + block1_in_integral04 + hsum_8x32(error04) + block3_out_integral + block3_in_integral04;
                errors[irow * error_width + icol + 5] = block1_out_integral + block1_in_integral05 + hsum_8x32(error05) + block3_out_integral + block3_in_integral05;
                errors[irow * error_width + icol + 6] = block1_out_integral + block1_in_integral06 + hsum_8x32(error06) + block3_out_integral + block3_in_integral06;
                errors[irow * error_width + icol + 7] = block1_out_integral + block1_in_integral07 + hsum_8x32(error07) + block3_out_integral + block3_in_integral07;
                errors[irow * error_width + icol + 8] = block1_out_integral + block1_in_integral08 + hsum_8x32(error08) + block3_out_integral + block3_in_integral08;
                errors[irow * error_width + icol + 9] = block1_out_integral + block1_in_integral09 + hsum_8x32(error09) + block3_out_integral + block3_in_integral09;
                errors[irow * error_width + icol + 10] = block1_out_integral + block1_in_integral10 + hsum_8x32(error10) + block3_out_integral + block3_in_integral10;
                errors[irow * error_width + icol + 11] = block1_out_integral + block1_in_integral11 + hsum_8x32(error11) + block3_out_integral + block3_in_integral11;
                errors[irow * error_width + icol + 12] = block1_out_integral + block1_in_integral12 + hsum_8x32(error12) + block3_out_integral + block3_in_integral12;
                errors[irow * error_width + icol + 13] = block1_out_integral + block1_in_integral13 + hsum_8x32(error13) + block3_out_integral + block3_in_integral13;
                errors[irow * error_width + icol + 14] = block1_out_integral + block1_in_integral14 + hsum_8x32(error14) + block3_out_integral + block3_in_integral14;
                errors[irow * error_width + icol + 15] = block1_out_integral + block1_in_integral15 + hsum_8x32(error15) + block3_out_integral + block3_in_integral15;
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
            for (int icol = 0; icol < error_width - 15; icol += 16)
            {
                __m256i error00, error00_block0, error00_block1, error00_block3, error00_block02, error00_block13;
                __m256i error01, error01_block0, error01_block1, error01_block3, error01_block02, error01_block13;
                __m256i error02, error02_block0, error02_block1, error02_block3, error02_block02, error02_block13;
                __m256i error03, error03_block0, error03_block1, error03_block3, error03_block02, error03_block13;
                __m256i error04, error04_block0, error04_block1, error04_block3, error04_block02, error04_block13;
                __m256i error05, error05_block0, error05_block1, error05_block3, error05_block02, error05_block13;
                __m256i error06, error06_block0, error06_block1, error06_block3, error06_block02, error06_block13;
                __m256i error07, error07_block0, error07_block1, error07_block3, error07_block02, error07_block13;
                __m256i error08, error08_block0, error08_block1, error08_block3, error08_block02, error08_block13;
                __m256i error09, error09_block0, error09_block1, error09_block3, error09_block02, error09_block13;
                __m256i error10, error10_block0, error10_block1, error10_block3, error10_block02, error10_block13;
                __m256i error11, error11_block0, error11_block1, error11_block3, error11_block02, error11_block13;
                __m256i error12, error12_block0, error12_block1, error12_block3, error12_block02, error12_block13;
                __m256i error13, error13_block0, error13_block1, error13_block3, error13_block02, error13_block13;
                __m256i error14, error14_block0, error14_block1, error14_block3, error14_block02, error14_block13;
                __m256i error15, error15_block0, error15_block1, error15_block3, error15_block02, error15_block13;
                error00 = error00_block0 = error00_block1 = error00_block3 = error00_block02 = error00_block13 = zero;
                error01 = error01_block0 = error01_block1 = error01_block3 = error01_block02 = error01_block13 = zero;
                error02 = error02_block0 = error02_block1 = error02_block3 = error02_block02 = error02_block13 = zero;
                error03 = error03_block0 = error03_block1 = error03_block3 = error03_block02 = error03_block13 = zero;
                error04 = error04_block0 = error04_block1 = error04_block3 = error04_block02 = error04_block13 = zero;
                error05 = error05_block0 = error05_block1 = error05_block3 = error05_block02 = error05_block13 = zero;
                error06 = error06_block0 = error06_block1 = error06_block3 = error06_block02 = error06_block13 = zero;
                error07 = error07_block0 = error07_block1 = error07_block3 = error07_block02 = error07_block13 = zero;
                error08 = error08_block0 = error08_block1 = error08_block3 = error08_block02 = error08_block13 = zero;
                error09 = error09_block0 = error09_block1 = error09_block3 = error09_block02 = error09_block13 = zero;
                error10 = error10_block0 = error10_block1 = error10_block3 = error10_block02 = error10_block13 = zero;
                error11 = error11_block0 = error11_block1 = error11_block3 = error11_block02 = error11_block13 = zero;
                error12 = error12_block0 = error12_block1 = error12_block3 = error12_block02 = error12_block13 = zero;
                error13 = error13_block0 = error13_block1 = error13_block3 = error13_block02 = error13_block13 = zero;
                error14 = error14_block0 = error14_block1 = error14_block3 = error14_block02 = error14_block13 = zero;
                error15 = error15_block0 = error15_block1 = error15_block3 = error15_block02 = error15_block13 = zero;
                /* calculation of block 1 integral part in-variant */
                int block1_in_integral_start00 = block1_in_integral_base + irow * integral_width + icol + 0;
                int block1_in_integral_start01 = block1_in_integral_base + irow * integral_width + icol + 1;
                int block1_in_integral_start02 = block1_in_integral_base + irow * integral_width + icol + 2;
                int block1_in_integral_start03 = block1_in_integral_base + irow * integral_width + icol + 3;
                int block1_in_integral_start04 = block1_in_integral_base + irow * integral_width + icol + 4;
                int block1_in_integral_start05 = block1_in_integral_base + irow * integral_width + icol + 5;
                int block1_in_integral_start06 = block1_in_integral_base + irow * integral_width + icol + 6;
                int block1_in_integral_start07 = block1_in_integral_base + irow * integral_width + icol + 7;
                int block1_in_integral_start08 = block1_in_integral_base + irow * integral_width + icol + 8;
                int block1_in_integral_start09 = block1_in_integral_base + irow * integral_width + icol + 9;
                int block1_in_integral_start10 = block1_in_integral_base + irow * integral_width + icol + 10;
                int block1_in_integral_start11 = block1_in_integral_base + irow * integral_width + icol + 11;
                int block1_in_integral_start12 = block1_in_integral_base + irow * integral_width + icol + 12;
                int block1_in_integral_start13 = block1_in_integral_base + irow * integral_width + icol + 13;
                int block1_in_integral_start14 = block1_in_integral_base + irow * integral_width + icol + 14;
                int block1_in_integral_start15 = block1_in_integral_base + irow * integral_width + icol + 15;
                error_t block1_in_integral00 = INTEGRAL(block1_in_integral_start00, height1, width1, integral_width);
                error_t block1_in_integral01 = INTEGRAL(block1_in_integral_start01, height1, width1, integral_width);
                error_t block1_in_integral02 = INTEGRAL(block1_in_integral_start02, height1, width1, integral_width);
                error_t block1_in_integral03 = INTEGRAL(block1_in_integral_start03, height1, width1, integral_width);
                error_t block1_in_integral04 = INTEGRAL(block1_in_integral_start04, height1, width1, integral_width);
                error_t block1_in_integral05 = INTEGRAL(block1_in_integral_start05, height1, width1, integral_width);
                error_t block1_in_integral06 = INTEGRAL(block1_in_integral_start06, height1, width1, integral_width);
                error_t block1_in_integral07 = INTEGRAL(block1_in_integral_start07, height1, width1, integral_width);
                error_t block1_in_integral08 = INTEGRAL(block1_in_integral_start08, height1, width1, integral_width);
                error_t block1_in_integral09 = INTEGRAL(block1_in_integral_start09, height1, width1, integral_width);
                error_t block1_in_integral10 = INTEGRAL(block1_in_integral_start10, height1, width1, integral_width);
                error_t block1_in_integral11 = INTEGRAL(block1_in_integral_start11, height1, width1, integral_width);
                error_t block1_in_integral12 = INTEGRAL(block1_in_integral_start12, height1, width1, integral_width);
                error_t block1_in_integral13 = INTEGRAL(block1_in_integral_start13, height1, width1, integral_width);
                error_t block1_in_integral14 = INTEGRAL(block1_in_integral_start14, height1, width1, integral_width);
                error_t block1_in_integral15 = INTEGRAL(block1_in_integral_start15, height1, width1, integral_width);

                int block3_in_integral_start00 = block3_in_integral_base + irow * integral_width + icol + 0;
                int block3_in_integral_start01 = block3_in_integral_base + irow * integral_width + icol + 1;
                int block3_in_integral_start02 = block3_in_integral_base + irow * integral_width + icol + 2;
                int block3_in_integral_start03 = block3_in_integral_base + irow * integral_width + icol + 3;
                int block3_in_integral_start04 = block3_in_integral_base + irow * integral_width + icol + 4;
                int block3_in_integral_start05 = block3_in_integral_base + irow * integral_width + icol + 5;
                int block3_in_integral_start06 = block3_in_integral_base + irow * integral_width + icol + 6;
                int block3_in_integral_start07 = block3_in_integral_base + irow * integral_width + icol + 7;
                int block3_in_integral_start08 = block3_in_integral_base + irow * integral_width + icol + 8;
                int block3_in_integral_start09 = block3_in_integral_base + irow * integral_width + icol + 9;
                int block3_in_integral_start10 = block3_in_integral_base + irow * integral_width + icol + 10;
                int block3_in_integral_start11 = block3_in_integral_base + irow * integral_width + icol + 11;
                int block3_in_integral_start12 = block3_in_integral_base + irow * integral_width + icol + 12;
                int block3_in_integral_start13 = block3_in_integral_base + irow * integral_width + icol + 13;
                int block3_in_integral_start14 = block3_in_integral_base + irow * integral_width + icol + 14;
                int block3_in_integral_start15 = block3_in_integral_base + irow * integral_width + icol + 15;
                error_t block3_in_integral00 = INTEGRAL(block3_in_integral_start00, height3, width3, integral_width);
                error_t block3_in_integral01 = INTEGRAL(block3_in_integral_start01, height3, width3, integral_width);
                error_t block3_in_integral02 = INTEGRAL(block3_in_integral_start02, height3, width3, integral_width);
                error_t block3_in_integral03 = INTEGRAL(block3_in_integral_start03, height3, width3, integral_width);
                error_t block3_in_integral04 = INTEGRAL(block3_in_integral_start04, height3, width3, integral_width);
                error_t block3_in_integral05 = INTEGRAL(block3_in_integral_start05, height3, width3, integral_width);
                error_t block3_in_integral06 = INTEGRAL(block3_in_integral_start06, height3, width3, integral_width);
                error_t block3_in_integral07 = INTEGRAL(block3_in_integral_start07, height3, width3, integral_width);
                error_t block3_in_integral08 = INTEGRAL(block3_in_integral_start08, height3, width3, integral_width);
                error_t block3_in_integral09 = INTEGRAL(block3_in_integral_start09, height3, width3, integral_width);
                error_t block3_in_integral10 = INTEGRAL(block3_in_integral_start10, height3, width3, integral_width);
                error_t block3_in_integral11 = INTEGRAL(block3_in_integral_start11, height3, width3, integral_width);
                error_t block3_in_integral12 = INTEGRAL(block3_in_integral_start12, height3, width3, integral_width);
                error_t block3_in_integral13 = INTEGRAL(block3_in_integral_start13, height3, width3, integral_width);
                error_t block3_in_integral14 = INTEGRAL(block3_in_integral_start14, height3, width3, integral_width);
                error_t block3_in_integral15 = INTEGRAL(block3_in_integral_start15, height3, width3, integral_width);

                for (int k = 0; k < overlap; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block0 l2norm part1
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_madd_epi16(diff00, diff00));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_madd_epi16(diff01, diff01));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_madd_epi16(diff02, diff02));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_madd_epi16(diff03, diff03));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_madd_epi16(diff04, diff04));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_madd_epi16(diff05, diff05));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_madd_epi16(diff06, diff06));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_madd_epi16(diff07, diff07));
                        error08_block0 = _mm256_add_epi32(error08_block0, _mm256_madd_epi16(diff08, diff08));
                        error09_block0 = _mm256_add_epi32(error09_block0, _mm256_madd_epi16(diff09, diff09));
                        error10_block0 = _mm256_add_epi32(error10_block0, _mm256_madd_epi16(diff10, diff10));
                        error11_block0 = _mm256_add_epi32(error11_block0, _mm256_madd_epi16(diff11, diff11));
                        error12_block0 = _mm256_add_epi32(error12_block0, _mm256_madd_epi16(diff12, diff12));
                        error13_block0 = _mm256_add_epi32(error13_block0, _mm256_madd_epi16(diff13, diff13));
                        error14_block0 = _mm256_add_epi32(error14_block0, _mm256_madd_epi16(diff14, diff14));
                        error15_block0 = _mm256_add_epi32(error15_block0, _mm256_madd_epi16(diff15, diff15));
                    }
                    if (overlap % 16 != 0)
                    {
                        // block0 l2norm part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i curr_out = loadu_8x16(&out[(orow + k) * ojump + (ocol)*3 + m]);
                        __m256i diff00 = _mm256_sub_epi16(curr_in_00, curr_out);
                        __m256i diff01 = _mm256_sub_epi16(curr_in_01, curr_out);
                        __m256i diff02 = _mm256_sub_epi16(curr_in_02, curr_out);
                        __m256i diff03 = _mm256_sub_epi16(curr_in_03, curr_out);
                        __m256i diff04 = _mm256_sub_epi16(curr_in_04, curr_out);
                        __m256i diff05 = _mm256_sub_epi16(curr_in_05, curr_out);
                        __m256i diff06 = _mm256_sub_epi16(curr_in_06, curr_out);
                        __m256i diff07 = _mm256_sub_epi16(curr_in_07, curr_out);
                        __m256i diff08 = _mm256_sub_epi16(curr_in_08, curr_out);
                        __m256i diff09 = _mm256_sub_epi16(curr_in_09, curr_out);
                        __m256i diff10 = _mm256_sub_epi16(curr_in_10, curr_out);
                        __m256i diff11 = _mm256_sub_epi16(curr_in_11, curr_out);
                        __m256i diff12 = _mm256_sub_epi16(curr_in_12, curr_out);
                        __m256i diff13 = _mm256_sub_epi16(curr_in_13, curr_out);
                        __m256i diff14 = _mm256_sub_epi16(curr_in_14, curr_out);
                        __m256i diff15 = _mm256_sub_epi16(curr_in_15, curr_out);
                        __m256i sqrdiff00 = _mm256_madd_epi16(diff00, diff00);
                        __m256i sqrdiff01 = _mm256_madd_epi16(diff01, diff01);
                        __m256i sqrdiff02 = _mm256_madd_epi16(diff02, diff02);
                        __m256i sqrdiff03 = _mm256_madd_epi16(diff03, diff03);
                        __m256i sqrdiff04 = _mm256_madd_epi16(diff04, diff04);
                        __m256i sqrdiff05 = _mm256_madd_epi16(diff05, diff05);
                        __m256i sqrdiff06 = _mm256_madd_epi16(diff06, diff06);
                        __m256i sqrdiff07 = _mm256_madd_epi16(diff07, diff07);
                        __m256i sqrdiff08 = _mm256_madd_epi16(diff08, diff08);
                        __m256i sqrdiff09 = _mm256_madd_epi16(diff09, diff09);
                        __m256i sqrdiff10 = _mm256_madd_epi16(diff10, diff10);
                        __m256i sqrdiff11 = _mm256_madd_epi16(diff11, diff11);
                        __m256i sqrdiff12 = _mm256_madd_epi16(diff12, diff12);
                        __m256i sqrdiff13 = _mm256_madd_epi16(diff13, diff13);
                        __m256i sqrdiff14 = _mm256_madd_epi16(diff14, diff14);
                        __m256i sqrdiff15 = _mm256_madd_epi16(diff15, diff15);
                        error00_block0 = _mm256_add_epi32(error00_block0, _mm256_blend_epi32(sqrdiff00, zero, 0b11110000));
                        error01_block0 = _mm256_add_epi32(error01_block0, _mm256_blend_epi32(sqrdiff01, zero, 0b11110000));
                        error02_block0 = _mm256_add_epi32(error02_block0, _mm256_blend_epi32(sqrdiff02, zero, 0b11110000));
                        error03_block0 = _mm256_add_epi32(error03_block0, _mm256_blend_epi32(sqrdiff03, zero, 0b11110000));
                        error04_block0 = _mm256_add_epi32(error04_block0, _mm256_blend_epi32(sqrdiff04, zero, 0b11110000));
                        error05_block0 = _mm256_add_epi32(error05_block0, _mm256_blend_epi32(sqrdiff05, zero, 0b11110000));
                        error06_block0 = _mm256_add_epi32(error06_block0, _mm256_blend_epi32(sqrdiff06, zero, 0b11110000));
                        error07_block0 = _mm256_add_epi32(error07_block0, _mm256_blend_epi32(sqrdiff07, zero, 0b11110000));
                        error08_block0 = _mm256_add_epi32(error08_block0, _mm256_blend_epi32(sqrdiff08, zero, 0b11110000));
                        error09_block0 = _mm256_add_epi32(error09_block0, _mm256_blend_epi32(sqrdiff09, zero, 0b11110000));
                        error10_block0 = _mm256_add_epi32(error10_block0, _mm256_blend_epi32(sqrdiff10, zero, 0b11110000));
                        error11_block0 = _mm256_add_epi32(error11_block0, _mm256_blend_epi32(sqrdiff11, zero, 0b11110000));
                        error12_block0 = _mm256_add_epi32(error12_block0, _mm256_blend_epi32(sqrdiff12, zero, 0b11110000));
                        error13_block0 = _mm256_add_epi32(error13_block0, _mm256_blend_epi32(sqrdiff13, zero, 0b11110000));
                        error14_block0 = _mm256_add_epi32(error14_block0, _mm256_blend_epi32(sqrdiff14, zero, 0b11110000));
                        error15_block0 = _mm256_add_epi32(error15_block0, _mm256_blend_epi32(sqrdiff15, zero, 0b11110000));

                        // block1 mulsum part1
                        // above_in = above-block input number 00
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, _mm256_blend_epi32(mul00, zero, 0b00001111));
                        error01_block1 = _mm256_add_epi32(error01_block1, _mm256_blend_epi32(mul01, zero, 0b00001111));
                        error02_block1 = _mm256_add_epi32(error02_block1, _mm256_blend_epi32(mul02, zero, 0b00001111));
                        error03_block1 = _mm256_add_epi32(error03_block1, _mm256_blend_epi32(mul03, zero, 0b00001111));
                        error04_block1 = _mm256_add_epi32(error04_block1, _mm256_blend_epi32(mul04, zero, 0b00001111));
                        error05_block1 = _mm256_add_epi32(error05_block1, _mm256_blend_epi32(mul05, zero, 0b00001111));
                        error06_block1 = _mm256_add_epi32(error06_block1, _mm256_blend_epi32(mul06, zero, 0b00001111));
                        error07_block1 = _mm256_add_epi32(error07_block1, _mm256_blend_epi32(mul07, zero, 0b00001111));
                        error08_block1 = _mm256_add_epi32(error08_block1, _mm256_blend_epi32(mul08, zero, 0b00001111));
                        error09_block1 = _mm256_add_epi32(error09_block1, _mm256_blend_epi32(mul09, zero, 0b00001111));
                        error10_block1 = _mm256_add_epi32(error10_block1, _mm256_blend_epi32(mul10, zero, 0b00001111));
                        error11_block1 = _mm256_add_epi32(error11_block1, _mm256_blend_epi32(mul11, zero, 0b00001111));
                        error12_block1 = _mm256_add_epi32(error12_block1, _mm256_blend_epi32(mul12, zero, 0b00001111));
                        error13_block1 = _mm256_add_epi32(error13_block1, _mm256_blend_epi32(mul13, zero, 0b00001111));
                        error14_block1 = _mm256_add_epi32(error14_block1, _mm256_blend_epi32(mul14, zero, 0b00001111));
                        error15_block1 = _mm256_add_epi32(error15_block1, _mm256_blend_epi32(mul15, zero, 0b00001111));
                        m += 16;
                    }
                    for (; m < (blocksize)*3 - 15; m += 16)
                    {
                        // block1 mulsum part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i above_in = loadu_8x16(&in[(arow + (blocksize - overlap) + k) * ijump + (acol)*3 + m]);
                        __m256i mul00 = _mm256_madd_epi16(curr_in_00, above_in);
                        __m256i mul01 = _mm256_madd_epi16(curr_in_01, above_in);
                        __m256i mul02 = _mm256_madd_epi16(curr_in_02, above_in);
                        __m256i mul03 = _mm256_madd_epi16(curr_in_03, above_in);
                        __m256i mul04 = _mm256_madd_epi16(curr_in_04, above_in);
                        __m256i mul05 = _mm256_madd_epi16(curr_in_05, above_in);
                        __m256i mul06 = _mm256_madd_epi16(curr_in_06, above_in);
                        __m256i mul07 = _mm256_madd_epi16(curr_in_07, above_in);
                        __m256i mul08 = _mm256_madd_epi16(curr_in_08, above_in);
                        __m256i mul09 = _mm256_madd_epi16(curr_in_09, above_in);
                        __m256i mul10 = _mm256_madd_epi16(curr_in_10, above_in);
                        __m256i mul11 = _mm256_madd_epi16(curr_in_11, above_in);
                        __m256i mul12 = _mm256_madd_epi16(curr_in_12, above_in);
                        __m256i mul13 = _mm256_madd_epi16(curr_in_13, above_in);
                        __m256i mul14 = _mm256_madd_epi16(curr_in_14, above_in);
                        __m256i mul15 = _mm256_madd_epi16(curr_in_15, above_in);
                        error00_block1 = _mm256_add_epi32(error00_block1, mul00);
                        error01_block1 = _mm256_add_epi32(error01_block1, mul01);
                        error02_block1 = _mm256_add_epi32(error02_block1, mul02);
                        error03_block1 = _mm256_add_epi32(error03_block1, mul03);
                        error04_block1 = _mm256_add_epi32(error04_block1, mul04);
                        error05_block1 = _mm256_add_epi32(error05_block1, mul05);
                        error06_block1 = _mm256_add_epi32(error06_block1, mul06);
                        error07_block1 = _mm256_add_epi32(error07_block1, mul07);
                        error08_block1 = _mm256_add_epi32(error08_block1, mul08);
                        error09_block1 = _mm256_add_epi32(error09_block1, mul09);
                        error10_block1 = _mm256_add_epi32(error10_block1, mul10);
                        error11_block1 = _mm256_add_epi32(error11_block1, mul11);
                        error12_block1 = _mm256_add_epi32(error12_block1, mul12);
                        error13_block1 = _mm256_add_epi32(error13_block1, mul13);
                        error14_block1 = _mm256_add_epi32(error14_block1, mul14);
                        error15_block1 = _mm256_add_epi32(error15_block1, mul15);
                    }
                }
                // we start with row overlap to ensure equivalence in handling block 3 (to other cases)
                for (int k = overlap; k < blocksize; k++)
                {
                    int m;
                    for (m = 0; m < overlap * 3 - 15; m += 16)
                    {
                        // block3 mulsum part1
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i left_in = loadu_8x16(&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);
                        __m256i mull08 = _mm256_madd_epi16(curr_in_08, left_in);
                        __m256i mull09 = _mm256_madd_epi16(curr_in_09, left_in);
                        __m256i mull10 = _mm256_madd_epi16(curr_in_10, left_in);
                        __m256i mull11 = _mm256_madd_epi16(curr_in_11, left_in);
                        __m256i mull12 = _mm256_madd_epi16(curr_in_12, left_in);
                        __m256i mull13 = _mm256_madd_epi16(curr_in_13, left_in);
                        __m256i mull14 = _mm256_madd_epi16(curr_in_14, left_in);
                        __m256i mull15 = _mm256_madd_epi16(curr_in_15, left_in);
                        error00_block3 = _mm256_add_epi32(mull00, error00_block3);
                        error01_block3 = _mm256_add_epi32(mull01, error01_block3);
                        error02_block3 = _mm256_add_epi32(mull02, error02_block3);
                        error03_block3 = _mm256_add_epi32(mull03, error03_block3);
                        error04_block3 = _mm256_add_epi32(mull04, error04_block3);
                        error05_block3 = _mm256_add_epi32(mull05, error05_block3);
                        error06_block3 = _mm256_add_epi32(mull06, error06_block3);
                        error07_block3 = _mm256_add_epi32(mull07, error07_block3);
                        error08_block3 = _mm256_add_epi32(mull08, error08_block3);
                        error09_block3 = _mm256_add_epi32(mull09, error09_block3);
                        error10_block3 = _mm256_add_epi32(mull10, error10_block3);
                        error11_block3 = _mm256_add_epi32(mull11, error11_block3);
                        error12_block3 = _mm256_add_epi32(mull12, error12_block3);
                        error13_block3 = _mm256_add_epi32(mull13, error13_block3);
                        error14_block3 = _mm256_add_epi32(mull14, error14_block3);
                        error15_block3 = _mm256_add_epi32(mull15, error15_block3);
                    }
                    if (overlap % 16 != 0)
                    {
                        // block3 mulsum part2
                        __m256i curr_in_00 = loadu_8x16(&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i curr_in_01 = loadu_8x16(&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i curr_in_02 = loadu_8x16(&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i curr_in_03 = loadu_8x16(&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i curr_in_04 = loadu_8x16(&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i curr_in_05 = loadu_8x16(&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i curr_in_06 = loadu_8x16(&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i curr_in_07 = loadu_8x16(&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i curr_in_08 = loadu_8x16(&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i curr_in_09 = loadu_8x16(&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i curr_in_10 = loadu_8x16(&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i curr_in_11 = loadu_8x16(&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i curr_in_12 = loadu_8x16(&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i curr_in_13 = loadu_8x16(&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i curr_in_14 = loadu_8x16(&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i curr_in_15 = loadu_8x16(&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i left_in = loadu_8x16(&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);
                        __m256i mull00 = _mm256_madd_epi16(curr_in_00, left_in);
                        __m256i mull01 = _mm256_madd_epi16(curr_in_01, left_in);
                        __m256i mull02 = _mm256_madd_epi16(curr_in_02, left_in);
                        __m256i mull03 = _mm256_madd_epi16(curr_in_03, left_in);
                        __m256i mull04 = _mm256_madd_epi16(curr_in_04, left_in);
                        __m256i mull05 = _mm256_madd_epi16(curr_in_05, left_in);
                        __m256i mull06 = _mm256_madd_epi16(curr_in_06, left_in);
                        __m256i mull07 = _mm256_madd_epi16(curr_in_07, left_in);
                        __m256i mull08 = _mm256_madd_epi16(curr_in_08, left_in);
                        __m256i mull09 = _mm256_madd_epi16(curr_in_09, left_in);
                        __m256i mull10 = _mm256_madd_epi16(curr_in_10, left_in);
                        __m256i mull11 = _mm256_madd_epi16(curr_in_11, left_in);
                        __m256i mull12 = _mm256_madd_epi16(curr_in_12, left_in);
                        __m256i mull13 = _mm256_madd_epi16(curr_in_13, left_in);
                        __m256i mull14 = _mm256_madd_epi16(curr_in_14, left_in);
                        __m256i mull15 = _mm256_madd_epi16(curr_in_15, left_in);
                        error00_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull00, zero, 0b11110000), error00_block3);
                        error01_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull01, zero, 0b11110000), error01_block3);
                        error02_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull02, zero, 0b11110000), error02_block3);
                        error03_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull03, zero, 0b11110000), error03_block3);
                        error04_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull04, zero, 0b11110000), error04_block3);
                        error05_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull05, zero, 0b11110000), error05_block3);
                        error06_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull06, zero, 0b11110000), error06_block3);
                        error07_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull07, zero, 0b11110000), error07_block3);
                        error08_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull08, zero, 0b11110000), error08_block3);
                        error09_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull09, zero, 0b11110000), error09_block3);
                        error10_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull10, zero, 0b11110000), error10_block3);
                        error11_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull11, zero, 0b11110000), error11_block3);
                        error12_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull12, zero, 0b11110000), error12_block3);
                        error13_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull13, zero, 0b11110000), error13_block3);
                        error14_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull14, zero, 0b11110000), error14_block3);
                        error15_block3 = _mm256_add_epi32(_mm256_blend_epi32(mull15, zero, 0b11110000), error15_block3);
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
                error08_block02 = error08_block0;
                error09_block02 = error09_block0;
                error10_block02 = error10_block0;
                error11_block02 = error11_block0;
                error12_block02 = error12_block0;
                error13_block02 = error13_block0;
                error14_block02 = error14_block0;
                error15_block02 = error15_block0;
                error00_block13 = _mm256_slli_epi32(_mm256_add_epi32(error00_block1, error00_block3), 1);
                error01_block13 = _mm256_slli_epi32(_mm256_add_epi32(error01_block1, error01_block3), 1);
                error02_block13 = _mm256_slli_epi32(_mm256_add_epi32(error02_block1, error02_block3), 1);
                error03_block13 = _mm256_slli_epi32(_mm256_add_epi32(error03_block1, error03_block3), 1);
                error04_block13 = _mm256_slli_epi32(_mm256_add_epi32(error04_block1, error04_block3), 1);
                error05_block13 = _mm256_slli_epi32(_mm256_add_epi32(error05_block1, error05_block3), 1);
                error06_block13 = _mm256_slli_epi32(_mm256_add_epi32(error06_block1, error06_block3), 1);
                error07_block13 = _mm256_slli_epi32(_mm256_add_epi32(error07_block1, error07_block3), 1);
                error08_block13 = _mm256_slli_epi32(_mm256_add_epi32(error08_block1, error08_block3), 1);
                error09_block13 = _mm256_slli_epi32(_mm256_add_epi32(error09_block1, error09_block3), 1);
                error10_block13 = _mm256_slli_epi32(_mm256_add_epi32(error10_block1, error10_block3), 1);
                error11_block13 = _mm256_slli_epi32(_mm256_add_epi32(error11_block1, error11_block3), 1);
                error12_block13 = _mm256_slli_epi32(_mm256_add_epi32(error12_block1, error12_block3), 1);
                error13_block13 = _mm256_slli_epi32(_mm256_add_epi32(error13_block1, error13_block3), 1);
                error14_block13 = _mm256_slli_epi32(_mm256_add_epi32(error14_block1, error14_block3), 1);
                error15_block13 = _mm256_slli_epi32(_mm256_add_epi32(error15_block1, error15_block3), 1);
                error00 = _mm256_sub_epi32(error00_block02, error00_block13);
                error01 = _mm256_sub_epi32(error01_block02, error01_block13);
                error02 = _mm256_sub_epi32(error02_block02, error02_block13);
                error03 = _mm256_sub_epi32(error03_block02, error03_block13);
                error04 = _mm256_sub_epi32(error04_block02, error04_block13);
                error05 = _mm256_sub_epi32(error05_block02, error05_block13);
                error06 = _mm256_sub_epi32(error06_block02, error06_block13);
                error07 = _mm256_sub_epi32(error07_block02, error07_block13);
                error08 = _mm256_sub_epi32(error08_block02, error08_block13);
                error09 = _mm256_sub_epi32(error09_block02, error09_block13);
                error10 = _mm256_sub_epi32(error10_block02, error10_block13);
                error11 = _mm256_sub_epi32(error11_block02, error11_block13);
                error12 = _mm256_sub_epi32(error12_block02, error12_block13);
                error13 = _mm256_sub_epi32(error13_block02, error13_block13);
                error14 = _mm256_sub_epi32(error14_block02, error14_block13);
                error15 = _mm256_sub_epi32(error15_block02, error15_block13);
                errors[irow * error_width + icol + 0] = block1_out_integral + block1_in_integral00 + hsum_8x32(error00) + block3_out_integral + block3_in_integral00;
                errors[irow * error_width + icol + 1] = block1_out_integral + block1_in_integral01 + hsum_8x32(error01) + block3_out_integral + block3_in_integral01;
                errors[irow * error_width + icol + 2] = block1_out_integral + block1_in_integral02 + hsum_8x32(error02) + block3_out_integral + block3_in_integral02;
                errors[irow * error_width + icol + 3] = block1_out_integral + block1_in_integral03 + hsum_8x32(error03) + block3_out_integral + block3_in_integral03;
                errors[irow * error_width + icol + 4] = block1_out_integral + block1_in_integral04 + hsum_8x32(error04) + block3_out_integral + block3_in_integral04;
                errors[irow * error_width + icol + 5] = block1_out_integral + block1_in_integral05 + hsum_8x32(error05) + block3_out_integral + block3_in_integral05;
                errors[irow * error_width + icol + 6] = block1_out_integral + block1_in_integral06 + hsum_8x32(error06) + block3_out_integral + block3_in_integral06;
                errors[irow * error_width + icol + 7] = block1_out_integral + block1_in_integral07 + hsum_8x32(error07) + block3_out_integral + block3_in_integral07;
                errors[irow * error_width + icol + 8] = block1_out_integral + block1_in_integral08 + hsum_8x32(error08) + block3_out_integral + block3_in_integral08;
                errors[irow * error_width + icol + 9] = block1_out_integral + block1_in_integral09 + hsum_8x32(error09) + block3_out_integral + block3_in_integral09;
                errors[irow * error_width + icol + 10] = block1_out_integral + block1_in_integral10 + hsum_8x32(error10) + block3_out_integral + block3_in_integral10;
                errors[irow * error_width + icol + 11] = block1_out_integral + block1_in_integral11 + hsum_8x32(error11) + block3_out_integral + block3_in_integral11;
                errors[irow * error_width + icol + 12] = block1_out_integral + block1_in_integral12 + hsum_8x32(error12) + block3_out_integral + block3_in_integral12;
                errors[irow * error_width + icol + 13] = block1_out_integral + block1_in_integral13 + hsum_8x32(error13) + block3_out_integral + block3_in_integral13;
                errors[irow * error_width + icol + 14] = block1_out_integral + block1_in_integral14 + hsum_8x32(error14) + block3_out_integral + block3_in_integral14;
                errors[irow * error_width + icol + 15] = block1_out_integral + block1_in_integral15 + hsum_8x32(error15) + block3_out_integral + block3_in_integral15;
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
