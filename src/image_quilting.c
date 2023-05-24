#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_quilting.h"

#define INTEGRAL(start, height, width, jumpsize) integral[start + height * jumpsize + width] - integral[start + width] - integral[start + height * jumpsize] + integral[start]
// from my earlier answer, with tuning for non-AVX CPUs removed
// static inline
static inline uint32_t hsum_epi32_avx(__m128i x)
{
    __m128i hi64  = _mm_unpackhi_epi64(x, x);           // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
    __m128i sum64 = _mm_add_epi32(hi64, x);
    __m128i hi32  = _mm_shuffle_epi32(sum64, _MM_SHUFFLE(2, 3, 0, 1));    // Swap the low two elements
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32);       // movd
}

// only needs AVX2
static inline uint32_t hsum_8x32(__m256i v)
{
    __m128i sum128 = _mm_add_epi32( 
                 _mm256_castsi256_si128(v),
                 _mm256_extracti128_si256(v, 1)); // silly GCC uses a longer AXV512VL instruction if AVX512 is enabled :/
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
void fill_error_matrix(image_t in_, image_t out_, int orow, int ocol, pixel_t *errors, pixel_t *integral, int lrow, int lcol, int arow, int acol, int blocksize, int overlap, int num_blocks)
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
        pixel_t block3_out_integral = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);

        int icol;
        for (int irow = 0; irow < error_height; irow++)
        {
            for (icol = 0; icol < error_width - 15; icol += 16)
            {
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

                pixel_t block3_in_integral00 = INTEGRAL(block3_in_integral_start00, height3, width3, integral_width);
                pixel_t block3_in_integral01 = INTEGRAL(block3_in_integral_start01, height3, width3, integral_width);
                pixel_t block3_in_integral02 = INTEGRAL(block3_in_integral_start02, height3, width3, integral_width);
                pixel_t block3_in_integral03 = INTEGRAL(block3_in_integral_start03, height3, width3, integral_width);
                pixel_t block3_in_integral04 = INTEGRAL(block3_in_integral_start04, height3, width3, integral_width);
                pixel_t block3_in_integral05 = INTEGRAL(block3_in_integral_start05, height3, width3, integral_width);
                pixel_t block3_in_integral06 = INTEGRAL(block3_in_integral_start06, height3, width3, integral_width);
                pixel_t block3_in_integral07 = INTEGRAL(block3_in_integral_start07, height3, width3, integral_width);
                pixel_t block3_in_integral08 = INTEGRAL(block3_in_integral_start08, height3, width3, integral_width);
                pixel_t block3_in_integral09 = INTEGRAL(block3_in_integral_start09, height3, width3, integral_width);
                pixel_t block3_in_integral10 = INTEGRAL(block3_in_integral_start10, height3, width3, integral_width);
                pixel_t block3_in_integral11 = INTEGRAL(block3_in_integral_start11, height3, width3, integral_width);
                pixel_t block3_in_integral12 = INTEGRAL(block3_in_integral_start12, height3, width3, integral_width);
                pixel_t block3_in_integral13 = INTEGRAL(block3_in_integral_start13, height3, width3, integral_width);
                pixel_t block3_in_integral14 = INTEGRAL(block3_in_integral_start14, height3, width3, integral_width);
                pixel_t block3_in_integral15 = INTEGRAL(block3_in_integral_start15, height3, width3, integral_width);

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

                for (int k = 0; k < blocksize; k++)
                {
                    for (int m = 0; m < overlap * 3; m += 8)
                    {
                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&in[(lrow + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);

                        // mullo until we found a better solution
                        error00 = _mm256_add_epi32(_mm256_mullo_epi32(m00_in, m_out), error00);
                        error01 = _mm256_add_epi32(_mm256_mullo_epi32(m01_in, m_out), error01);
                        error02 = _mm256_add_epi32(_mm256_mullo_epi32(m02_in, m_out), error02);
                        error03 = _mm256_add_epi32(_mm256_mullo_epi32(m03_in, m_out), error03);
                        error04 = _mm256_add_epi32(_mm256_mullo_epi32(m04_in, m_out), error04);
                        error05 = _mm256_add_epi32(_mm256_mullo_epi32(m05_in, m_out), error05);
                        error06 = _mm256_add_epi32(_mm256_mullo_epi32(m06_in, m_out), error06);
                        error07 = _mm256_add_epi32(_mm256_mullo_epi32(m07_in, m_out), error07);
                        error08 = _mm256_add_epi32(_mm256_mullo_epi32(m08_in, m_out), error08);
                        error09 = _mm256_add_epi32(_mm256_mullo_epi32(m09_in, m_out), error09);
                        error10 = _mm256_add_epi32(_mm256_mullo_epi32(m10_in, m_out), error10);
                        error11 = _mm256_add_epi32(_mm256_mullo_epi32(m11_in, m_out), error11);
                        error12 = _mm256_add_epi32(_mm256_mullo_epi32(m12_in, m_out), error12);
                        error13 = _mm256_add_epi32(_mm256_mullo_epi32(m13_in, m_out), error13);
                        error14 = _mm256_add_epi32(_mm256_mullo_epi32(m14_in, m_out), error14);
                        error15 = _mm256_add_epi32(_mm256_mullo_epi32(m15_in, m_out), error15);
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
        pixel_t block1_out_integral = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);

        int icol;
        for (int irow = 0; irow < error_height; irow++)
        {
            for (icol = 0; icol < error_width - 15; icol += 16)
            {
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

                pixel_t block1_in_integral00 = INTEGRAL(block1_in_integral_start00, height1, width1, integral_width);
                pixel_t block1_in_integral01 = INTEGRAL(block1_in_integral_start01, height1, width1, integral_width);
                pixel_t block1_in_integral02 = INTEGRAL(block1_in_integral_start02, height1, width1, integral_width);
                pixel_t block1_in_integral03 = INTEGRAL(block1_in_integral_start03, height1, width1, integral_width);
                pixel_t block1_in_integral04 = INTEGRAL(block1_in_integral_start04, height1, width1, integral_width);
                pixel_t block1_in_integral05 = INTEGRAL(block1_in_integral_start05, height1, width1, integral_width);
                pixel_t block1_in_integral06 = INTEGRAL(block1_in_integral_start06, height1, width1, integral_width);
                pixel_t block1_in_integral07 = INTEGRAL(block1_in_integral_start07, height1, width1, integral_width);
                pixel_t block1_in_integral08 = INTEGRAL(block1_in_integral_start08, height1, width1, integral_width);
                pixel_t block1_in_integral09 = INTEGRAL(block1_in_integral_start09, height1, width1, integral_width);
                pixel_t block1_in_integral10 = INTEGRAL(block1_in_integral_start10, height1, width1, integral_width);
                pixel_t block1_in_integral11 = INTEGRAL(block1_in_integral_start11, height1, width1, integral_width);
                pixel_t block1_in_integral12 = INTEGRAL(block1_in_integral_start12, height1, width1, integral_width);
                pixel_t block1_in_integral13 = INTEGRAL(block1_in_integral_start13, height1, width1, integral_width);
                pixel_t block1_in_integral14 = INTEGRAL(block1_in_integral_start14, height1, width1, integral_width);
                pixel_t block1_in_integral15 = INTEGRAL(block1_in_integral_start15, height1, width1, integral_width);

                __m256i error00_l2 = zero;
                __m256i error01_l2 = zero;
                __m256i error02_l2 = zero;
                __m256i error03_l2 = zero;
                __m256i error04_l2 = zero;
                __m256i error05_l2 = zero;
                __m256i error06_l2 = zero;
                __m256i error07_l2 = zero;
                __m256i error08_l2 = zero;
                __m256i error09_l2 = zero;
                __m256i error10_l2 = zero;
                __m256i error11_l2 = zero;
                __m256i error12_l2 = zero;
                __m256i error13_l2 = zero;
                __m256i error14_l2 = zero;
                __m256i error15_l2 = zero;

                __m256i error00_integral = zero;
                __m256i error01_integral = zero;
                __m256i error02_integral = zero;
                __m256i error03_integral = zero;
                __m256i error04_integral = zero;
                __m256i error05_integral = zero;
                __m256i error06_integral = zero;
                __m256i error07_integral = zero;
                __m256i error08_integral = zero;
                __m256i error09_integral = zero;
                __m256i error10_integral = zero;
                __m256i error11_integral = zero;
                __m256i error12_integral = zero;
                __m256i error13_integral = zero;
                __m256i error14_integral = zero;
                __m256i error15_integral = zero;

                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < overlap * 3; m += 8)
                    {
                        // block2 l2norm
                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol + blocksize - overlap) * 3 + m]);

                        error00_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m00_in, m_out), _mm256_sub_epi32(m00_in, m_out)), error00_l2);
                        error01_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m01_in, m_out), _mm256_sub_epi32(m01_in, m_out)), error01_l2);
                        error02_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m02_in, m_out), _mm256_sub_epi32(m02_in, m_out)), error02_l2);
                        error03_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m03_in, m_out), _mm256_sub_epi32(m03_in, m_out)), error03_l2);
                        error04_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m04_in, m_out), _mm256_sub_epi32(m04_in, m_out)), error04_l2);
                        error05_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m05_in, m_out), _mm256_sub_epi32(m05_in, m_out)), error05_l2);
                        error06_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m06_in, m_out), _mm256_sub_epi32(m06_in, m_out)), error06_l2);
                        error07_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m07_in, m_out), _mm256_sub_epi32(m07_in, m_out)), error07_l2);
                        error08_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m08_in, m_out), _mm256_sub_epi32(m08_in, m_out)), error08_l2);
                        error09_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m09_in, m_out), _mm256_sub_epi32(m09_in, m_out)), error09_l2);
                        error10_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m10_in, m_out), _mm256_sub_epi32(m10_in, m_out)), error10_l2);
                        error11_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m11_in, m_out), _mm256_sub_epi32(m11_in, m_out)), error11_l2);
                        error12_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m12_in, m_out), _mm256_sub_epi32(m12_in, m_out)), error12_l2);
                        error13_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m13_in, m_out), _mm256_sub_epi32(m13_in, m_out)), error13_l2);
                        error14_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m14_in, m_out), _mm256_sub_epi32(m14_in, m_out)), error14_l2);
                        error15_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m15_in, m_out), _mm256_sub_epi32(m15_in, m_out)), error15_l2);
                    }
                    for (int m = 0; m < ((blocksize - overlap) * 3); m += 8)
                    {

                        // block 1 -> mul_sum
                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + acol * 3 + m]);

                        // mullo until we found a better solution
                        error00_integral = _mm256_add_epi32(_mm256_mullo_epi32(m00_in, m_out), error00_integral);
                        error01_integral = _mm256_add_epi32(_mm256_mullo_epi32(m01_in, m_out), error01_integral);
                        error02_integral = _mm256_add_epi32(_mm256_mullo_epi32(m02_in, m_out), error02_integral);
                        error03_integral = _mm256_add_epi32(_mm256_mullo_epi32(m03_in, m_out), error03_integral);
                        error04_integral = _mm256_add_epi32(_mm256_mullo_epi32(m04_in, m_out), error04_integral);
                        error05_integral = _mm256_add_epi32(_mm256_mullo_epi32(m05_in, m_out), error05_integral);
                        error06_integral = _mm256_add_epi32(_mm256_mullo_epi32(m06_in, m_out), error06_integral);
                        error07_integral = _mm256_add_epi32(_mm256_mullo_epi32(m07_in, m_out), error07_integral);
                        error08_integral = _mm256_add_epi32(_mm256_mullo_epi32(m08_in, m_out), error08_integral);
                        error09_integral = _mm256_add_epi32(_mm256_mullo_epi32(m09_in, m_out), error09_integral);
                        error10_integral = _mm256_add_epi32(_mm256_mullo_epi32(m10_in, m_out), error10_integral);
                        error11_integral = _mm256_add_epi32(_mm256_mullo_epi32(m11_in, m_out), error11_integral);
                        error12_integral = _mm256_add_epi32(_mm256_mullo_epi32(m12_in, m_out), error12_integral);
                        error13_integral = _mm256_add_epi32(_mm256_mullo_epi32(m13_in, m_out), error13_integral);
                        error14_integral = _mm256_add_epi32(_mm256_mullo_epi32(m14_in, m_out), error14_integral);
                        error15_integral = _mm256_add_epi32(_mm256_mullo_epi32(m15_in, m_out), error15_integral);
                    }
                }
                errors[irow * error_width + icol + 0] = hsum_8x32(error00_l2) + block1_out_integral - 2 * hsum_8x32(error00_integral) + block1_in_integral00;
                errors[irow * error_width + icol + 1] = hsum_8x32(error01_l2) + block1_out_integral - 2 * hsum_8x32(error01_integral) + block1_in_integral01;
                errors[irow * error_width + icol + 2] = hsum_8x32(error02_l2) + block1_out_integral - 2 * hsum_8x32(error02_integral) + block1_in_integral02;
                errors[irow * error_width + icol + 3] = hsum_8x32(error03_l2) + block1_out_integral - 2 * hsum_8x32(error03_integral) + block1_in_integral03;
                errors[irow * error_width + icol + 4] = hsum_8x32(error04_l2) + block1_out_integral - 2 * hsum_8x32(error04_integral) + block1_in_integral04;
                errors[irow * error_width + icol + 5] = hsum_8x32(error05_l2) + block1_out_integral - 2 * hsum_8x32(error05_integral) + block1_in_integral05;
                errors[irow * error_width + icol + 6] = hsum_8x32(error06_l2) + block1_out_integral - 2 * hsum_8x32(error06_integral) + block1_in_integral06;
                errors[irow * error_width + icol + 7] = hsum_8x32(error07_l2) + block1_out_integral - 2 * hsum_8x32(error07_integral) + block1_in_integral07;
                errors[irow * error_width + icol + 8] = hsum_8x32(error08_l2) + block1_out_integral - 2 * hsum_8x32(error08_integral) + block1_in_integral08;
                errors[irow * error_width + icol + 9] = hsum_8x32(error09_l2) + block1_out_integral - 2 * hsum_8x32(error09_integral) + block1_in_integral09;
                errors[irow * error_width + icol + 10] = hsum_8x32(error10_l2) + block1_out_integral - 2 * hsum_8x32(error10_integral) + block1_in_integral10;
                errors[irow * error_width + icol + 11] = hsum_8x32(error11_l2) + block1_out_integral - 2 * hsum_8x32(error11_integral) + block1_in_integral11;
                errors[irow * error_width + icol + 12] = hsum_8x32(error12_l2) + block1_out_integral - 2 * hsum_8x32(error12_integral) + block1_in_integral12;
                errors[irow * error_width + icol + 13] = hsum_8x32(error13_l2) + block1_out_integral - 2 * hsum_8x32(error13_integral) + block1_in_integral13;
                errors[irow * error_width + icol + 14] = hsum_8x32(error14_l2) + block1_out_integral - 2 * hsum_8x32(error14_integral) + block1_in_integral14;
                errors[irow * error_width + icol + 15] = hsum_8x32(error15_l2) + block1_out_integral - 2 * hsum_8x32(error15_integral) + block1_in_integral15;
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
        pixel_t block1_out_integral = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        // precalculation for block 3
        int block3_out_integral_start = (lrow + overlap) * integral_width + lcol + (blocksize - overlap);
        pixel_t block3_out_integral = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);

        for (int irow = 0; irow < error_height; irow++)
        {
            for (int icol = 0; icol < error_width - 15; icol += 16)
            {

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

                pixel_t block1_in_integral00 = INTEGRAL(block1_in_integral_start00, height1, width1, integral_width);
                pixel_t block1_in_integral01 = INTEGRAL(block1_in_integral_start01, height1, width1, integral_width);
                pixel_t block1_in_integral02 = INTEGRAL(block1_in_integral_start02, height1, width1, integral_width);
                pixel_t block1_in_integral03 = INTEGRAL(block1_in_integral_start03, height1, width1, integral_width);
                pixel_t block1_in_integral04 = INTEGRAL(block1_in_integral_start04, height1, width1, integral_width);
                pixel_t block1_in_integral05 = INTEGRAL(block1_in_integral_start05, height1, width1, integral_width);
                pixel_t block1_in_integral06 = INTEGRAL(block1_in_integral_start06, height1, width1, integral_width);
                pixel_t block1_in_integral07 = INTEGRAL(block1_in_integral_start07, height1, width1, integral_width);
                pixel_t block1_in_integral08 = INTEGRAL(block1_in_integral_start08, height1, width1, integral_width);
                pixel_t block1_in_integral09 = INTEGRAL(block1_in_integral_start09, height1, width1, integral_width);
                pixel_t block1_in_integral10 = INTEGRAL(block1_in_integral_start10, height1, width1, integral_width);
                pixel_t block1_in_integral11 = INTEGRAL(block1_in_integral_start11, height1, width1, integral_width);
                pixel_t block1_in_integral12 = INTEGRAL(block1_in_integral_start12, height1, width1, integral_width);
                pixel_t block1_in_integral13 = INTEGRAL(block1_in_integral_start13, height1, width1, integral_width);
                pixel_t block1_in_integral14 = INTEGRAL(block1_in_integral_start14, height1, width1, integral_width);
                pixel_t block1_in_integral15 = INTEGRAL(block1_in_integral_start15, height1, width1, integral_width);

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

                pixel_t block3_in_integral00 = INTEGRAL(block3_in_integral_start00, height3, width3, integral_width);
                pixel_t block3_in_integral01 = INTEGRAL(block3_in_integral_start01, height3, width3, integral_width);
                pixel_t block3_in_integral02 = INTEGRAL(block3_in_integral_start02, height3, width3, integral_width);
                pixel_t block3_in_integral03 = INTEGRAL(block3_in_integral_start03, height3, width3, integral_width);
                pixel_t block3_in_integral04 = INTEGRAL(block3_in_integral_start04, height3, width3, integral_width);
                pixel_t block3_in_integral05 = INTEGRAL(block3_in_integral_start05, height3, width3, integral_width);
                pixel_t block3_in_integral06 = INTEGRAL(block3_in_integral_start06, height3, width3, integral_width);
                pixel_t block3_in_integral07 = INTEGRAL(block3_in_integral_start07, height3, width3, integral_width);
                pixel_t block3_in_integral08 = INTEGRAL(block3_in_integral_start08, height3, width3, integral_width);
                pixel_t block3_in_integral09 = INTEGRAL(block3_in_integral_start09, height3, width3, integral_width);
                pixel_t block3_in_integral10 = INTEGRAL(block3_in_integral_start10, height3, width3, integral_width);
                pixel_t block3_in_integral11 = INTEGRAL(block3_in_integral_start11, height3, width3, integral_width);
                pixel_t block3_in_integral12 = INTEGRAL(block3_in_integral_start12, height3, width3, integral_width);
                pixel_t block3_in_integral13 = INTEGRAL(block3_in_integral_start13, height3, width3, integral_width);
                pixel_t block3_in_integral14 = INTEGRAL(block3_in_integral_start14, height3, width3, integral_width);
                pixel_t block3_in_integral15 = INTEGRAL(block3_in_integral_start15, height3, width3, integral_width);

                __m256i error00_l2 = zero;
                __m256i error01_l2 = zero;
                __m256i error02_l2 = zero;
                __m256i error03_l2 = zero;
                __m256i error04_l2 = zero;
                __m256i error05_l2 = zero;
                __m256i error06_l2 = zero;
                __m256i error07_l2 = zero;
                __m256i error08_l2 = zero;
                __m256i error09_l2 = zero;
                __m256i error10_l2 = zero;
                __m256i error11_l2 = zero;
                __m256i error12_l2 = zero;
                __m256i error13_l2 = zero;
                __m256i error14_l2 = zero;
                __m256i error15_l2 = zero;

                __m256i error00_integral = zero;
                __m256i error01_integral = zero;
                __m256i error02_integral = zero;
                __m256i error03_integral = zero;
                __m256i error04_integral = zero;
                __m256i error05_integral = zero;
                __m256i error06_integral = zero;
                __m256i error07_integral = zero;
                __m256i error08_integral = zero;
                __m256i error09_integral = zero;
                __m256i error10_integral = zero;
                __m256i error11_integral = zero;
                __m256i error12_integral = zero;
                __m256i error13_integral = zero;
                __m256i error14_integral = zero;
                __m256i error15_integral = zero;

                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < overlap * 3; m += 8)
                    {

                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + ocol * 3 + m]);

                        error00_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m00_in, m_out), _mm256_sub_epi32(m00_in, m_out)), error00_l2);
                        error01_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m01_in, m_out), _mm256_sub_epi32(m01_in, m_out)), error01_l2);
                        error02_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m02_in, m_out), _mm256_sub_epi32(m02_in, m_out)), error02_l2);
                        error03_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m03_in, m_out), _mm256_sub_epi32(m03_in, m_out)), error03_l2);
                        error04_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m04_in, m_out), _mm256_sub_epi32(m04_in, m_out)), error04_l2);
                        error05_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m05_in, m_out), _mm256_sub_epi32(m05_in, m_out)), error05_l2);
                        error06_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m06_in, m_out), _mm256_sub_epi32(m06_in, m_out)), error06_l2);
                        error07_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m07_in, m_out), _mm256_sub_epi32(m07_in, m_out)), error07_l2);
                        error08_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m08_in, m_out), _mm256_sub_epi32(m08_in, m_out)), error08_l2);
                        error09_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m09_in, m_out), _mm256_sub_epi32(m09_in, m_out)), error09_l2);
                        error10_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m10_in, m_out), _mm256_sub_epi32(m10_in, m_out)), error10_l2);
                        error11_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m11_in, m_out), _mm256_sub_epi32(m11_in, m_out)), error11_l2);
                        error12_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m12_in, m_out), _mm256_sub_epi32(m12_in, m_out)), error12_l2);
                        error13_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m13_in, m_out), _mm256_sub_epi32(m13_in, m_out)), error13_l2);
                        error14_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m14_in, m_out), _mm256_sub_epi32(m14_in, m_out)), error14_l2);
                        error15_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m15_in, m_out), _mm256_sub_epi32(m15_in, m_out)), error15_l2);

                        //------------------
                        __m256i m00_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 0) * 3 + m]);
                        __m256i m01_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 1) * 3 + m]);
                        __m256i m02_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 2) * 3 + m]);
                        __m256i m03_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 3) * 3 + m]);
                        __m256i m04_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 4) * 3 + m]);
                        __m256i m05_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 5) * 3 + m]);
                        __m256i m06_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 6) * 3 + m]);
                        __m256i m07_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 7) * 3 + m]);
                        __m256i m08_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 8) * 3 + m]);
                        __m256i m09_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 9) * 3 + m]);
                        __m256i m10_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 10) * 3 + m]);
                        __m256i m11_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 11) * 3 + m]);
                        __m256i m12_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 12) * 3 + m]);
                        __m256i m13_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 13) * 3 + m]);
                        __m256i m14_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 14) * 3 + m]);
                        __m256i m15_in_2 = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + blocksize - overlap + 15) * 3 + m]);
                        __m256i m_out_2 = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + (ocol + blocksize - overlap) * 3 + m]);

                        error00_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m00_in_2, m_out_2), _mm256_sub_epi32(m00_in_2, m_out_2)), error00_l2);
                        error01_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m01_in_2, m_out_2), _mm256_sub_epi32(m01_in_2, m_out_2)), error01_l2);
                        error02_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m02_in_2, m_out_2), _mm256_sub_epi32(m02_in_2, m_out_2)), error02_l2);
                        error03_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m03_in_2, m_out_2), _mm256_sub_epi32(m03_in_2, m_out_2)), error03_l2);
                        error04_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m04_in_2, m_out_2), _mm256_sub_epi32(m04_in_2, m_out_2)), error04_l2);
                        error05_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m05_in_2, m_out_2), _mm256_sub_epi32(m05_in_2, m_out_2)), error05_l2);
                        error06_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m06_in_2, m_out_2), _mm256_sub_epi32(m06_in_2, m_out_2)), error06_l2);
                        error07_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m07_in_2, m_out_2), _mm256_sub_epi32(m07_in_2, m_out_2)), error07_l2);
                        error08_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m08_in_2, m_out_2), _mm256_sub_epi32(m08_in_2, m_out_2)), error08_l2);
                        error09_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m09_in_2, m_out_2), _mm256_sub_epi32(m09_in_2, m_out_2)), error09_l2);
                        error10_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m10_in_2, m_out_2), _mm256_sub_epi32(m10_in_2, m_out_2)), error10_l2);
                        error11_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m11_in_2, m_out_2), _mm256_sub_epi32(m11_in_2, m_out_2)), error11_l2);
                        error12_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m12_in_2, m_out_2), _mm256_sub_epi32(m12_in_2, m_out_2)), error12_l2);
                        error13_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m13_in_2, m_out_2), _mm256_sub_epi32(m13_in_2, m_out_2)), error13_l2);
                        error14_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m14_in_2, m_out_2), _mm256_sub_epi32(m14_in_2, m_out_2)), error14_l2);
                        error15_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m15_in_2, m_out_2), _mm256_sub_epi32(m15_in_2, m_out_2)), error15_l2);
                    }
                    for (int m = 0; m < (blocksize - 2 * overlap) * 3; m += 8)
                    {

                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol + overlap) * 3 + m]);

                        error00_integral = _mm256_add_epi32(_mm256_mullo_epi32(m00_in, m_out), error00_integral);
                        error01_integral = _mm256_add_epi32(_mm256_mullo_epi32(m01_in, m_out), error01_integral);
                        error02_integral = _mm256_add_epi32(_mm256_mullo_epi32(m02_in, m_out), error02_integral);
                        error03_integral = _mm256_add_epi32(_mm256_mullo_epi32(m03_in, m_out), error03_integral);
                        error04_integral = _mm256_add_epi32(_mm256_mullo_epi32(m04_in, m_out), error04_integral);
                        error05_integral = _mm256_add_epi32(_mm256_mullo_epi32(m05_in, m_out), error05_integral);
                        error06_integral = _mm256_add_epi32(_mm256_mullo_epi32(m06_in, m_out), error06_integral);
                        error07_integral = _mm256_add_epi32(_mm256_mullo_epi32(m07_in, m_out), error07_integral);
                        error08_integral = _mm256_add_epi32(_mm256_mullo_epi32(m08_in, m_out), error08_integral);
                        error09_integral = _mm256_add_epi32(_mm256_mullo_epi32(m09_in, m_out), error09_integral);
                        error10_integral = _mm256_add_epi32(_mm256_mullo_epi32(m10_in, m_out), error10_integral);
                        error11_integral = _mm256_add_epi32(_mm256_mullo_epi32(m11_in, m_out), error11_integral);
                        error12_integral = _mm256_add_epi32(_mm256_mullo_epi32(m12_in, m_out), error12_integral);
                        error13_integral = _mm256_add_epi32(_mm256_mullo_epi32(m13_in, m_out), error13_integral);
                        error14_integral = _mm256_add_epi32(_mm256_mullo_epi32(m14_in, m_out), error14_integral);
                        error15_integral = _mm256_add_epi32(_mm256_mullo_epi32(m15_in, m_out), error15_integral);
                    }
                }
                for (int k = 0; k < blocksize - overlap; k++)
                {
                    for (int m = 0; m < overlap * 3; m += 8)
                    {

                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);

                        error00_integral = _mm256_add_epi32(_mm256_mullo_epi32(m00_in, m_out), error00_integral);
                        error01_integral = _mm256_add_epi32(_mm256_mullo_epi32(m01_in, m_out), error01_integral);
                        error02_integral = _mm256_add_epi32(_mm256_mullo_epi32(m02_in, m_out), error02_integral);
                        error03_integral = _mm256_add_epi32(_mm256_mullo_epi32(m03_in, m_out), error03_integral);
                        error04_integral = _mm256_add_epi32(_mm256_mullo_epi32(m04_in, m_out), error04_integral);
                        error05_integral = _mm256_add_epi32(_mm256_mullo_epi32(m05_in, m_out), error05_integral);
                        error06_integral = _mm256_add_epi32(_mm256_mullo_epi32(m06_in, m_out), error06_integral);
                        error07_integral = _mm256_add_epi32(_mm256_mullo_epi32(m07_in, m_out), error07_integral);
                        error08_integral = _mm256_add_epi32(_mm256_mullo_epi32(m08_in, m_out), error08_integral);
                        error09_integral = _mm256_add_epi32(_mm256_mullo_epi32(m09_in, m_out), error09_integral);
                        error10_integral = _mm256_add_epi32(_mm256_mullo_epi32(m10_in, m_out), error10_integral);
                        error11_integral = _mm256_add_epi32(_mm256_mullo_epi32(m11_in, m_out), error11_integral);
                        error12_integral = _mm256_add_epi32(_mm256_mullo_epi32(m12_in, m_out), error12_integral);
                        error13_integral = _mm256_add_epi32(_mm256_mullo_epi32(m13_in, m_out), error13_integral);
                        error14_integral = _mm256_add_epi32(_mm256_mullo_epi32(m14_in, m_out), error14_integral);
                        error15_integral = _mm256_add_epi32(_mm256_mullo_epi32(m15_in, m_out), error15_integral);
                    }
                }
                errors[irow * error_width + icol + 0] = hsum_8x32(error00_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error00_integral) + block1_in_integral00 + block3_in_integral00;
                errors[irow * error_width + icol + 1] = hsum_8x32(error01_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error01_integral) + block1_in_integral01 + block3_in_integral01;
                errors[irow * error_width + icol + 2] = hsum_8x32(error02_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error02_integral) + block1_in_integral02 + block3_in_integral02;
                errors[irow * error_width + icol + 3] = hsum_8x32(error03_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error03_integral) + block1_in_integral03 + block3_in_integral03;
                errors[irow * error_width + icol + 4] = hsum_8x32(error04_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error04_integral) + block1_in_integral04 + block3_in_integral04;
                errors[irow * error_width + icol + 5] = hsum_8x32(error05_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error05_integral) + block1_in_integral05 + block3_in_integral05;
                errors[irow * error_width + icol + 6] = hsum_8x32(error06_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error06_integral) + block1_in_integral06 + block3_in_integral06;
                errors[irow * error_width + icol + 7] = hsum_8x32(error07_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error07_integral) + block1_in_integral07 + block3_in_integral07;
                errors[irow * error_width + icol + 8] = hsum_8x32(error08_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error08_integral) + block1_in_integral08 + block3_in_integral08;
                errors[irow * error_width + icol + 9] = hsum_8x32(error09_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error09_integral) + block1_in_integral09 + block3_in_integral09;
                errors[irow * error_width + icol + 10] = hsum_8x32(error10_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error10_integral) + block1_in_integral10 + block3_in_integral10;
                errors[irow * error_width + icol + 11] = hsum_8x32(error11_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error11_integral) + block1_in_integral11 + block3_in_integral11;
                errors[irow * error_width + icol + 12] = hsum_8x32(error12_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error12_integral) + block1_in_integral12 + block3_in_integral12;
                errors[irow * error_width + icol + 13] = hsum_8x32(error13_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error13_integral) + block1_in_integral13 + block3_in_integral13;
                errors[irow * error_width + icol + 14] = hsum_8x32(error14_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error14_integral) + block1_in_integral14 + block3_in_integral14;
                errors[irow * error_width + icol + 15] = hsum_8x32(error15_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error15_integral) + block1_in_integral15 + block3_in_integral15;
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
        pixel_t block1_out_integral = INTEGRAL(block1_out_integral_start, height1, width1, integral_width);
        // precalculation for block 3
        int block3_out_integral_start = (lrow + overlap) * integral_width + lcol + (blocksize - overlap);
        pixel_t block3_out_integral = INTEGRAL(block3_out_integral_start, height3, width3, integral_width);

        int icol;
        for (int irow = 0; irow < error_height; irow++)
        {
            for (icol = 0; icol < error_width - 15; icol += 16)
            {
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

                pixel_t block1_in_integral00 = INTEGRAL(block1_in_integral_start00, height1, width1, integral_width);
                pixel_t block1_in_integral01 = INTEGRAL(block1_in_integral_start01, height1, width1, integral_width);
                pixel_t block1_in_integral02 = INTEGRAL(block1_in_integral_start02, height1, width1, integral_width);
                pixel_t block1_in_integral03 = INTEGRAL(block1_in_integral_start03, height1, width1, integral_width);
                pixel_t block1_in_integral04 = INTEGRAL(block1_in_integral_start04, height1, width1, integral_width);
                pixel_t block1_in_integral05 = INTEGRAL(block1_in_integral_start05, height1, width1, integral_width);
                pixel_t block1_in_integral06 = INTEGRAL(block1_in_integral_start06, height1, width1, integral_width);
                pixel_t block1_in_integral07 = INTEGRAL(block1_in_integral_start07, height1, width1, integral_width);
                pixel_t block1_in_integral08 = INTEGRAL(block1_in_integral_start08, height1, width1, integral_width);
                pixel_t block1_in_integral09 = INTEGRAL(block1_in_integral_start09, height1, width1, integral_width);
                pixel_t block1_in_integral10 = INTEGRAL(block1_in_integral_start10, height1, width1, integral_width);
                pixel_t block1_in_integral11 = INTEGRAL(block1_in_integral_start11, height1, width1, integral_width);
                pixel_t block1_in_integral12 = INTEGRAL(block1_in_integral_start12, height1, width1, integral_width);
                pixel_t block1_in_integral13 = INTEGRAL(block1_in_integral_start13, height1, width1, integral_width);
                pixel_t block1_in_integral14 = INTEGRAL(block1_in_integral_start14, height1, width1, integral_width);
                pixel_t block1_in_integral15 = INTEGRAL(block1_in_integral_start15, height1, width1, integral_width);

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

                pixel_t block3_in_integral00 = INTEGRAL(block3_in_integral_start00, height3, width3, integral_width);
                pixel_t block3_in_integral01 = INTEGRAL(block3_in_integral_start01, height3, width3, integral_width);
                pixel_t block3_in_integral02 = INTEGRAL(block3_in_integral_start02, height3, width3, integral_width);
                pixel_t block3_in_integral03 = INTEGRAL(block3_in_integral_start03, height3, width3, integral_width);
                pixel_t block3_in_integral04 = INTEGRAL(block3_in_integral_start04, height3, width3, integral_width);
                pixel_t block3_in_integral05 = INTEGRAL(block3_in_integral_start05, height3, width3, integral_width);
                pixel_t block3_in_integral06 = INTEGRAL(block3_in_integral_start06, height3, width3, integral_width);
                pixel_t block3_in_integral07 = INTEGRAL(block3_in_integral_start07, height3, width3, integral_width);
                pixel_t block3_in_integral08 = INTEGRAL(block3_in_integral_start08, height3, width3, integral_width);
                pixel_t block3_in_integral09 = INTEGRAL(block3_in_integral_start09, height3, width3, integral_width);
                pixel_t block3_in_integral10 = INTEGRAL(block3_in_integral_start10, height3, width3, integral_width);
                pixel_t block3_in_integral11 = INTEGRAL(block3_in_integral_start11, height3, width3, integral_width);
                pixel_t block3_in_integral12 = INTEGRAL(block3_in_integral_start12, height3, width3, integral_width);
                pixel_t block3_in_integral13 = INTEGRAL(block3_in_integral_start13, height3, width3, integral_width);
                pixel_t block3_in_integral14 = INTEGRAL(block3_in_integral_start14, height3, width3, integral_width);
                pixel_t block3_in_integral15 = INTEGRAL(block3_in_integral_start15, height3, width3, integral_width);

                __m256i error00_l2 = zero;
                __m256i error01_l2 = zero;
                __m256i error02_l2 = zero;
                __m256i error03_l2 = zero;
                __m256i error04_l2 = zero;
                __m256i error05_l2 = zero;
                __m256i error06_l2 = zero;
                __m256i error07_l2 = zero;
                __m256i error08_l2 = zero;
                __m256i error09_l2 = zero;
                __m256i error10_l2 = zero;
                __m256i error11_l2 = zero;
                __m256i error12_l2 = zero;
                __m256i error13_l2 = zero;
                __m256i error14_l2 = zero;
                __m256i error15_l2 = zero;

                __m256i error00_integral = zero;
                __m256i error01_integral = zero;
                __m256i error02_integral = zero;
                __m256i error03_integral = zero;
                __m256i error04_integral = zero;
                __m256i error05_integral = zero;
                __m256i error06_integral = zero;
                __m256i error07_integral = zero;
                __m256i error08_integral = zero;
                __m256i error09_integral = zero;
                __m256i error10_integral = zero;
                __m256i error11_integral = zero;
                __m256i error12_integral = zero;
                __m256i error13_integral = zero;
                __m256i error14_integral = zero;
                __m256i error15_integral = zero;

                for (int k = 0; k < overlap; k++)
                {
                    for (int m = 0; m < overlap * 3; m += 8)
                    {

                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&out[(orow + k) * ojump + ocol * 3 + m]);

                        error00_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m00_in, m_out), _mm256_sub_epi32(m00_in, m_out)), error00_l2);
                        error01_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m01_in, m_out), _mm256_sub_epi32(m01_in, m_out)), error01_l2);
                        error02_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m02_in, m_out), _mm256_sub_epi32(m02_in, m_out)), error02_l2);
                        error03_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m03_in, m_out), _mm256_sub_epi32(m03_in, m_out)), error03_l2);
                        error04_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m04_in, m_out), _mm256_sub_epi32(m04_in, m_out)), error04_l2);
                        error05_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m05_in, m_out), _mm256_sub_epi32(m05_in, m_out)), error05_l2);
                        error06_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m06_in, m_out), _mm256_sub_epi32(m06_in, m_out)), error06_l2);
                        error07_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m07_in, m_out), _mm256_sub_epi32(m07_in, m_out)), error07_l2);
                        error08_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m08_in, m_out), _mm256_sub_epi32(m08_in, m_out)), error08_l2);
                        error09_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m09_in, m_out), _mm256_sub_epi32(m09_in, m_out)), error09_l2);
                        error10_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m10_in, m_out), _mm256_sub_epi32(m10_in, m_out)), error10_l2);
                        error11_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m11_in, m_out), _mm256_sub_epi32(m11_in, m_out)), error11_l2);
                        error12_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m12_in, m_out), _mm256_sub_epi32(m12_in, m_out)), error12_l2);
                        error13_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m13_in, m_out), _mm256_sub_epi32(m13_in, m_out)), error13_l2);
                        error14_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m14_in, m_out), _mm256_sub_epi32(m14_in, m_out)), error14_l2);
                        error15_l2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_sub_epi32(m15_in, m_out), _mm256_sub_epi32(m15_in, m_out)), error15_l2);
                    }
                    for (int m = 0; m < (blocksize - overlap) * 3; m += 8)
                    {

                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k) * ijump + (icol + overlap + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&in[(arow + (blocksize - overlap) + k) * ijump + (acol + overlap) * 3 + m]);

                        error00_integral = _mm256_add_epi32(_mm256_mullo_epi32(m00_in, m_out), error00_integral);
                        error01_integral = _mm256_add_epi32(_mm256_mullo_epi32(m01_in, m_out), error01_integral);
                        error02_integral = _mm256_add_epi32(_mm256_mullo_epi32(m02_in, m_out), error02_integral);
                        error03_integral = _mm256_add_epi32(_mm256_mullo_epi32(m03_in, m_out), error03_integral);
                        error04_integral = _mm256_add_epi32(_mm256_mullo_epi32(m04_in, m_out), error04_integral);
                        error05_integral = _mm256_add_epi32(_mm256_mullo_epi32(m05_in, m_out), error05_integral);
                        error06_integral = _mm256_add_epi32(_mm256_mullo_epi32(m06_in, m_out), error06_integral);
                        error07_integral = _mm256_add_epi32(_mm256_mullo_epi32(m07_in, m_out), error07_integral);
                        error08_integral = _mm256_add_epi32(_mm256_mullo_epi32(m08_in, m_out), error08_integral);
                        error09_integral = _mm256_add_epi32(_mm256_mullo_epi32(m09_in, m_out), error09_integral);
                        error10_integral = _mm256_add_epi32(_mm256_mullo_epi32(m10_in, m_out), error10_integral);
                        error11_integral = _mm256_add_epi32(_mm256_mullo_epi32(m11_in, m_out), error11_integral);
                        error12_integral = _mm256_add_epi32(_mm256_mullo_epi32(m12_in, m_out), error12_integral);
                        error13_integral = _mm256_add_epi32(_mm256_mullo_epi32(m13_in, m_out), error13_integral);
                        error14_integral = _mm256_add_epi32(_mm256_mullo_epi32(m14_in, m_out), error14_integral);
                        error15_integral = _mm256_add_epi32(_mm256_mullo_epi32(m15_in, m_out), error15_integral);
                    }
                }
                for (int k = 0; k < blocksize - overlap; k++)
                {
                    for (int m = 0; m < overlap * 3; m += 8)
                    {
                        __m256i m00_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 0) * 3 + m]);
                        __m256i m01_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 1) * 3 + m]);
                        __m256i m02_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 2) * 3 + m]);
                        __m256i m03_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 3) * 3 + m]);
                        __m256i m04_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 4) * 3 + m]);
                        __m256i m05_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 5) * 3 + m]);
                        __m256i m06_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 6) * 3 + m]);
                        __m256i m07_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 7) * 3 + m]);
                        __m256i m08_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 8) * 3 + m]);
                        __m256i m09_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 9) * 3 + m]);
                        __m256i m10_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 10) * 3 + m]);
                        __m256i m11_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 11) * 3 + m]);
                        __m256i m12_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 12) * 3 + m]);
                        __m256i m13_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 13) * 3 + m]);
                        __m256i m14_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 14) * 3 + m]);
                        __m256i m15_in = _mm256_loadu_si256((__m256i_u *)&in[(irow + k + overlap) * ijump + (icol + 15) * 3 + m]);
                        __m256i m_out = _mm256_loadu_si256((__m256i_u *)&in[(lrow + overlap + k) * ijump + (lcol + blocksize - overlap) * 3 + m]);

                        error00_integral = _mm256_add_epi32(_mm256_mullo_epi32(m00_in, m_out), error00_integral);
                        error01_integral = _mm256_add_epi32(_mm256_mullo_epi32(m01_in, m_out), error01_integral);
                        error02_integral = _mm256_add_epi32(_mm256_mullo_epi32(m02_in, m_out), error02_integral);
                        error03_integral = _mm256_add_epi32(_mm256_mullo_epi32(m03_in, m_out), error03_integral);
                        error04_integral = _mm256_add_epi32(_mm256_mullo_epi32(m04_in, m_out), error04_integral);
                        error05_integral = _mm256_add_epi32(_mm256_mullo_epi32(m05_in, m_out), error05_integral);
                        error06_integral = _mm256_add_epi32(_mm256_mullo_epi32(m06_in, m_out), error06_integral);
                        error07_integral = _mm256_add_epi32(_mm256_mullo_epi32(m07_in, m_out), error07_integral);
                        error08_integral = _mm256_add_epi32(_mm256_mullo_epi32(m08_in, m_out), error08_integral);
                        error09_integral = _mm256_add_epi32(_mm256_mullo_epi32(m09_in, m_out), error09_integral);
                        error10_integral = _mm256_add_epi32(_mm256_mullo_epi32(m10_in, m_out), error10_integral);
                        error11_integral = _mm256_add_epi32(_mm256_mullo_epi32(m11_in, m_out), error11_integral);
                        error12_integral = _mm256_add_epi32(_mm256_mullo_epi32(m12_in, m_out), error12_integral);
                        error13_integral = _mm256_add_epi32(_mm256_mullo_epi32(m13_in, m_out), error13_integral);
                        error14_integral = _mm256_add_epi32(_mm256_mullo_epi32(m14_in, m_out), error14_integral);
                        error15_integral = _mm256_add_epi32(_mm256_mullo_epi32(m15_in, m_out), error15_integral);
                    }
                }
                errors[irow * error_width + icol + 0] = hsum_8x32(error00_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error00_integral) + block1_in_integral00 + block3_in_integral00;
                errors[irow * error_width + icol + 1] = hsum_8x32(error01_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error01_integral) + block1_in_integral01 + block3_in_integral01;
                errors[irow * error_width + icol + 2] = hsum_8x32(error02_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error02_integral) + block1_in_integral02 + block3_in_integral02;
                errors[irow * error_width + icol + 3] = hsum_8x32(error03_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error03_integral) + block1_in_integral03 + block3_in_integral03;
                errors[irow * error_width + icol + 4] = hsum_8x32(error04_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error04_integral) + block1_in_integral04 + block3_in_integral04;
                errors[irow * error_width + icol + 5] = hsum_8x32(error05_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error05_integral) + block1_in_integral05 + block3_in_integral05;
                errors[irow * error_width + icol + 6] = hsum_8x32(error06_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error06_integral) + block1_in_integral06 + block3_in_integral06;
                errors[irow * error_width + icol + 7] = hsum_8x32(error07_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error07_integral) + block1_in_integral07 + block3_in_integral07;
                errors[irow * error_width + icol + 8] = hsum_8x32(error08_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error08_integral) + block1_in_integral08 + block3_in_integral08;
                errors[irow * error_width + icol + 9] = hsum_8x32(error09_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error09_integral) + block1_in_integral09 + block3_in_integral09;
                errors[irow * error_width + icol + 10] = hsum_8x32(error10_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error10_integral) + block1_in_integral10 + block3_in_integral10;
                errors[irow * error_width + icol + 11] = hsum_8x32(error11_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error11_integral) + block1_in_integral11 + block3_in_integral11;
                errors[irow * error_width + icol + 12] = hsum_8x32(error12_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error12_integral) + block1_in_integral12 + block3_in_integral12;
                errors[irow * error_width + icol + 13] = hsum_8x32(error13_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error13_integral) + block1_in_integral13 + block3_in_integral13;
                errors[irow * error_width + icol + 14] = hsum_8x32(error14_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error14_integral) + block1_in_integral14 + block3_in_integral14;
                errors[irow * error_width + icol + 15] = hsum_8x32(error15_l2) + block1_out_integral + block3_out_integral - 2 * hsum_8x32(error15_integral) + block1_in_integral15 + block3_in_integral15;
            }
        }
    }
}

/* Here the integral image of the input image is calculated, where the left and upper border consists of zeros to facilitate the calculation of the integral */
void generate_integral_image(image_t in, pixel_t *out)
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
            pixel_t temp = out[i * outjump + j - 1] + out[(i - 1) * outjump + j] - out[(i - 1) * outjump + j - 1];
            // offset -1 in row and col as we start with (1, 1)
            temp += in.data[(i - 1) * injump + 3 * j - 3] * in.data[(i - 1) * injump + 3 * j - 3];
            temp += in.data[(i - 1) * injump + 3 * j - 2] * in.data[(i - 1) * injump + 3 * j - 2];
            temp += in.data[(i - 1) * injump + 3 * j - 1] * in.data[(i - 1) * injump + 3 * j - 1];
            out[i * outjump + j] = temp;
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
    pixel_t tol_range = min_error + (min_error * tol_nom) / tol_den;
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
    pixel_t errors[errorlen];

    // Integral image pads with 0 by 1 row and 1 col to simplify calculations
    pixel_t integral[(in.height + 1) * (in.width + 1)];

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
