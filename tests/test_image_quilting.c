#include <stdio.h>
#include <stdlib.h>

#include "acutest.h"
#include "image_quilting.h"

/* preparation of data for generate_integral_image */
static pixel_t data_gen_integral[] = {
    1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
    5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8,
    9, 9, 9, 10, 10, 10, 11, 11, 11, 255, 255, 255};
static const image_t in_gen_integral = {data_gen_integral, 4, 3, 3};
static const error_t out_gen_integral[] = {
    0, 0, 0, 0, 0,
    0, 3, 15, 42, 90,
    0, 78, 198, 372, 612,
    0, 321, 741, 1278, 196593U};

/* preparation of data for fill_error_matrix (should work with vector assumuption)*/
#define BLOCKSIZE 32
#define OVERLAP 8
#define NUM_BLOCKS 3 /* -> thus we get a 3x3 block image */
#define INSIZE 200
#define OUTSIZE (NUM_BLOCKS * BLOCKSIZE - (NUM_BLOCKS - 1) * OVERLAP)
#define ERRSIZE (INSIZE - BLOCKSIZE + 1)
static pixel_t in[INSIZE * INSIZE * 3]; /* 3 because of the colors */
static pixel_t out[OUTSIZE * OUTSIZE * 3];
static error_t errors_exp[ERRSIZE * ERRSIZE];
static error_t errors_res[ERRSIZE * ERRSIZE];
static error_t integral[(INSIZE + 1) * (INSIZE + 1)];
static image_t in_img = {in, INSIZE, INSIZE, 3};
static image_t out_img = {out, OUTSIZE, OUTSIZE, 3};

void fill_error_matrix(image_t in_, image_t out_, int orow, int ocol, error_t *errors, error_t *integral, int lrow, int lcol, int arow, int acol, int blocksize, int overlap, int num_blocks);
void generate_integral_image(image_t in, error_t *out);

void prepare_in_img()
{
    for (int i = 0; i < INSIZE; i++)
        for (int j = 0; j < INSIZE * 3; j++)
            // I found doing that not random is quite useful for visual checks
            in[i * (INSIZE * 3) + j] = MIN(i + j / 3, 255);
}

/* basically tell what blocks should be put at that position, array stops with -1
 * havoc all overlapping regions with random data.
 * cannot be called before prepare_in_img!
 */
void prepare_out_img(int rows[], int cols[])
{
    for (int i = 0; i < NUM_BLOCKS; i++)
    {
        for (int j = 0; j < NUM_BLOCKS; j++)
        {
            int inrow = rows[i * NUM_BLOCKS + j];
            int incol = cols[i * NUM_BLOCKS + j];
            if (inrow == -1 || incol == -1)
                return; /* reached end of list */
            int outrow = i * (BLOCKSIZE - OVERLAP);
            int outcol = j * (BLOCKSIZE - OVERLAP);
            for (int k = 0; k < BLOCKSIZE; k++)
                for (int m = 0; m < BLOCKSIZE * 3; m++)
                    out[(outrow + k) * (OUTSIZE * 3) + outcol * 3 + m] = in[(inrow + k) * (INSIZE * 3) + incol * 3 + m];
            if (j != 0)
                for (int k = BLOCKSIZE - OVERLAP; k < BLOCKSIZE; k++)
                    for (int m = 0; m < OVERLAP * 3; m++)
                        out[(outrow + k) * (OUTSIZE * 3) + outcol * 3 + m] = rand() % 256;
            if (i != 0)
                for (int k = 0; k < OVERLAP; k++)
                    for (int m = (BLOCKSIZE - OVERLAP) * 3; m < BLOCKSIZE * 3; m++)
                        out[(outrow + k) * (OUTSIZE * 3) + outcol * 3 + m] = rand() % 256;

        }
    }
}

void basic_fill_error_implementation(image_t img_in, image_t img_out, int orow, int ocol, error_t *errors)
{
    for (int i = 0; i < ERRSIZE; i++)
    {
        for (int j = 0; j < ERRSIZE; j++)
        {
            /* have this mental picture in your head
             * 0 . .  1 1 1  2 . .
             * . . .  . . .  2 . .
             * . . .  . . .  2 . .
             */
            errors[i * ERRSIZE + j] = 0;
            if (orow != 0 && ocol != 0)
            { // add 1 & 2 -> subtract 0
                slice_t slice_in = slice_image(img_in, i, j, i + OVERLAP, j + OVERLAP);
                slice_t slice_out = slice_image(img_out, orow, ocol, orow + OVERLAP, ocol + OVERLAP);
                errors[i * ERRSIZE + j] -= l2norm(slice_in, slice_out);
            }
            if (orow != 0)
            { // add 1
                slice_t slice_in = slice_image(img_in, i, j, i + OVERLAP, j + BLOCKSIZE);
                slice_t slice_out = slice_image(img_out, orow, ocol, orow + OVERLAP, ocol + BLOCKSIZE);
                errors[i * ERRSIZE + j] += l2norm(slice_in, slice_out);
            }
            if (ocol != 0)
            { // add 2
                slice_t slice_in = slice_image(img_in, i, j, i + BLOCKSIZE, j + OVERLAP);
                slice_t slice_out = slice_image(img_out, orow, ocol, orow + BLOCKSIZE, ocol + OVERLAP);
                errors[i * ERRSIZE + j] += l2norm(slice_in, slice_out);
            }
        }
    }
}

/* test case where next block with col = 0 */
void test_fill_error_matrix_left_border()
{
    prepare_in_img();
    generate_integral_image(in_img, integral);
    // fill the first row (numbers where randomly selected)
    int rows[] = {16, 59, 74, -1};
    int cols[] = {46, 11, 37, -1};
    prepare_out_img(rows, cols);
    // second row of blocks, first col of blocks
    int orow = 1 * (BLOCKSIZE - OVERLAP);
    int ocol = 0 * (BLOCKSIZE - OVERLAP);

    // INTMIN -> should not be used (hopefully undefined behavior gets caught)
    basic_fill_error_implementation(in_img, out_img, orow, ocol, errors_exp);
    fill_error_matrix(in_img, out_img, orow, ocol, errors_res, integral, INT_MIN, INT_MIN, rows[0], cols[0], BLOCKSIZE, OVERLAP, NUM_BLOCKS);
    int count = 0;
    for (int i = 0; i < ERRSIZE; i++)
        for (int j = 0; j < ERRSIZE; j++)
        {
            if (!TEST_CHECK(errors_exp[i * ERRSIZE + j] == errors_res[i * ERRSIZE + j]))
            {
                count++;
                TEST_MSG(
                    "index (%d, %d), expected %d, got %d\n",
                    i, j, errors_exp[i * ERRSIZE + j], errors_res[i * ERRSIZE + j]);
            }
        }
    if (count)
        TEST_MSG("encountered %d errors\n", count);
}

/* test case where next block with row = 0 */
void test_fill_error_matrix_top_border()
{
    prepare_in_img();
    generate_integral_image(in_img, integral);
    // fill the first row up to last entry
    int rows[] = {16, 59, -1};
    int cols[] = {46, 11, -1};
    prepare_out_img(rows, cols);
    // first row of blocks, third col of blocks
    int orow = 0 * (BLOCKSIZE - OVERLAP);
    int ocol = 2 * (BLOCKSIZE - OVERLAP);

    basic_fill_error_implementation(in_img, out_img, orow, ocol, errors_exp);
    // INTMIN -> should not be used (hopefully undefined behavior gets caught)
    fill_error_matrix(in_img, out_img, orow, ocol, errors_res, integral, rows[1], cols[1], INT_MIN, INT_MIN, BLOCKSIZE, OVERLAP, NUM_BLOCKS);
    int count = 0;
    for (int i = 0; i < ERRSIZE; i++)
        for (int j = 0; j < ERRSIZE; j++)
        {
            TEST_ASSERT(errors_exp[i * ERRSIZE + j] == errors_res[i * ERRSIZE + j]);
            if (!TEST_CHECK(errors_exp[i * ERRSIZE + j] == errors_res[i * ERRSIZE + j]))
            {
                count++;
                TEST_MSG(
                    "index (%d, %d), expected %d, got %d\n",
                    i, j, errors_exp[i * ERRSIZE + j], errors_res[i * ERRSIZE + j]);
            }
        }
    if (count)
        TEST_MSG("encountered %d errors\n", count);
}

/* test case where next block with row != 0 && col = num_blocks-1 */
void test_fill_error_matrix_right_border()
{
    prepare_in_img();
    generate_integral_image(in_img, integral);
    // fill up first row, and second row up (and without) to last
    int rows[] = {16, 59, 15, 66, 62, -1};
    int cols[] = {46, 11, 71, 43, 23, -1};
    prepare_out_img(rows, cols);
    // second row of blocks, third col of blocks
    int orow = 1 * (BLOCKSIZE - OVERLAP);
    int ocol = 2 * (BLOCKSIZE - OVERLAP);

    basic_fill_error_implementation(in_img, out_img, orow, ocol, errors_exp);
    // INTMIN -> should not be used (hopefully undefined behavior gets caught)
    fill_error_matrix(in_img, out_img, orow, ocol, errors_res, integral, rows[4], cols[4], rows[2], cols[2], BLOCKSIZE, OVERLAP, NUM_BLOCKS);
    int count = 0;
    for (int i = 0; i < ERRSIZE; i++)
        for (int j = 0; j < ERRSIZE; j++)
        {
            if (!TEST_CHECK(errors_exp[i * ERRSIZE + j] == errors_res[i * ERRSIZE + j]))
            {
                count++;
                TEST_MSG(
                    "index (%d, %d), expected %d, got %d\n",
                    i, j, errors_exp[i * ERRSIZE + j], errors_res[i * ERRSIZE + j]);
            }
        }
    if (count)
        TEST_MSG("encountered %d errors\n", count);
}

/* test case where next block with row != 0 && col != num_blocks-1 && col != 0 */
void test_fill_error_matrix_no_border()
{
    prepare_in_img();
    generate_integral_image(in_img, integral);
    // fill up first row, first element of second one
    int rows[] = {16, 59, 15, 66, -1};
    int cols[] = {46, 11, 71, 43, -1};
    prepare_out_img(rows, cols);
    // second row of blocks, second col of blocks
    int orow = 1 * (BLOCKSIZE - OVERLAP);
    int ocol = 1 * (BLOCKSIZE - OVERLAP);

    basic_fill_error_implementation(in_img, out_img, orow, ocol, errors_exp);
    fill_error_matrix(in_img, out_img, orow, ocol, errors_res, integral, rows[3], cols[3], rows[1], cols[1], BLOCKSIZE, OVERLAP, NUM_BLOCKS);
    int count = 0;
    for (int i = 0; i < ERRSIZE; i++)
        for (int j = 0; j < ERRSIZE; j++)
        {
            if (!TEST_CHECK(errors_exp[i * ERRSIZE + j] == errors_res[i * ERRSIZE + j]))
            {
                count++;
                TEST_MSG(
                    "index (%d, %d), expected %d, got %d\n",
                    i, j, errors_exp[i * ERRSIZE + j], errors_res[i * ERRSIZE + j]);
            }
        }
    if (count)
        TEST_MSG("encountered %d errors\n", count);
}

void test_generate_integral_image()
{
    error_t output[(in_gen_integral.height + 1) * (in_gen_integral.width + 1)];
    generate_integral_image(in_gen_integral, output);
    for (int i = 0; i < (in_gen_integral.height + 1) * (in_gen_integral.width + 1); i++)
    {
        if (!TEST_CHECK(out_gen_integral[i] == output[i]))
            TEST_MSG(
                "index %i, expected %d, got %d, diff %d\n",
                i, out_gen_integral[i], output[i], out_gen_integral[i] - output[i]);
    }
}

void test_find_singleton()
{
    error_t errors[] = {1};
    coord xy = find(errors, 1, 1, 100, 1);
    TEST_CHECK(xy.row == 0 && xy.col == 0);
}

void test_find_tolerance()
{
    error_t errors[] = {
        10, 20, 30,
        40, 50, 60,
        70, 80, 12};
    coord c;
    int pos0, pos8;
    for (int i = 0; i < 10; i++)
    {
        c = find(errors, 3, 3, 3, 10);
        pos0 = (c.row == 0 && c.col == 0);
        pos8 = (c.row == 2 && c.col == 2);
        TEST_CHECK(pos0 || pos8);
    }
}

TEST_LIST = {
    {"find_singleton", test_find_singleton},
    {"find_tolerance", test_find_tolerance},
    {"generate_integral_image", test_generate_integral_image},
    {"fill_error_matrix_left_border", test_fill_error_matrix_left_border},
    {"fill_error_matrix_top_border", test_fill_error_matrix_top_border},
    {"fill_error_matrix_right_border", test_fill_error_matrix_right_border},
    {"fill_error_matrix_no_border", test_fill_error_matrix_no_border},
    {NULL, NULL}};