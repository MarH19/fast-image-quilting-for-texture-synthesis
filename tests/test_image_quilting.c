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
static const pixel_t out_gen_integral[] = {
    0, 0, 0, 0, 0,
    0, 3, 15, 42, 90,
    0, 78, 198, 372, 612,
    0, 321, 741, 1278, 196593};

/* preparation of data for fill_error_matrix */
#define SIZE 100
#define OVERLAP 8
#define BLOCKSIZE 24
static pixel_t in[SIZE * SIZE * 3]; /* 3 because of the colors */

void fill_error_matrix(image_t in_, image_t out_, int orow, int ocol, pixel_t *errors, pixel_t *integral, int lrow, int lcol, int arow, int acol, int blocksize, int overlap, int num_blocks);
void generate_integral_image(image_t in, pixel_t *out);

void prepare_in_img()
{
    for (int i = 0; i < SIZE * SIZE * 3; i++)
    {
        in[i] = rand() % 256;
    }
}

void test_generate_integral_image()
{
    pixel_t output[(in_gen_integral.height + 1)*(in_gen_integral.width + 1)];
    generate_integral_image(in_gen_integral, output);
    for (int i = 0; i < (in_gen_integral.height + 1) * (in_gen_integral.width + 1); i++)
    {
        if (!TEST_CHECK(out_gen_integral[i] == output[i]))
            TEST_MSG(
                "index %i, expected %u, got %u, diff %d\n",
                i, out_gen_integral[i], output[i], out_gen_integral[i] - output[i]);
    }
}

void test_find_singleton()
{
    pixel_t errors[] = {1};
    coord xy = find(errors, 1, 1, 100, 1);
    TEST_CHECK(xy.row == 0 && xy.col == 0);
}

void test_find_tolerance()
{
    pixel_t errors[] = {
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
    {NULL, NULL}};