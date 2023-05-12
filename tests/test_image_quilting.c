#include <stdio.h>
#include <stdlib.h>

#include "acutest.h"
#include "image_quilting.h"

static pixel_t data1[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
static const int width = 3;
static const int size = sizeof(data1) / sizeof(pixel_t);
static const int heigth = size / width;
static const int channels = 1;
static image_t in = {data1, width, heigth, channels};
static const int jumpsize = width * channels;
static const int width_slice = 1;
static const int heigth_slice = 1;
static pixel_t data_out[] = {0};

static slice_t slice_out = {data_out, heigth_slice, width_slice, channels, jumpsize};
static slice_t slice_in = {NULL, heigth_slice, width_slice, channels, jumpsize};
static const int blocksize = 1;
void test_calc_errors_plus()
{
    int errorlen = (in.height - blocksize + 1) * (in.width - blocksize + 1) * sizeof(pixel_t);
    pixel_t *errors = (pixel_t *)malloc(errorlen);
    calc_errors(in, blocksize, slice_in, slice_out, errors, 1);
    int error_jumpsize = in.width - blocksize + 1;
    for (int i = 0; i < in.height - blocksize + 1; i++)
    {
        for (int j = 0; j < in.width - blocksize + 1; j++)
        {
            // check that the value is 1
            if(!TEST_CHECK(errors[i * error_jumpsize + j] == 1))
                TEST_MSG("index %i, %i, expected %u, got %u\n",
                i,j, 1, errors[i * error_jumpsize + j]);
        }
    }
    free(errors);
}

void test_calc_errors_neg() {}

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
    {"calc errors plus",test_calc_errors_plus},
    {NULL, NULL}};