#include <stdio.h>
#include <stdlib.h>

#include "acutest.h"
#include "image_quilting.h"

void calc_errors(image_t img, int blocksize, slice_t in_slice, slice_t out_slice, pixel_t *errors, int add);

static pixel_t data[] = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9};
static pixel_t data_out[] = {1}; // for slice_out
static pixel_t data_in[] = {1};  // for slice_in (should not be used)
static pixel_t exp_res[] = {
    0, 1, 4,
    9, 16, 25,
    36, 49, 64};

static const int size = sizeof(data) / sizeof(pixel_t);
static const int width = 3;
static const int heigth = size / width;
static const int channels = 1;
static const int jumpsize = width * channels;
static const int width_slice = 1;
static const int heigth_slice = 1;

static image_t img = {data, width, heigth, channels};
static slice_t slice_out = {data_out, heigth_slice, width_slice, channels, jumpsize};
static slice_t slice_in = {data_in, heigth_slice, width_slice, channels, jumpsize};
static const int blocksize = 1;

void test_calc_errors_plus()
{
    int errorlen = (img.height - blocksize + 1) * (img.width - blocksize + 1);
    pixel_t *errors = (pixel_t *)calloc(errorlen, sizeof(pixel_t));
    calc_errors(img, blocksize, slice_in, slice_out, errors, 1);
    for (int i = 0; i < errorlen; i++)
        if (!TEST_CHECK(errors[i] == exp_res[i]))
            TEST_MSG("index %d, expected %d, got %d\n",
                     i, exp_res[i], errors[i]);
    free(errors);
}

static pixel_t exp_res_neg[] = {
    0, -1, -4,
    -9, -16, -25,
    -36, -49, -64};

void test_calc_errors_neg()
{
    int errorlen = (img.height - blocksize + 1) * (img.width - blocksize + 1);
    pixel_t *errors = (pixel_t *)calloc(errorlen, sizeof(pixel_t));
    calc_errors(img, blocksize, slice_in, slice_out, errors, 0);
    for (int i = 0; i < errorlen; i++)
        if (!TEST_CHECK(errors[i] == exp_res_neg[i]))
            TEST_MSG("index %d, expected %d, got %d\n",
                     i, exp_res_neg[i], errors[i]);
    free(errors);
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
    {"calc_errors_plus", test_calc_errors_plus},
    {"calc_errors_minus", test_calc_errors_neg},
    {NULL, NULL}};