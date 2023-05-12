#include <stdio.h>

#include "acutest.h"
#include "image_quilting.h"

static pixel_t data[] =
    {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        0, 1, 2};
static const int size = sizeof(data) / sizeof(pixel_t);
static const int width = 3;
static const int height = size / width;
static const int channels = 1;
static image_t img = {data, width, height, channels};
static slice_t sli = {data, width, height, channels, channels *width};

void test_imwrite_imread()
{
    char *path = "temp.png";
    imwrite(img, path);
    image_t res = imread(path);
    remove(path);
    for (int i = 0; i < size; i++)
        if (!TEST_CHECK(data[i] == res.data[i]))
            TEST_MSG(
                "index %i, expected %u, got %u, diff %d\n",
                i, data[i], res.data[i], data[i] - res.data[i]);
    imfree(res);
}

// slices are really dependent so far on having correct input
void test_slice_image_total()
{
    slice_t res = slice_image(img, 0, 0, height, width);
    TEST_CHECK(res.data == data);
    TEST_CHECK(res.width == sli.width);
    TEST_CHECK(res.height == sli.height);
    TEST_CHECK(res.channels == sli.channels);
    TEST_CHECK(res.jumpsize == sli.jumpsize);
}

void test_slice_image_subsection()
{
    slice_t res = slice_image(img, 1, 1, height - 1, width);
    TEST_CHECK(res.data == (data + 1 * sli.jumpsize + 1));
}

void test_slice_slice()
{
    slice_t inp = slice_image(img, 1, 1, height - 1, width);
    slice_t res = slice_slice(inp, 1, 1, inp.height, inp.width);
    TEST_CHECK(res.data == inp.data + 1 * sli.jumpsize + 1);
}

TEST_LIST = {
    {"imwrite_imread", test_imwrite_imread},
    {"slice_image_total", test_slice_image_total},
    {"slice_image_subsection", test_slice_image_subsection},
    {"slice_slice", test_slice_slice},
    {NULL, NULL}};