#include "acutest.h"
#include "image_quilting.h"
// only use this readonly (there is no way to express this in C)
static pixel_t data1[] = {7.0, 6.0, 1.0, 0.0, 0.0, 3.0, 4.0, 2.0, 2.0, 8.0, 8.0, 8.0, 2.0, 3.0, 1.0};
static pixel_t data2[] = {5.0, 7.0, 4.0, 1.0, 3.0, 0.0, 2.0, 1.0, 4.0, 5.0, 7.0, 6.0, 0.0, 0.0, 0.0};
static pixel_t exp_res[] = {7.0, 7.0, 4.0, 1.0, 3.0, 0.0, 4.0, 1.0, 4.0, 8.0, 7.0, 6.0, 2.0, 3.0, 0.0};
static const int size = sizeof(data1) / sizeof(pixel_t);
static const int width = 3;
static const int height = size / width;
static const int channels = 1;
static const int jumpsize = width * channels;

/* a bit of a complicated test just to check execution
 *  data1  data2  error  path   out
 *  7 6 1  5 7 4  4 1 9  0 x 0  7 7 4
 *  0 0 3  1 3 0  1 9 9  x 0 0  1 3 0
 *  4 2 2  2 1 4  4 1 4  0 x 0  4 1 4
 *  8 8 8  5 7 6  9 1 4  0 x 0  8 7 6
 *  2 3 1  0 0 0  4 9 1  0 0 x  2 3 0
 */
// void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int horizontal)
void test_dpcut_normal()
{
    slice_t s1 = {data1, width, height, channels, jumpsize};
    slice_t s2 = {data2, width, height, channels, jumpsize};

    pixel_t res[size];
    slice_t out = {res, width, height, channels, jumpsize};
    dpcut(s1, s2, out, 0);
    for (int i = 0; i < size; i++)
    {
        if (!TEST_CHECK(IS_CLOSE(res[i], exp_res[i])))
            TEST_MSG(
                "index %i, expected %e, got %e, diff %e\n",
                i, exp_res[i], res[i], res[i] - exp_res[i]);
    }
}

void test_dpcut_transpose()
{
    pixel_t tdata1[size], tdata2[size], texp_res[size];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            tdata1[j * height + i] = data1[i * width + j];
            tdata2[j * height + i] = data2[i * width + j];
            texp_res[j * height + i] = exp_res[i * width + j];
        }
    }
    slice_t ts1 = {tdata1, height, width, channels, height * channels};
    slice_t ts2 = {tdata2, height, width, channels, height * channels};
    pixel_t res[size];
    slice_t out = {res, height, width, channels, height * channels};
    dpcut(ts1, ts2, out, 1);
    for (int i = 0; i < size; i++)
    {
        if (!TEST_CHECK(IS_CLOSE(res[i], texp_res[i])))
            TEST_MSG(
                "index %i, expected %e, got %e, diff %e\n",
                i, texp_res[i], res[i], res[i] - texp_res[i]);
    }
}

/* test the same thing but now with three channels instead of one */
void test_dpcut_channels()
{
    pixel_t cdata1[3 * size], cdata2[3 * size], cexp_res[3 * size], res[3 * size];
    // create the same pattern now with three colors
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
        {
            cdata1[i * (3 * width) + j * 3 + 0] = data1[i * width + j];
            cdata1[i * (3 * width) + j * 3 + 1] = data1[i * width + j];
            cdata1[i * (3 * width) + j * 3 + 2] = 0.0;
            cdata2[i * (3 * width) + j * 3 + 0] = data2[i * width + j];
            cdata2[i * (3 * width) + j * 3 + 1] = data2[i * width + j];
            cdata2[i * (3 * width) + j * 3 + 2] = 0.0;
            cexp_res[i * (3 * width) + j * 3 + 0] = exp_res[i * width + j];
            cexp_res[i * (3 * width) + j * 3 + 1] = exp_res[i * width + j];
            cexp_res[i * (3 * width) + j * 3 + 2] = 0.0;
        }

    slice_t ts1 = {cdata1, width, height, 3, width * 3};
    slice_t ts2 = {cdata2, width, height, 3, width * 3};
    slice_t out = {res, width, height, 3, width * 3};
    dpcut(ts1, ts2, out, 0);
    for (int i = 0; i < 3 * size; i++)
        if (!TEST_CHECK(IS_CLOSE(res[i], cexp_res[i])))
            TEST_MSG(
                "index %i, expected %e, got %e, diff %e\n",
                i, cexp_res[i], res[i], res[i] - cexp_res[i]);
}

void test_dpcut_singleton()
{
    pixel_t d1[1] = {5};
    pixel_t d2[1] = {7};
    pixel_t r[1];
    slice_t s1 = {d1, 1, 1, 1, 1};
    slice_t s2 = {d2, 1, 1, 1, 1};
    slice_t out = {r, 1, 1, 1, 1};
    dpcut(s1, s2, out, 0);
    TEST_CHECK(r[0] == d2[0]);
    dpcut(s1, s2, out, 1);
    TEST_CHECK(r[0] == d2[0]);
}

TEST_LIST = {
    {"dpcut_normal", test_dpcut_normal},
    {"dpcut_transpose", test_dpcut_transpose},
    {"dpcut_channels", test_dpcut_channels},
    {"dpcut_singleton", test_dpcut_singleton},
    {NULL, NULL}};