#include "acutest.h"
#include "image_quilting.h"
// only use this readonly (there is no way to express this in C)
static pixel_t data1[] = {7.0, 6.0, 1.0, 0.0, 0.0, 3.0, 4.0, 2.0, 2.0, 8.0, 8.0, 8.0, 2.0, 3.0, 1.0};
static pixel_t data2[] = {5.0, 7.0, 4.0, 1.0, 3.0, 0.0, 2.0, 1.0, 4.0, 5.0, 7.0, 6.0, 0.0, 0.0, 0.0};
static pixel_t exp_res[] = {7.0, 7.0, 4.0, 1.0, 3.0, 0.0, 4.0, 1.0, 4.0, 8.0, 7.0, 6.0, 2.0, 3.0, 0.0};
const static int size = sizeof(data1) / sizeof(pixel_t);
const static int width  = 3;
const static int height = size / width;
const static int channels = 1;
const static int jumpsize = width * channels;


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
            tdata1[j*height + i] = data1[i*width + j];
            tdata2[j*height + i] = data2[i*width + j];
            texp_res[j*height + i] = exp_res[i*width + j];
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


TEST_LIST = {
    {"dpcut_normal", test_dpcut_normal},
    {"dpcut_transpose", test_dpcut_transpose},
    {NULL, NULL}};