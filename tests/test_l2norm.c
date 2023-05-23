#include <stdlib.h>
#include <math.h>

#include "acutest.h"
#include "image_quilting.h"

void test_l2norm_duplicate(void)
{
    pixel_t arr[] = {1, 2, 3, 4};
    slice_t s1 = {arr, 2, 2, 1, 2};
    error_t res = l2norm(s1, s1);
    TEST_CHECK(res == 0);
}

void test_l2norm_jumpsize(void)
{
    // we want to generate 2x2 arrays
    pixel_t arr1[] = {1, 2, 3, 4, 5, 6, 4, 2};
    pixel_t arr2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    slice_t s1 = {arr1, 2, 2, 1, 6};
    slice_t s2 = {arr2, 2, 2, 1, 14};
    error_t expected_res = 25;
    error_t res = l2norm(s1, s2);
    if (!TEST_CHECK(res == expected_res))
        TEST_MSG("[ERROR] %u != %u\n", res, expected_res);
}

void test_l2norm_reverse(void)
{
    pixel_t arr1[] = {1, 2, 3, 4};
    pixel_t arr2[] = {1, 0, 0, 0};
    slice_t s1 = {arr1, 2, 2, 1, 2};
    slice_t s2 = {arr2, 2, 2, 1, 2};
    error_t res1 = l2norm(s1, s2);
    error_t res2 = l2norm(s2, s1);
    if (!TEST_CHECK(res1 == res2))
        TEST_MSG("[ERROR] %u != %u\n", res1, res2);
    error_t expected_res = 29;
    if (!TEST_CHECK(res1 == expected_res))
        TEST_MSG("[ERROR] %u != %u\n", res1, expected_res);
}

void test_l2norm_channels(void)
{
    pixel_t arr1[] = {0, 1, 2, 2, 1, 0};
    pixel_t arr2[] = {1, 1, 1, 2, 2, 2};
    slice_t s1 = {arr1, 2, 1, 3, 6};
    slice_t s2 = {arr2, 2, 1, 3, 6};
    error_t res = l2norm(s1, s2);
    error_t expected_res = 7;
    if (!TEST_CHECK(res == expected_res))
        TEST_MSG("[ERROR] %u != %u\n", res, expected_res);
}

TEST_LIST = {
    {"l2norm_duplicate", test_l2norm_duplicate},
    {"l2norm_jumpsize", test_l2norm_jumpsize},
    {"l2norm_reverse", test_l2norm_reverse},
    {"l2norm_channels", test_l2norm_channels},
    {NULL, NULL}};

// gcc tests/test_l2norm.c src/L2norm.c -I./include -o test -lm