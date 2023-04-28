#include <stdlib.h>

#include "image_quilting.h"
#include "acutest.h"

#define REL_TOL 1e-9
#define ABS_TOL 0.0
#define ABS(a) ((a) < 0) ? -1 * (a) : (a)
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define IS_CLOSE(a, b) ABS(a - b) <= MAX(REL_TOL *MAX(ABS(a), ABS(b)), ABS_TOL)

/*
typedef struct
{
    pixel_t *data;
    int width;
    int height;
    int channels;
    int jumpsize;
} slice_t;
*/
void test_l2norm_duplicate(void)
{
    pixel_t arr[] = {1.0, 2.0, 3.0, 4.0};
    slice_t s1 = {arr, 2, 2, 1, 2};
    double res = l2norm(s1, s1);
    TEST_CHECK(res == 0.0);
}

void test_l2norm_jumpsize(void)
{
    // we want to generate 2x2 arrays
    pixel_t arr1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    pixel_t arr2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    slice_t s1 = {arr1, 2, 2, 1, 4};
    slice_t s2 = {arr2, 2, 2, 1, 12};
    double res = l2norm(s1, s2);
    TEST_ASSERT(IS_CLOSE(res, 116.0));
}

void test_l2norm_reverse(void)
{
    pixel_t arr1[] = {1.0, 2.0, 3.0, 4.0};
    pixel_t arr2[] = {1.0, 0.0, 0.0, 0.0};
    slice_t s1 = {arr1, 2, 2, 1, 2};
    slice_t s2 = {arr2, 2, 2, 1, 2};
    double res1 = l2norm(s1, s2);
    double res2 = l2norm(s2, s1);
    TEST_CHECK(res1 == res2);
    TEST_ASSERT(IS_CLOSE(res1, 29.0));
}

void test_l2norm_channels(void)
{
    pixel_t arr1[] = {0.0, 1.0, 2.0, 2.0, 1.0, 0.0};
    pixel_t arr2[] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    slice_t s1 = {arr1, 2, 1, 3, 6};
    slice_t s2 = {arr2, 2, 1, 3, 6};
    double res = l2norm(s1, s2);
    TEST_ASSERT(IS_CLOSE(res, 5.0));
}

void test_fail(void)
{
    int a, b;

    /* This condition is designed to fail so you can see what the failed test
     * output looks like. */
    a = 1;
    b = 2;
    TEST_CHECK(a + b == 5);

    /* Here is TEST_CHECK_ in action. */
    TEST_CHECK_(a + b == 5, "%d + %d == 5", a, b);

    /* We may also show more information about the failure. */
    if (!TEST_CHECK(a + b == 5))
    {
        TEST_MSG("a: %d", a);
        TEST_MSG("b: %d", b);
    }

    /* The macro TEST_MSG() only outputs something when the preceding
     * condition fails, so we can avoid the 'if' statement. */
    TEST_CHECK(a + b == 3);
    TEST_MSG("a: %d", a);
    TEST_MSG("b: %d", b);
}

static void helper(void)
{
    /* Kill the current test with a condition which is never true. */
    TEST_ASSERT(1 == 2);

    /* This never happens because the test is aborted above. */
    TEST_CHECK(1 + 2 == 2 + 1);
}

void test_abort(void)
{
    helper();

    /* This test never happens because the test is aborted inside the helper()
     * function. */
    TEST_CHECK(1 * 2 == 2 * 1);
}

void test_crash(void)
{
    int *invalid = ((int *)NULL) + 0xdeadbeef;

    *invalid = 42;
    TEST_CHECK_(1 == 1, "This should never execute, due to a write into "
                        "an invalid address.");
}

TEST_LIST = {
    {"l2norm_duplicate", test_l2norm_duplicate},
    {"l2norm_jumpsize", test_l2norm_jumpsize},
    {"l2norm_reverse", test_l2norm_reverse},
    {"l2norm_channels", test_l2norm_channels},
    {NULL, NULL}};

// gcc tests/test_l2norm.c src/L2norm.c -I./include -o test -lm