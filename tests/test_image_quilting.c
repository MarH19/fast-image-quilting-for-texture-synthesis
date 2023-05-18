#include <stdio.h>
#include <stdlib.h>

#include "acutest.h"
#include "image_quilting.h"

void fill_error_matrix(image_t in_, image_t out_, int orow, int ocol, pixel_t *errors, pixel_t *integral, int lrow, int lcol, int arow, int acol, int blocksize, int overlap, int num_blocks);

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
    {NULL, NULL}};