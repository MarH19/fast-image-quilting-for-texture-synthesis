#include <stdio.h>
#include <stdlib.h>

#include "acutest.h"
#include "image_quilting.h"

void test_find_singleton()
{
    pixel_t errors[] = {1.0};
    coord xy = find(errors, 1, 1, 100.0);
    TEST_CHECK(xy.row == 0 && xy.col == 0);
}

void test_find_tolerance()
{
    pixel_t errors[] = {
        10.0, 20.0, 30.0,
        40.0, 50.0, 60.0,
        70.0, 80.0, 12.0};
    coord c;
    int pos0, pos8;
    for (int i = 0; i < 10; i++) {
        c = find(errors, 3, 3, 0.3);
        pos0 = (c.row == 0 && c.col == 0);
        pos8 = (c.row == 2 && c.col == 2);
        TEST_CHECK(pos0 || pos8);
    }
}

TEST_LIST = {
    {"find_singleton", test_find_singleton},
    {"find_tolerance", test_find_tolerance},
    {NULL, NULL}};