#include "image_quilting.h"
#include <stdio.h>
int main()
{   
    printf("test print");
    int blocksize = 35;
    int num_blocks = 10;
    int overlap = 5;
    pixel_t tolerance = 0.3;

    image_t in = imread("floor.jpg");

    image_t out = image_quilting(in, blocksize, num_blocks, overlap, tolerance);

    imwrite(out, "out.jpg");
}
// gcc test.c image_quilting.c dpcut.c imageio.c L2norm.c  -o test -lm