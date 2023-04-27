#include "image_quilting.h"
#include <stdio.h>

int main()
{
    int blocksize = 35;
    int num_blocks = 10;
    int overlap = blocksize / 6;
    pixel_t tolerance = 0.3;

    printf("num_blocks: %d\n", num_blocks);
    image_t in = imread("image.jpg");

    image_t out = image_quilting(in, blocksize, num_blocks, overlap, tolerance);

    imwrite(out, "out.jpg");
}
// gcc test.c image_quilting.c dpcut.c imageio.c L2norm.c  -o test -lm