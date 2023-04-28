#include "image_quilting.h"

int main()
{
    int blocksize = 35;
    int num_blocks = 10;
    int overlap = blocksize / 6;
    pixel_t tolerance = 0.3;

    image_t in = imread("data/floor.jpg");
    image_t out = image_quilting(in, blocksize, num_blocks, overlap, tolerance);
    imwrite(out, "data/out.jpg");
}
// gcc -I./include src/test.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c -lm
// gcc src/test.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c  -o test -lm
// make tests
// ./bin/tests/test.out