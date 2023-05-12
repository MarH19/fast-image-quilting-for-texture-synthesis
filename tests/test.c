#include "image_quilting.h"

int main()
{
    int blocksize = 35;
    int num_blocks = 10;
    int overlap = blocksize / 6;
    pixel_t tol_nom = 1;
    pixel_t tol_den = 3;

    image_t in = imread("data/floor.jpg");
    image_t out = image_quilting(in, blocksize, num_blocks, overlap, tol_nom, tol_den);
    imwrite(out, "data/out.jpg");
    imfree(in);
    imfree(out);
}
// gcc -I./include src/test.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c -lm
// gcc src/test.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c  -o test -lm
// make tests
// ./bin/tests/test.out