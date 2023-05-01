#include "image_quilting.h"
#define FREQUENCY 2.3e9

int main()
{
    // TO COUNT FLOPS:
    // nb = n_blocks, ih = in_height, iw = in_width, ov = overlap, bs = blocksize
    // = nb*(nb-1) * 2 * ((ih - bs + 1) * (iw -bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs -1) * 5) + (nb-1)^2 * (ov * ov * 3 * 3 + 1) + (nb * nb -1) * ((ih -bs + 1)*(iw - bs + 1) * 3 + 2)
    
    
    // measure different input sizes with following fixed parameters
    int blocksize = 35;
    int num_blocks = 30;
    int overlap = blocksize / 6;
    pixel_t tolerance = 0.3;
    
    image_t in2 = imread("data/internet/100KB");
    image_t in3 = imread("data/internet/500KB");
    image_t in4 = imread("data/internet/1MB");
    image_t in5 = imread("data/internet/2.5MB");
    image_t arr [] = {in2,in3,in4,in5};
    FILE *fp = freopen("benchmark.csv", "w", stdout);
    image_t t;
    int n;
    for(int i=0;i<4;i++){
        t = arr[i];
        n = in.width * in.height
        double r = rdtsc(in, blocksize, num_blocks, overlap, tolerance);
        // save n, time (in seconds), cycles into the csv
        printf("%d,%lf,%lf\n",n, r / FREQUENCY,r);
    }
    fclose(fp);
    // image_t out = image_quilting(in, blocksize, num_blocks, overlap, tolerance);
    // returns amount of cycles spent for computation
    
    
    //imwrite(out, "data/out.jpg");
}
// */
// gcc -O3 -mfma -fno-tree-vectorize -ffp-contract=fast timing.c image_quilting.c dpcut.c imageio.c L2norm.c -o timing -lm


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