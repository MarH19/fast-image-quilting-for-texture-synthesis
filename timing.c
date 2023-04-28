#include "image_quilting.h"
#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "tsc_x86.h"

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define CALIBRATE
#define FREQUENCY 2.3e9

double rdtsc(image_t in, int blocksize, int num_blocks, int overlap, pixel_t tolerance)
{
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i)
        {
            image_quilting(in, blocksize, num_blocks, overlap, tolerance);
        }
        cycles = stop_tsc(start);

        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i)
    {
        image_quilting(in, blocksize, num_blocks, overlap, tolerance);
    }

    cycles = stop_tsc(start) / num_runs;
    return (double)cycles;
}
int main()
{
    // TO COUNT FLOPS:
    // nb = n_blocks, ih = in_height, iw = in_width, ov = overlap, bs = blocksize
    // = nb*(nb-1) * 2 * ((ih - bs + 1) * (iw -bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs -1) * 5) + (nb-1)^2 * (ov * ov * 3 * 3 + 1) + (nb * nb -1) * ((ih -bs + 1)*(iw - bs + 1) * 3 + 2)

    int blocksize = 35;
    int num_blocks = 30;
    int overlap = blocksize / 6;
    pixel_t tolerance = 0.3;

    image_t in = imread("floor.jpg");
    // image_t out = image_quilting(in, blocksize, num_blocks, overlap, tolerance);
    double r = rdtsc(in, blocksize, num_blocks, overlap, tolerance);
    printf("Time (in seconds): %lf \n", r / FREQUENCY);
    // imwrite(out, "out.jpg");
}
// gcc -O3 -mfma -fno-tree-vectorize -ffp-contract=fast timing.c image_quilting.c dpcut.c imageio.c L2norm.c -o timing -lm