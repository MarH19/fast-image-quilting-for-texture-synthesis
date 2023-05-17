#include <stdio.h>
#include "image_quilting.h"
#include "tsc_x86.h"
#define NUM_RUNS 5

/* just a quick runtime check such that we can easily see improvements in runtime */
double rdtsc(image_t in, int blocksize, int num_blocks, int overlap, pixel_t nom, pixel_t denom)
{
    int i, num_runs;
    myInt64 cycles, start;
    num_runs = NUM_RUNS;

    start = start_tsc();
    for (i = 0; i < num_runs; ++i)
        image_quilting(in, blocksize, num_blocks, overlap, nom, denom);

    cycles = stop_tsc(start) / num_runs;
    return (double)cycles;
}
int main()
{
    image_t in = imread("data/floor.jpg");
    int blocksize = 32;
    int num_blocks = 20;
    int overlap = 8;
    pixel_t nom = 3;
    pixel_t denom = 10;
    double runtime = rdtsc(in, blocksize, num_blocks, 8, nom, denom);
    printf("%.0f cycles, avg of %d runs\n", runtime, NUM_RUNS);
    return 0;
}