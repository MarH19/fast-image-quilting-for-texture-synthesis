#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef WIN32
#include <sys/time.h>
#endif

#include <time.h>
#include "tsc_x86.h"
#include "image_quilting.h"

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define CALIBRATE

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

    // measure different input sizes with following fixed parameters
    // TODO check if nb is valid for the selected images
    int nb = 30;
    int bs = 35;
    int ov = bs / 6;
    pixel_t tolerance = 0.3;
    image_t in2 = imread("data/internet/100KB.jpg");
    image_t in3 = imread("data/internet/500KB.jpg");
    image_t in4 = imread("data/internet/1MB.jpg");
    image_t in5 = imread("data/internet/2.5MB.jpg");
    image_t arr[] = {in2, in3, in4, in5};
    FILE *fp1 = freopen("output/benchmark_in.csv", "w", stdout);
    int flops1;
    double r;
    for (int i = 0; i < 4; i++)
    {
        image_t t = arr[i];
        int ih = t.height;
        int iw = t.width;
        int n = t.width * t.height;
        flops1 = nb * (nb - 1) * 2 * ((ih - bs + 1) * (iw - bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs - 1) * 5) + pow((nb - 1),2) * (ov * ov * 3 * 3 + 1) + (nb * nb - 1) * ((ih - bs + 1) * (iw - bs + 1) * 3 + 2);
        r = rdtsc(t, bs, nb, ov, tolerance);
        // save flops,n,time (in seconds), cycles into the csv
        printf("%d,%d,%lf,%lf\n", flops1, n, r / FREQUENCY, r);
    }
    fclose(fp1);

    // measure different blocksizes
    // TODO check if nb is valid for the selected image
    image_t in = imread("data/internet/500KB");
    int ih = in.height;
    int iw = in.width;
    FILE *fp2 = freopen("output/benchmark_blocksize.csv", "w", stdout);
    int flops2;
    double r2;
    for (int i = 10; i < 100; i += 10)
    {
        bs = i;
        ov = bs / 6;
        flops2 = nb * (nb - 1) * 2 * ((ih - bs + 1) * (iw - bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs - 1) * 5) + pow((nb - 1),2) * (ov * ov * 3 * 3 + 1) + (nb * nb - 1) * ((ih - bs + 1) * (iw - bs + 1) * 3 + 2);
        r2 = rdtsc(in, i, nb, ov, tolerance);
        // save flops,blocksize, time (in seconds), cycles into the csv
        printf("%d,%d,%lf,%lf\n", flops2, bs, r2 / FREQUENCY, r2);
    }
    fclose(fp2);
    // image_t out = image_quilting(in, blocksize, num_blocks, overlap, tolerance);
    // imwrite(out, "data/out.jpg");
}

// gcc -O3 -mfma -fno-tree-vectorize -ffp-contract=fast test.c src/timing.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c -o benchmark -lm

// gcc -I./include benchmark/test.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c -lm
// gcc -O3 -mfma -fno-tree-vectorize -ffp-contract=fast benchmark/test.c src/timing.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c  -o benchmark -lm
// make benchmark
// ./bin/benchmark/test.out