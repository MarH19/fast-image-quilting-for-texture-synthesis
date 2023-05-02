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
#define MIN(a, b) ((a) < (b) ? (a) : (b))

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
double flop_counter(int nb_i,int ih_i, int iw_i, int ov_i, int bs_i){
    double nb = (double) nb_i;
    double ih = (double) ih_i;
    double iw = (double) iw_i;
    double ov = (double) ov_i;
    double bs = (double) bs_i;
    return nb * (nb - 1) * 2 * ((ih - bs + 1) * (iw - bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs - 1) * 5) + pow((nb - 1),2) * (ov * ov * 3 * 3 + 1) + (nb * nb - 1) * ((ih - bs + 1) * (iw - bs + 1) * 3 + 2);

}

int main()
{
    // TO COUNT FLOPS:
    // nb = n_blocks, ih = in_height, iw = in_width, ov = overlap, bs = blocksize
    // = nb*(nb-1) * 2 * ((ih - bs + 1) * (iw -bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs -1) * 5) + (nb-1)^2 * (ov * ov * 3 * 3 + 1) + (nb * nb -1) * ((ih -bs + 1)*(iw - bs + 1) * 3 + 2)

    // measure different input sizes with following fixed parameters
    // TODO check if nb is valid for the selected images
    /*
    int nb = 30;
    
    int ov = bs / 6;
    pixel_t tolerance = 0.3;
    
    image_t in2 = imread("data/internet/100KB.jpg");
    image_t in3 = imread("data/internet/500KB.jpg");
    image_t in4 = imread("data/internet/1MB.jpg");
    image_t in5 = imread("data/internet/2.5MB.jpg");
    image_t arr[] = {in2, in3, in4, in5};
    FILE *fp1 = freopen("output/benchmark_in.csv", "w", stdout);
    double flops1;
    double r;
    for (int i = 0; i < 4; i++)
    {
        image_t t = arr[i];
        int ih = t.height;
        int iw = t.width;
        int bs = 35;
        int n = t.width * t.height;
        flops1 = flop_counter();
        r = rdtsc(t, bs, nb, ov, tolerance);
        // save flops,n,time (in seconds), cycles into the csv
        printf("%d,%d,%lf,%lf\n", flops1, n, r / FREQUENCY, r);
    }
    fclose(fp1);
    */
    // new variant 
    
    int num_blocks = 2;
    pixel_t tolerance = 0.3;
    int block_size;
    FILE *fp_in = freopen("output/benchmark_in", "w", stdout);
    while(block_size < 4000){
        image_t in = imread("data/floor_51.jpg"); // initial image
        double ih = in.height;
        double iw = in.width;
        block_size = MIN(ih,iw)-1;
        double flops = flop_counter(num_blocks,ih,iw,block_size/2,block_size);
        double r = rdtsc(in, block_size, num_blocks, block_size/2, tolerance);
        // save flops,block_size,time (in seconds), cycles into the csv
        printf("%lf,%d,%lf,%lf\n", flops, block_size, r / FREQUENCY, r);
        // printf("performance (F/C):%lf \n",flops/r);
        image_t out = image_quilting(in, block_size, num_blocks, block_size/2, tolerance);
        imwrite(out, "data/floor_51.jpg");

    }
    fclose(fp_in);

    



    // measure different blocksizes
    image_t in2 = imread("data/floor.jpg");
    int ih = in2.height;
    int iw = in2.width;
    int nb = 30;
    FILE *fp2 = freopen("output/benchmark_blocksize", "w", stdout);
    double flops2;
    double r2;
    for (int i = 10; i <= 80; i += 10)
    {
        int bs = i;
        int ov = bs / 6;
        flops2 = flop_counter(nb,ih,iw,ov,bs);
        r2 = rdtsc(in2, bs, nb, ov, tolerance);
        // save flops,blocksize, time (in seconds), cycles into the csv
        printf("%lf,%d,%lf,%lf\n", flops2, bs, r2 / FREQUENCY, r2);
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