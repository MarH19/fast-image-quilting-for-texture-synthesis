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
#include <assert.h>
double rdtsc(image_t in, int blocksize, int num_blocks, int overlap, pixel_t nom, pixel_t denom)
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
            image_quilting(in, blocksize, num_blocks, overlap, nom, denom);
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
        image_quilting(in, blocksize, num_blocks, overlap, nom, denom);
    }

    cycles = stop_tsc(start) / num_runs;
    return (double)cycles;
}

double flop_counter(int nb_i, int ih_i, int iw_i, int ov_i, int bs_i)
{
    double nb = (double)nb_i;
    double ih = (double)ih_i;
    double iw = (double)iw_i;
    double ov = (double)ov_i;
    double bs = (double)bs_i;
    return nb*nb * (2 * ((ov * bs) * 10 + ov * 2 + (bs -1) * 5) + (ih - bs + 1) * (iw -bs + 1) * ((ov * (bs-ov) * 3 * 2 + 7) + (ov * (bs- 2 * ov) * 3 * 2 + 7))) + 2 * nb*nb * (ov * ov * 3 * 3 + 1) + (nb * nb -1) * ((ih -bs + 1)*(iw - bs + 1) * 2 + 2);

}
image_t generate_image(int size){
        int n = size * size * 3;
        pixel_t *image = malloc(sizeof(pixel_t) * n);
        assert(image);
        for (int i = 0; i < n; i++)
        {
            image[i] = (pixel_t)(rand() % 256);
        }
        image_t in_image = {image, size, size, 3};

        return in_image;
    }

int main()
{
    // TO COUNT FLOPS:
    // nb = n_blocks, ih = in_height, iw = in_width, ov = overlap, bs = blocksize
    // = nb*(nb-1) * 2 * ((ih - bs + 1) * (iw -bs + 1) * (1 + ov * bs * 3 * 3 + 1) + (ov * bs) * 10 + ov * 2 + (bs -1) * 5) + (nb-1)^2 * (ov * ov * 3 * 3 + 1) + (nb * nb -1) * ((ih -bs + 1)*(iw - bs + 1) * 3 + 2)

    FILE *fp_in = fopen("output/benchmark_inputsize_opt2", "w");
    int inp_sizes[] = {63,79,95,111,127,143,159,175,191,207,223,239};
    int block_size = 48;
    int overlap = 8;
    pixel_t nom = 3;
    pixel_t denom = 10;
    pixel_t tolerance = 0.3;
    int num_blocks = 8;

    for (int i=0;i<12;i++){
        image_t in = generate_image(inp_sizes[i]);
        double flops = flop_counter(num_blocks,inp_sizes[i],inp_sizes[i],overlap,block_size);
        double r = rdtsc(in, block_size, num_blocks, overlap, nom,denom);
        // save flops,input size,time (in seconds), cycles into the csv
        fprintf(fp_in,"%lf,%d,%lf,%lf\n", flops, inp_sizes[i], r / FREQUENCY, r);
        printf("%lf,%d,%lf,%lf\n", flops, inp_sizes[i], r / FREQUENCY, r);
    }    
    fclose(fp_in);
    printf("finished inputsize test \n");
    
    // measure different blocksizes
    FILE *fp_bs = fopen("output/benchmark_blocksize_opt2", "w");
    image_t in2 = generate_image(239);
    int ih2 = in2.height;
    int iw2 = in2.width;
    int nb = 10;
    double flops2;
    double r2;
    int block_sizes[] = {48,64,80,96,112,128,144};
    int overlap_sizes [] = {8,16,24,32,40,48,56}; 
    for (int i=0; i < 7; i++)
    {
        int bs = block_sizes[i];
        int ov = overlap_sizes[i];
        flops2 = flop_counter(nb,ih2,iw2,ov,bs);
        r2 = rdtsc(in2, bs, nb, ov, nom,denom);
        // save flops,blocksize,overlap, time (in seconds), cycles into the csv
        fprintf(fp_bs,"%lf,%d,%d,%lf,%lf\n", flops2, bs, ov, r2 / FREQUENCY, r2);
        printf("%lf,%d,%d,%lf,%lf\n", flops2, bs, ov, r2 / FREQUENCY, r2);
    }
    fclose(fp_bs);
    printf("finished blocksize test \n");
    
    // measure different num of blocks
    FILE *fp_nb = fopen("output/benchmark_numblocks_opt2", "w");
    image_t in3 = generate_image(239);
    int ih3 = in3.height;
    int iw3 = in3.width;
    int bs = 48;
    int ov = 8;
    double flops3;
    double r3;
    for (int i=5; i <= 30; i += 5)
    {
        flops3 = flop_counter(i,ih3,iw3,ov,bs);
        r3 = rdtsc(in3, bs, i, ov, nom,denom);
        // save flops,number of blocks, time (in seconds), cycles into the csv
        fprintf(fp_nb,"%lf,%d,%lf,%lf\n", flops3, i, r3 / FREQUENCY, r3);
        printf("%lf,%d,%lf,%lf\n", flops3, i, r3 / FREQUENCY, r3);
    }
    fclose(fp_nb);
    printf("finished num of blocks test \n");
    
}

// gcc -O3 -mfma -fno-tree-vectorize -ffp-contract=fast test.c src/timing.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c -o benchmark -lm

// gcc -I./include benchmark/test.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c -lm
// gcc -O3 -mfma -fno-tree-vectorize -ffp-contract=fast benchmark/test.c src/timing.c src/image_quilting.c src/dpcut.c src/imageio.c src/L2norm.c  -o benchmark -lm
// make benchmarks
// ./bin/benchmark/test.out