#include "image_quilting.h"
#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "tsc_x86.h"
#include "timing.h"
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
