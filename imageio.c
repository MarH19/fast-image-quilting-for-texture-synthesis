#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stdio.h>
#include "image_quilting.h"
#include "stb_image.h"
#include "stb_image_write.h"

// ======================================================================================================== Start Piero
#include <time.h>
// ======================================================================================================== End Piero

extern double l2norm(slice_t inp_slice, slice_t out_slice);
extern void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int c);
extern pixel_t *transpose(pixel_t *mat, int width, int height);

image_t imread(char *path)
{
    int width, height, channels, n;
    unsigned char *uimage;

    uimage = stbi_load(path, &width, &height, &channels, 0);
    assert(uimage);

    n = width * height * channels;
    pixel_t *image = malloc(sizeof(pixel_t) * n);
    assert(image);
    for (int i = 0; i < n; i++)
        image[i] = (pixel_t)uimage[i];

    image_t in_image = {image, width, height, channels};

    stbi_image_free(uimage);
    return in_image;
}

void imwrite(image_t image, char *path)
{
    int n = image.width * image.height * image.channels;

    unsigned char *out_image = malloc(sizeof(unsigned char) * n);
    assert(out_image);
    for (int i = 0; i < n; i++)
    {
        out_image[i] = (unsigned char)image.data[i];
    }
    stbi_write_png(path, image.width, image.height, image.channels, out_image, image.width * image.channels);
    free(out_image);
    out_image = NULL;
}

void imfree(image_t image)
{
    free(image.data);
}

/* slice_image: start inclusive, end exclusive */
slice_t slice_image(image_t image, int start_row, int start_col, int end_row, int end_col)
{
    int jumpsize;
    pixel_t *startptr;
    slice_t slice;

    jumpsize = image.width * image.channels;
    startptr = image.data + start_row * jumpsize + image.channels * start_col;

    slice.height = end_row - start_row;
    slice.width = end_col - start_col;
    slice.channels = image.channels;
    slice.data = startptr;
    slice.jumpsize = jumpsize;

    return slice;
}

slice_t slice_slice(slice_t sin, int start_row, int start_col, int end_row, int end_col)
{
    pixel_t *startptr;
    slice_t slice;

    startptr = sin.data + start_row * sin.jumpsize + sin.channels * start_col;
    slice.height = end_row - start_row;
    slice.width = end_col - start_col;
    slice.channels = sin.channels;
    slice.data = startptr;
    slice.jumpsize = sin.jumpsize;
    return slice;
}

// ======================================================================================================== Start Piero

/*
Copy block of input image to block of output image
*/
void slice_cpy(slice_t in, slice_t out)
{
    for (int i = 0; i < out.height; i++)
    {
        for (int j = 0; j < out.channels * out.width; j++)
        {
            out.data[i * out.jumpsize + j] = in.data[i * in.jumpsize + j];
        }
    }
}

/*
Calculates the error between the slice of an input block (i.e., for each block of the input image) and the slice of
a fixed output block
*/
void calc_errors(image_t in, int blocksize, slice_t in_slice, slice_t out_slice, pixel_t *errors, int add)
{
    int error_jumpsize = in.width - blocksize + 1;

    if (add)
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {

                in_slice.data = in.data + in_slice.jumpsize * i + j * in_slice.channels;
                errors[i * error_jumpsize + j] += l2norm(in_slice, out_slice); // if add > 0, the operation is plus (+)
            }
        }
    }
    else
    {
        for (int i = 0; i < in.height - blocksize + 1; i++)
        {
            for (int j = 0; j < in.width - blocksize + 1; j++)
            {

                in_slice.data = in.data + in_slice.jumpsize * i + j * in_slice.channels;
                errors[i * error_jumpsize + j] -= l2norm(in_slice, out_slice); // if add == 0, the operation is minus (-)
            }
        }
    }
}

/*
The function find finds the coordinates of a candidate block that is within a certain range (specified by tolerance)
from the best fitting block
 */
coord find(pixel_t *errors, int height, int width, pixel_t tolerance)
{
    pixel_t min_error = INFINITY;

    // search for the minimum error in the errors array and store it in a variable.
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (errors[i * width + j] < min_error)
            {
                min_error = errors[i * width + j];
            }
        }
    }

    // Count how many canditates exist in order to know the size of the array of candidates
    pixel_t tol_range = (1.0 + tolerance) * min_error;
    int nr_candidates = 0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (errors[i * width + j] <= tol_range)
            {
                nr_candidates++;
            }
        }
    }

    // Create the array with candidates and populate it with
    coord candidates[nr_candidates];
    int idx = 0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {

            if (errors[i * width + j] <= tol_range)
            {
                coord candid = {i, j};
                candidates[idx] = candid;
                idx++;
            }
        }
    }

    // Choose randomly a candidate
    int random_idx = rand() % nr_candidates;
    coord random_candidate = candidates[random_idx];

    // return coordinates of the random candidate
    return random_candidate;
}

// Command to compile the code:   gcc dpcut.c imageio.c L2norm.c  -o imageio -lm
int main()
{

    // floor
    int blocksize = 35;
    int num_blocks = 10;
    int overlap = blocksize / 6;
    pixel_t tolerance = 0.3;

    // Open image
    image_t in = imread("floor.jpg");

    // printf("Height: %d\n", in.height);
    // printf("Width: %d\n", in.width);

    int out_size = num_blocks * blocksize - (num_blocks - 1) * overlap;
    int n = out_size * out_size * 3;
    double *out_image = (double *)calloc(n, sizeof(double));
    assert(out_image);
    image_t out = {out_image, out_size, out_size, 3};
    int errorlen = (in.height - blocksize + 1) * (in.width - blocksize + 1) * sizeof(pixel_t);
    pixel_t *errors = (pixel_t *)malloc(errorlen);

    for (int row = 0; row < num_blocks; row++)
    {
        for (int col = 0; col < num_blocks; col++)
        {
            int si = row * (blocksize - overlap);
            int sj = col * (blocksize - overlap);
            memset(errors, 0, errorlen);

            // Very first case, so pick one at random
            if (row == 0 && col == 0)
            {
                int row_idx = rand() % (in.height - blocksize + 1);
                int col_idx = rand() % (in.width - blocksize + 1);

                slice_t out_block = slice_image(out, si, sj, si + blocksize, sj + blocksize);
                slice_t in_block = slice_image(in, row_idx, col_idx, row_idx + blocksize, col_idx + blocksize);

                slice_cpy(in_block, out_block);
                continue;
            }
            if (row != 0) // safe to check top
            {
                slice_t out_slice = slice_image(out, si, sj, si + overlap, sj + blocksize);
                slice_t in_slice = slice_image(in, 0, 0, overlap, blocksize);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }
            if (col != 0) // safe to check left
            {
                slice_t out_slice = slice_image(out, si, sj, si + blocksize, sj + overlap);
                slice_t in_slice = slice_image(in, 0, 0, blocksize, overlap);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 1);
            }
            if (col != 0 && row != 0) // overlap so remove common
            {
                slice_t out_slice = slice_image(out, si, sj, si + overlap, sj + overlap);
                slice_t in_slice = slice_image(in, 0, 0, overlap, overlap);
                calc_errors(in, blocksize, in_slice, out_slice, errors, 0);
            }

            // Search for a random candidate block among the best matching blocks (determined by the tolerance)
            coord random_candidate = find(errors, in.height - blocksize + 1, in.width - blocksize + 1, tolerance);
            int rand_row = random_candidate.row;
            int rand_col = random_candidate.col;

            slice_t out_block = slice_image(out, si, sj, si + blocksize, sj + blocksize);
            slice_t in_block = slice_image(in, rand_row, rand_col, rand_row + blocksize, rand_col + blocksize);
            // slice_cpy(in_block, out_block);

            if (row != 0)
            {

                dpcut(out_block, in_block, out_block, 1);
            }

            if (col != 0)
            {

                dpcut(out_block, in_block, out_block, 0);
            }
        }
    }

    free(errors);
    // Save the NO CUT version
    imwrite(out, "out/out.jpg");

    return 0;
    // ======================================================================================================== End Piero
}
// gcc dpcut.c imageio.c L2norm.c  -o imageio -lm