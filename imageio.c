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
    pixel_t * startptr;
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
void cpy_block(slice_t in_block, slice_t out_block){

    for (int i = 0; i < out_block.height; i++){
        for (int j = 0; j < out_block.channels * out_block.width; j++){

            out_block.data[i * out_block.jumpsize + j] = in_block.data[i * in_block.jumpsize + j];

        }
    }
}

/* 
Calculates the error between the slice of an input block (i.e., for each block of the input image) and the slice of
a fixed output block
*/
void calc_errors(image_t in, int bsize, slice_t in_slice, slice_t out_slice, pixel_t *errors, int add){

    int error_jumpsize = in.width - bsize;

    if(add){
        for (int i = 0; i < in.height - bsize; i++) {
            for (int j = 0; j < in.width - bsize; j++) {
                
                in_slice.data = in.data + in_slice.jumpsize*i + j*in_slice.channels;
                errors[i*error_jumpsize + j] += l2norm(in_slice, out_slice); // if add > 0, the operation is plus (+)


            }
        }
    }
    else{
        for (int i = 0; i < in.height - bsize; i++) {
            for (int j = 0; j < in.width - bsize; j++) {
            
                in_slice.data = in.data + in_slice.jumpsize*i + j*in_slice.channels;
                errors[i*error_jumpsize + j] -= l2norm(in_slice, out_slice); // if add == 0, the operation is minus (-)

            }
        }
    }
}

// Coordinates of a candidate block
typedef struct{
    int row;
    int col;
} coord;

/* 
The function find finds the coordinates of a candidate block that is within a certain range (specified by tolerance)
from the best fitting block
 */
coord find(pixel_t *errors, int height, int width, pixel_t tolerance){

    pixel_t min_error = INFINITY;

    // search for the minimum error in the errors array and store it in a variable.
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {

            if(errors[i*width + j] < min_error){
                min_error = errors[i*width + j];
            }
        }
    }

    // Count how many canditates exist in order to know the size of the array of candidates
    pixel_t tol_range = (1.0 + tolerance) * min_error;
    int nr_candidates = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {

            if(errors[i*width + j] <= tol_range){
                nr_candidates++;
            }
        }
    }

    // Create the array with candidates and populate it with
    coord * candidates = (coord*)malloc(nr_candidates * sizeof(coord));
    int idx = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {

            if(errors[i*width + j] <= tol_range){
                coord candid = {i,j};
                candidates[idx] = candid;
                idx++;
            }
        }
    }

    // Choose randomly a candidate
    srand(time(0));
    int random_idx = rand() % nr_candidates;
    coord random_candidate = candidates[random_idx];
    free(candidates);
    candidates = NULL;

    // return coordinates of the random candidate
    return random_candidate;
}



// Command to compile the code:   gcc imageio.c L2norm.c  -o imageio -lm
int main()
{   

    //floor
    int bsize = 35;
    int num_blocks_out = 10;
    int ovsize = floor(bsize/6);
    pixel_t tolerance = 0.3;

    //Open image
    image_t in = imread("floor.jpg");

    //printf("Height: %d\n", in.height);
    //printf("Width: %d\n", in.width);

    int out_size = num_blocks_out * (bsize - ovsize - 1);
    int n = out_size * out_size * 3;
    double *out_image = (double *)calloc(n, sizeof(double));
    assert(out_image);
    image_t out = {out_image, out_size, out_size, 3};


    for (int i = 0; i < num_blocks_out - 1; i++) { 
        for (int j = 0; j < num_blocks_out - 1; j++) { 
            
            int si = i*bsize - i*ovsize;
            int sj = j*bsize - j*ovsize;


            pixel_t *errors = (pixel_t *)calloc((out.height-bsize)*(out.width-bsize), sizeof(pixel_t));

            // Very first case, so pick one at random
            if(i == 0 && j == 0){
                srand(time(0));
                int row_idx = rand() % (in.height - bsize);
                int col_idx = rand() % (in.width - bsize);

                slice_t out_block = slice_image(out, si, sj, si+bsize, sj+bsize);
                slice_t in_block = slice_image(in, row_idx, col_idx, row_idx+bsize, col_idx+bsize);

                cpy_block(in_block, out_block);
                free(errors);
                errors = NULL;
                continue;
            }
            // Top row, so only check left edges
            else if(i == 0){

                slice_t out_slice = slice_image(out, si, sj, si+bsize, sj+ovsize);
                slice_t in_slice = slice_image(in, 0, 0, bsize, ovsize);
                calc_errors(in, bsize, in_slice, out_slice, errors, 1);
            }
            // Left col, so only check above
            else if(j == 0){
                
                slice_t out_slice = slice_image(out, si, sj, si+ovsize, sj+bsize);
                slice_t in_slice = slice_image(in, 0, 0, ovsize, bsize);
                calc_errors(in, bsize, in_slice, out_slice, errors, 1);

            }
            // Typical case, check above and left
            else{
                slice_t out_slice = slice_image(out, si, sj, si+bsize, sj+ovsize);
                slice_t in_slice = slice_image(in, 0, 0, bsize, ovsize);
                calc_errors(in, bsize, in_slice, out_slice, errors, 1);

                out_slice = slice_image(out, si, sj, si+ovsize, sj+bsize);
                in_slice = slice_image(in, 0, 0, ovsize, bsize);
                calc_errors(in, bsize, in_slice, out_slice, errors, 1);

                // Remove the overlapping region in above and left cases
                out_slice = slice_image(out, si, sj, si+ovsize, sj+ovsize);
                in_slice = slice_image(in, 0, 0, ovsize, ovsize);
                calc_errors(in, bsize, in_slice, out_slice, errors, 0);

            }

            // Search for a random candidate block among the best matching blocks (determined by the tolerance)
            coord random_candidate = find(errors, in.height-bsize, in.width-bsize, tolerance);
            int rand_row = random_candidate.row;
            int rand_col = random_candidate.col;

            slice_t out_block = slice_image(out, si, sj, si+bsize, sj+bsize);
            slice_t in_block = slice_image(in, rand_row, rand_col, rand_row+bsize, rand_col+bsize);
            cpy_block(in_block, out_block);
            free(errors);
            errors = NULL;

        }
    }

    // Save the NO CUT version
    imwrite(out, "out.jpg");

    return 0;

// ======================================================================================================== End Piero


    /*
    image_t in_image = imread("image.jpg");
    printf("Image width: %d\n", in_image.width);
    printf("Image height: %d\n", in_image.height);
    printf("Number of channels: %d\n", in_image.channels);
    slice_t slice = slice_image(in_image, 50, 50, 70, 70);
    printf("%d %d %d\n", slice.width, slice.height, slice.channels);
    int n = slice.width * slice.height * slice.channels;
    printf("%d\n", n);
    unsigned char *out_image = malloc(sizeof(unsigned char) * n);
    assert(out_image);
    for (int row = 0; row < slice.height; row++) {
        for (int col = 0; col < slice.width*slice.channels; col++) {
            out_image[row*(slice.width*slice.channels) + col] = (unsigned char) slice.data[row * slice.jumpsize + col];
        }
    }
    stbi_write_png("out1.jpg", slice.width, slice.height, slice.channels, out_image, slice.channels*slice.width);
    printf("%f\n", l2norm(slice, slice));
    pixel_t arr1[] = {1.0, 1.0, 1.0, 1.0};
    pixel_t arr2[] = {0.0, 0.0, 0.0, 0.0};
    slice_t s1 = {arr1, 2, 2, 1, 2};
    slice_t s2 = {arr2, 2, 2, 1, 2};
    printf("%f\n", l2norm(s1, s2));

    imwrite(in_image, "out.png");
    imfree(in_image);
    
    return 0;
    */
}

// gcc test.c -o test -lm
