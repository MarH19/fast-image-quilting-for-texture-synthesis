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
