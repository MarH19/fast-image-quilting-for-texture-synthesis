#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stdio.h>
#include "image_quilting.h"
#include "stb_image.h"
#include "stb_image_write.h"

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

int main()
{
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
}

// gcc test.c -o test -lm
