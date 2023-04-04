#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stdio.h>
#include "image_quilting.h"
#include "stb_image.h"
#include "stb_image_write.h"

extern double l2norm(slice_t inp_slice, slice_t out_slice);
extern void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int c);
extern pixel_t* transpose(pixel_t *mat, int width, int height);

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

int main()
{
    image_t in_image = imread("image.jpg");

    printf("Image width: %d\n", in_image.width);
    printf("Image height: %d\n", in_image.height);
    printf("Number of channels: %d\n", in_image.channels);
    slice_t imslice1 = slice_image(in_image, 65, 50, 85, 70);
    slice_t imslice2 = slice_image(in_image, 15, 10, 35, 30);
    printf("%d %d %d\n", imslice1.width, imslice1.height, imslice1.channels);
    int n = imslice1.width * imslice1.height * imslice1.channels;
    printf("%d\n", n);

    // create the out array
    pixel_t out[1200] = {0.0};
    slice_t imslice_out = {out, 20, 20, 3, 60};
    dpcut(imslice1, imslice2, imslice_out, 0);

    unsigned char *out_image = malloc(sizeof(unsigned char) * n);
    assert(out_image);

    // write out
    for (int row = 0; row < imslice1.height; row++)
    {
        for (int col = 0; col < imslice1.width * imslice1.channels; col++)
        {
            out_image[row * (imslice1.width * imslice1.channels) + col] = (unsigned char)imslice_out.data[row * imslice_out.jumpsize + col];
        }
    }
    stbi_write_png("out1.jpg", imslice1.width, imslice1.height, imslice1.channels, out_image, imslice1.channels * imslice1.width);

    // write in1
    out_image = malloc(sizeof(unsigned char) * n);
    assert(out_image);
    for (int row = 0; row < imslice1.height; row++)
    {
        for (int col = 0; col < imslice1.width * imslice1.channels; col++)
        {
            out_image[row * (imslice1.width * imslice1.channels) + col] = (unsigned char)imslice1.data[row * imslice1.jumpsize + col];
        }
    }
    stbi_write_png("in1.jpg", imslice1.width, imslice1.height, imslice1.channels, out_image, imslice1.channels * imslice1.width);

    // write in2
    out_image = malloc(sizeof(unsigned char) * n);
    assert(out_image);
    for (int row = 0; row < imslice2.height; row++)
    {
        for (int col = 0; col < imslice2.width * imslice2.channels; col++)
        {
            out_image[row * (imslice2.width * imslice2.channels) + col] = (unsigned char)imslice2.data[row * imslice2.jumpsize + col];
        }
    }
    stbi_write_png("in2.jpg", imslice2.width, imslice2.height, imslice2.channels, out_image, imslice2.channels * imslice2.width);

    // // test l2norm
    // printf("%f\n", l2norm(slice, slice));
    // pixel_t arr1[] = {1.0, 1.0, 1.0, 1.0};
    // pixel_t arr2[] = {0.0, 0.0, 0.0, 0.0};
    // slice_t s1 = {arr1, 2, 2, 1, 2};
    // slice_t s2 = {arr2, 2, 2, 1, 2};
    // printf("%f\n\n", l2norm(s1, s2));

    // test dpcut
    // pixel_t matrix1[72] = {0.0};
    // pixel_t matrix2[72] = {1.0};
    // pixel_t out[72] = {0.0};

    // for (int i = 0; i < 72; i++){
    //     matrix2[i] = 1;
    // }
    // for (int i = 0; i < 6; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         for (int k = 0; k < 3; k++)
    //         {
    //             if (j == 1)
    //             {
    //                 matrix1[i * 4 * 3 + j*3 + k] = 0.5;
    //             }
    //             else if (j == 2)
    //             {
    //                 matrix1[i * 4 * 3 + j*3 + k] = 0.8;
    //             } else
    //             {
    //                 matrix1[i * 4 * 3 + j*3 + k] = 1.0;
    //             }
    //         }
    //     }
    // }
    
    // slice_t slice1 = {matrix1, 4, 6, 3, 12};
    // slice_t slice2 = {matrix2, 4, 6, 3, 12};
    // slice_t slice_out = {out, 4, 6, 3, 12};
    // for (int i = 0; i < slice1.height; i++)
    // {
    //     for (int j = 0; j < slice1.width*3; j += 3)
    //     {
    //         printf("[%f, %f, %f] ", slice1.data[i * slice1.jumpsize + j], slice1.data[i * slice1.jumpsize + j + 1], slice1.data[i * slice1.jumpsize + j + 2]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // for (int i = 0; i < slice2.height; i++)
    // {
    //     for (int j = 0; j < slice2.width * 3; j += 3)
    //     {
    //         printf("[%f, %f, %f] ", slice2.data[i * slice2.jumpsize + j], slice2.data[i * slice2.jumpsize + j + 1], slice2.data[i * slice2.jumpsize + j + 2]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // // for printing out array
    // dpcut(slice1, slice2, slice_out, 0);
    // for (int i = 0; i < slice_out.height; i++)
    // {
    //     for (int j = 0; j < slice_out.width * 3; j += 3)
    //     {
    //         printf("[%f, %f, %f] ", slice_out.data[i * slice_out.jumpsize + j], slice_out.data[i * slice_out.jumpsize + j + 1], slice_out.data[i * slice_out.jumpsize + j + 2]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    
    // exit(-1);

    // // for printing dp array
    // pixel_t *dp = dpcut(slice1, slice2, slice_out, 0);
    // for (int i = 0; i < slice1.height; i++)
    // {
    //     for (int j = 0; j < slice1.width; j++)
    //     {
    //         printf("%f ", dp[i * slice1.width + j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    
    // exit(-1);

    
    // // for printing errors array
    // pixel_t *errors = dpcut(slice1, slice2, slice_out, 0);
    // for (int i = 0; i < slice1.height; i++)
    // {
    //     for (int j = 0; j < slice1.width; j++)
    //     {
    //         printf("%f ", errors[i * slice1.width + j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    
    // exit(-1);
    //free(errors);

    // print input image
    imwrite(in_image, "out.png");
    imfree(in_image);

    return 0;
}

// gcc imageio.c L2norm.c dpcut.c -o test -lm

