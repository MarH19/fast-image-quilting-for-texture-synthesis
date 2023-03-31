#define STB_IMAGE_IMPLEMENTATION
#include <stdio.h>
#include "imread.h"
#include "stb_image.h"

image_int8 imread(char *path)
{
    int width, height, channels;
    unsigned char *image = stbi_load(path, &width, &height, &channels, 0);

    image_int8 in_image = {image, width, height, channels};

    stbi_image_free(image);

    return in_image;
}

int main()
{
    image_int8 in_image = imread("image.jpg");
    printf("Image width: %d\n", in_image.width);
    printf("Image height: %d\n", in_image.height);
    printf("Number of channels: %d\n", in_image.channels);

    return 0;
}

// gcc test.c -o test -lm
