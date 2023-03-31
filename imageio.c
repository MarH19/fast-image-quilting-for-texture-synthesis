#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stdio.h>
#include "image_quilting.h"
#include "stb_image.h"
#include "stb_image_write.h"

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
        image[i] = (pixel_t) uimage[i];

    image_t in_image = {image, width, height, channels};

    stbi_image_free(uimage);
    return in_image;
}

void imwrite(image_t image, char *path) {
    int n = image.width * image.height * image.channels;

    unsigned char * out_image = malloc(sizeof(unsigned char) * n);
    assert(out_image);
    for (int i = 0; i < n; i++) {
        out_image[i] = (unsigned char) image.data[i];
    }
    stbi_write_png(path, image.width, image.height, image.channels, out_image, image.width * image.channels);
    free(out_image);
    out_image = NULL;

}

void imfree(image_t image) {
    free(image.data);
}

int main()
{
    image_t in_image = imread("image.jpg");
    printf("Image width: %d\n", in_image.width);
    printf("Image height: %d\n", in_image.height);
    printf("Number of channels: %d\n", in_image.channels);
    imwrite(in_image, "out.png");
    imfree(in_image);

    return 0;
}

// gcc test.c -o test -lm
