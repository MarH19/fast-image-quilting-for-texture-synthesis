typedef double pixel_t;

typedef struct
{
    pixel_t *data;
    int width;
    int height;
    int channels;
} image_t;

// struct that allows for better acces of slices
typedef struct{
    pixel_t * data;
    int width;
    int height;
    int channels;
    int jumpsize;
} slice_t;