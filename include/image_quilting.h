#include <stdint.h>

#ifndef IMAGE_QUILTING_HEADER
#define IMAGE_QUILTING_HEADER
typedef uint32_t pixel_t;

typedef struct
{
    pixel_t *data;
    int width;
    int height;
    int channels;
} image_t;

// struct that allows for better acces of slices
typedef struct
{
    pixel_t *data;
    int width;
    int height;
    int channels;
    int jumpsize;
} slice_t;

typedef struct
{
    int row;
    int col;
} coord;

/* imageio.c */
image_t imread(char *path);              /* read image at path */
void imwrite(image_t image, char *path); /* write image to path*/
void imfree(image_t image);              /* free memory of the image */
void slice_cpy(slice_t in, slice_t out); /* copy slice from in into out */
slice_t slice_image(image_t image, int start_row, int start_col, int end_row, int end_col);
slice_t slice_slice(slice_t sin, int start_row, int start_col, int end_row, int end_col);

/* l2norm.c */
pixel_t l2norm(slice_t, slice_t); /* calculate L2 norm between the two slices */

/* dpcut.c */
/* dpcut: calulate minimal cut between slice_1 and slice_2 and put it into out */
void dpcut(slice_t slice_1, slice_t slice_2, slice_t out, int left2right);

/* image_quilting.c */
image_t image_quilting(image_t in, int blocksize, int num_blocks, int overlap, pixel_t tol_nom, pixel_t tol_den);
coord find(pixel_t *errors, int height, int width, pixel_t tol_nom, pixel_t tol_den);

#define REL_TOL 1e-9
#define ABS_TOL 0.0
#define ABS(a) (((a) < 0) ? -1 * (a) : (a))

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#define IS_CLOSE(a, b) (ABS(a - b) <= MAX(REL_TOL * MAX(ABS(a), ABS(b)), ABS_TOL))

#endif /* #ifndef IMAGE_QUILTING_HEADER */