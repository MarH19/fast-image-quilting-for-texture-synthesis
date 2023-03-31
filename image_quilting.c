#include "image_quilting.h"
extern image_t imread(char *path);
/* What do we need to do?
 *  1. load the image
 *  2. calculate all the parameters
 *  3. malloc the output image
 *  4. call function for random patch
 *  5. copy patches to something
 *  6. find next random patch
 *  7. 
 */

int image_quilting(char *input_file, char *output_file)
{
    image_t img = imread(input_file);
}