#include "image_quilting.h"

double error;

double l2norm(slice_t inp_slice, slice_t out_slice){
    int i = inp_slice.data[0];
    int j = inp_slice.data[0];
    
    
    for(i;i<inp_slice.height;i++){
        for(j;j<(inp_slice.channels * inp_slice.width);j++){

            double inp_data = inp_slice.data[i*inp_slice.jumpsize+j];
            double out_data = out_slice.data[i*inp_slice.jumpsize+j];

            error += sqrt(pow(inp_data-out_data, 2));

        }
        
    }

    return error;
}