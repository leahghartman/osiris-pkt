#ifndef OS_SPEC_VALIDATE_CUDA_H
#define OS_SPEC_VALIDATE_CUDA_H

#include "os-chunk-pool-cuda.h"

  void validate_wrapper( t_chunk ***chunks, int *d_nchunks, int *my_nx_p_, int *block_start_tiles,
                         int *block_start_chunks, int *block_num_par, int chunk_size,
                         bool over_, int move_num0, int move_num1, int move_num2,
                         int *err, int my_aid, int nblocks);


#endif
