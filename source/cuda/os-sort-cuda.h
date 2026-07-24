#ifndef OS_SORT_CUDA_H
#define OS_SORT_CUDA_H

#include "os-vdf-cuda.h"
#include "os-chunk-pool-cuda.h"

void sort_wrapper( t_chunk ***chunks, t_chunk ***chunks_ext,
                   int *nchunks_max, int *nchunks,
                   int *nchunks_max_ext, int *nchunks_ext,
                   int *my_nx_p,
                   int *block_start_tiles, int *block_start_chunks,
                   int *block_num_par, int chunk_size, t_f0<int> *ihole,
                   int ihole_size, int *n_hole, t_f1<int> *neighbor_til_id, t_f2<int> *shift,
                   int *n_send_ext,
                   int *n_recv_loc, int *num_par, int nblocks, int *err );

void compactify_wrapper( t_chunk ***chunks, t_chunk ***chunks_ext,
                         int *nchunks_max, int *nchunks,
                         int *nchunks_max_ext, int *nchunks_ext,
                         int chunk_size,
                         t_f0<int> *ihole, int *n_hole, int *n_recv_ext, int *num_par,
                         int *n_recv_loc, int nblocks, int *err );

#endif
