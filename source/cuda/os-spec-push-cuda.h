#ifndef OS_SPEC_PUSH_CUDA_H
#define OS_SPEC_PUSH_CUDA_H

#include "os-param-cuda.h"
#include "os-vdf-cuda.h"
#include "os-current-cuda.h"
#include "os-emf-cuda.h"
#include "os-sys-cuda.h"
#include "os-chunk-pool-cuda.h"

void zero_jay_wrapper( const int *const prefix_size, const int avg_block_tile, 
                       const t_current_c *const jay, const int nblocks,
                       const cudaStream_t stream );

void shift_positions_wrapper( t_chunk ***chunks, int *d_block_start_tiles,
                              int *d_block_start_chunks, int *d_block_num_par,
                              int nmove0, int nmove1, int nmove2,
                              int nblocks, int chunk_size );

void dudt_boris_wrapper( t_chunk ***chunks, int *block_start_tiles,
                         int *block_start_chunks, int *block_num_par,
                         p_k_part rqm, t_emf_c *emf, int interpolation, bool grid_center,
                         p_double *energy, bool report_energy, p_double dt,
                         int nblocks, int chunk_size, cudaStream_t stream,
                         size_t shmem_size );

void advance_deposit_wrapper( t_chunk ***chunks, int *d_block_start_tiles,
                              int *d_block_start_chunks, int *d_block_num_par,
                              p_double *dx, t_current_c *jay, int interpolation,
                              p_double dt, int nblocks,
                              int chunk_size, int *n_hole, int *n_comm_ext,
                              int *n_recv_loc, cudaStream_t stream, size_t shmem_size );

#endif
