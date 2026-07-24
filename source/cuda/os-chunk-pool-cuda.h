#ifndef OS_CHUNK_POOL_CUDA_H
#define OS_CHUNK_POOL_CUDA_H

#include <stdio.h>
#include "os-param-cuda.h"
#include "os-vdf-cuda.h"
#include "os-mem-transfer-manager-cuda.h"

// #define DEBUG_CHUNKS 1 // uncomment to enable debug functions
// #define CHUNK_POOL_DEBUG_NUM 123

typedef char t_chunk;

class t_chunk_pool {
public:
  int n_tils;
  int chunk_pool_size; // num chunks allocated in pool (max num of chunks avail on this node)
  int chunk_size;      // number of particles in each chunk

  // giant contiguous buffer to store all the chunks.
  t_chunk *chunks;

  // queue of available chunks, changes every get / return call
  t_chunk **chunk_queue;

  // first and last entries of the circular buffer. h_state is a contiguous 2-element array,
  // h_start and h_end point to these entries.
  int *h_state, *h_start, *h_end;
  size_t state_size;

  // device start / end buffers for accessing state on device -- must be synced
  // calls to get/return chunks functions on host / device. d_start, d_end, and d_end_tmp point
  // to the entries of d_state.
  int *d_state, *d_start, *d_end, *d_end_tmp;

  // offsets of where data is stored within each chunk, memory taken by a chunk
  int x_offset, p_offset, q_offset, ix_offset;
  size_t chunk_mem;

  // this is a chunk pool that lives on the device -- all pointers to device memory
  // are the same. we need this object to pass into device functions.
  t_chunk_pool *d_chunk_pool;

  t_chunk_pool( int chunk_pool_size, int chunk_size, int n_tils, int my_aid );
  void cleanup( void );

  // state sync functions -- both are called from host.
  __host__ void copy_state_to_device( void );
  __host__ void copy_state_from_device( void );

  __host__ int get_chunks_host( int nchunks, t_chunk **dest, t_mem_transfer_manager *mem_transfer_manager );
  __host__ int return_chunks_host( int nchunks, t_chunk **source );
  __device__ int get_chunks_device( int nchunks, t_chunk **dest );

  __host__ __device__ void init_vdfs( t_f1<p_k_part> *x_, t_f1<p_k_part> *p_,
                                                t_f0<p_k_part> *q_, t_f1<int> *ix_ );

  __host__ __device__ void set_data_pointers( p_k_part **x, p_k_part **p,
                                                          p_k_part **q, int **ix, t_chunk **chunk );

  __host__ __device__ void set_vdf_pointers( t_f1<p_k_part> *x_, t_f1<p_k_part> *p_,
                                                t_f0<p_k_part> *q_, t_f1<int> *ix_, t_chunk **chunk );
// #ifdef DEBUG_CHUNKS
//   __host__ void check_chunks( t_chunk **chunks, int nchunks );
// #endif


};

__global__ void init_chunks( t_chunk **chunk_queue, t_chunk *chunks,
                            int chunk_pool_size, size_t chunk_mem );

__host__ __device__ void unpack_chunk( p_k_part **x, p_k_part **p,
                                     p_k_part **q, int **ix, t_chunk **chunk, int chunk_size  );


// debug functions -- normally disabled 
#ifdef DEBUG_CHUNKS
// __host__ void check_chunks( t_chunk **chunks, int nchunks );
__global__ void test_chunk_pool_device( t_chunk_pool *chunk_pool );
void test_chunk_pool_host( t_chunk_pool *chunk_pool );
#endif

#endif
