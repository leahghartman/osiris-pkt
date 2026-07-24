#ifndef MEM_TRANSFER_MANAGER_H
#define MEM_TRANSFER_MANAGER_H

#include "os-param-cuda.h"
#include "os-sys-cuda.h"

// Author: Jacob Pierce, May 2022

#include <stdio.h>


// available copy modes
typedef enum {
  MEM_TRANSFER_MEMCPY,
  MEM_TRANSFER_CHUNKCOPY
} t_mem_transfer_mode;


// ********************************************
// bufferedCudaMemcpy
struct t_bufferedMemcpyHostMeta {
  void *src;
  void *dst;
  size_t count;
};

struct t_bufferedMemcpyDeviceMeta {
  void *src;
  void *dst;
  size_t count;
};

// ********************************************
// bufferedChunkCopy
struct t_bufferedChunkCopyHostMeta {
  p_k_part *x;
  p_k_part *p;
  p_k_part *q;
  int *ix;
  char *packed_data;
  int npar;
};

struct t_bufferedChunkCopyDeviceMeta {
  char *packed_data; // within d_com_buf
  t_chunk **chunks;
  int npar, npar_present;
};

struct t_bufferedChunkCopyGlobalMeta {
  int chunk_size;
};


// ********************************************
/* t_mem_transfer_manager
This class implements a function bufferedCudaMemcpy. This function stores the data and pointer addresses to be
copied, but only writes them if the internal communication buffer is filled. There is also a flush function to
manually flush the buffer when necessary. This function is not thread-safe (on the host) and does not work if
any of the pointer addresses change between the first load and the flush. This is agnostic to the specific
data represented in calls to memcpy.
*/
class t_mem_transfer_manager {

private :


  int nblocks;
  int block_size;
  int nblocks_requested;

  size_t data_buf_size; // max amount of data that can be held (both host and device)
  size_t hcopy_meta_buf_size;
  size_t dcopy_meta_buf_size;
  size_t global_meta_size;

  // these are char * so that we can use them for pointer arithmetic.
  char *h_data_buf;
  char *d_data_buf;

  char *h_hcopy_meta_buf;

  char *h_dcopy_meta_buf;
  char *d_dcopy_meta_buf;

  void *global_meta;

  char *h_data_cur;
  // char *h_hcopy_meta_cur;
  // char *h_dcopy_meta_cur;

  int ncopies_buffered;

  t_mem_transfer_mode copy_mode;
  cudaMemcpyKind kind;
  bool copy_not_append;

  void check_memcpy_kind( cudaMemcpyKind kind );

public :

  t_mem_transfer_manager( int data_buf_size_mb );
  void cleanup( void );

  // void realloc_hcopy_meta_buf( size_t size_new );
  // void realloc_dcopy_meta_bufs( size_t size_new );

  void checkCopyMode( t_mem_transfer_mode copy_mode_, cudaMemcpyKind kind_, int npar_present=p_none_present );

  void bufferedCudaMemcpy( void *dst, void *src, size_t count, cudaMemcpyKind kind_ );

  void bufferedChunkCopy( p_k_part *x, p_k_part *p, p_k_part *q, int *ix, t_chunk **chunks,
                          int npar, int npar_present, int chunk_size, cudaMemcpyKind kind_ );

  void flush( void );

  void flushBufferedCudaMemcpy( void );
  void flushBufferedChunkCopy( void );

};


__global__ void bufferedCudaMemcpyKernel( int *data_buf, t_bufferedMemcpyDeviceMeta *copy_meta_buf );
__global__ void bufferedChunkCopyKernel( char *data_buf, t_bufferedChunkCopyDeviceMeta *copy_meta_buf,
                                         const int chunk_size, cudaMemcpyKind kind );
__global__ void bufferedChunkAppendKernel( char *data_buf, t_bufferedChunkCopyDeviceMeta *copy_meta_buf,
                                         const int chunk_size, cudaMemcpyKind kind );




#endif // MEM_TRANSFER_MANAGER_H
