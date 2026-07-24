#include "os-mem-transfer-manager-cuda.h"
#include "os-sys-cuda.h"
#include "os-chunk-pool-cuda.h"
#include "os-param-cuda.h"

#include <iostream>
#include <assert.h>


// --------------------------------------------------------------------------------------
// Reallocation helper functions 
// --------------------------------------------------------------------------------------
void safeRealloc( void **ptr, const size_t size_new, size_t *size_old, int force ) {

  // if( ! force and *size_old >= size_new ) return;

  // allocate and point host meta bufs
  *ptr = realloc( *ptr, size_new );

  if( ! *ptr ) {
    fprintf(stderr, "ERROR: unable to realloc ptr.\n" );
    abort_program();
  }

  // *size_old = size_new;

}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void safeCudaRealloc( void **ptr, const size_t size_new, size_t *size_old, int force, int copy ) {

  // if( not force and *size_old >= size_new ) return;

  // realloc is not available in CUDA. do it manually here.
  void *ptr_tmp;
  chkErr(cudaMalloc(reinterpret_cast<void **>(&ptr_tmp), size_new));


  if( *ptr ) {
    if( copy ) cudaMemcpy( ptr_tmp, *ptr, *size_old, cudaMemcpyDeviceToDevice );
    cudaFree( *ptr );
  }

  *ptr = ptr_tmp;
  // *size_old = size_new;

}


/****************************** Class member functions *********************************/


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
t_mem_transfer_manager::t_mem_transfer_manager( int data_buf_size_mb ) {


  // default: 1 block per send, 1024 threads
  nblocks = -1;
  block_size = 1024;
  ncopies_buffered = 0;

  data_buf_size = static_cast<size_t>(data_buf_size_mb) * (1<<20);

  const size_t default_meta_buf_size = 1024 * sizeof(t_bufferedMemcpyHostMeta);

  // allocate h_data_buf
  h_data_buf = (char *) malloc( data_buf_size );

  if( ! h_data_buf ) {
    fprintf(stderr, "ERROR: unable to malloc host communication buffer.\n" );
    abort_program();
  }

  // allocate d_data_buf
  chkErrCstm(cudaMalloc(reinterpret_cast<void **>(&d_data_buf), data_buf_size),
      "ERROR: unable to cudaMalloc device communication buffer. Either decrease chunk_pool_size or com_buf_size.");

  // allocate meta bufs
  h_hcopy_meta_buf = NULL;
  h_dcopy_meta_buf = NULL;
  d_dcopy_meta_buf = NULL;
  global_meta = NULL;

  hcopy_meta_buf_size = default_meta_buf_size;
  dcopy_meta_buf_size = default_meta_buf_size;
  global_meta_size = sizeof(t_bufferedChunkCopyGlobalMeta);

  safeRealloc( (void **) &h_hcopy_meta_buf, default_meta_buf_size, &hcopy_meta_buf_size, true );
  safeRealloc( (void **) &h_dcopy_meta_buf, default_meta_buf_size, &dcopy_meta_buf_size, true );
  safeCudaRealloc( (void **) &d_dcopy_meta_buf, default_meta_buf_size, &dcopy_meta_buf_size, true, false );
  safeRealloc( (void **) &global_meta, sizeof(t_bufferedChunkCopyGlobalMeta), &global_meta_size, true );

  h_data_cur = h_data_buf;
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::cleanup( void ) {

  free( h_data_buf );
  cudaFree( d_data_buf );

  free( h_hcopy_meta_buf );
  free( h_dcopy_meta_buf );
  cudaFree( d_dcopy_meta_buf );

  free( global_meta );

}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::flush( void ) {

  switch( copy_mode ) {

    case( MEM_TRANSFER_MEMCPY ) :
      flushBufferedCudaMemcpy();

    case( MEM_TRANSFER_CHUNKCOPY ) :
      flushBufferedChunkCopy();

  }

}




// --------------------------------------------------------------------------------------
// set the data transfer kind if no data hash been loaded yet; otherwise,
// verify that the request was consistent with the kind initiated.
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::checkCopyMode( t_mem_transfer_mode copy_mode_, 
                                            cudaMemcpyKind kind_, int npar_present ) {

  // handle kind
  if( ncopies_buffered == 0 ) {

    copy_mode = copy_mode_;
    kind = kind_;
    copy_not_append = npar_present==p_none_present;

  } else {

    if( copy_mode != copy_mode_ ) {
      fprintf(stderr, "ERROR: requested bufferedCudaMemcpy copy_mode is inconsistent. Between flushes, all writes must be unidirectional.\n" );
      abort_program();
    }

    if( kind != kind_ ) {
      fprintf(stderr, "ERROR: requested bufferedCudaMemcpy kind is inconsistent. Between flushes, all kinds must be equivalent.\n" );
      abort_program();
    }

    if( copy_not_append != (npar_present==p_none_present) ) {
      fprintf(stderr, "ERROR: copy_not_append is inconsistent. Between flushes, we must have all copy or all append.\n" );
      abort_program();
    }

  }
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::bufferedCudaMemcpy( void *dst, void *src,
                                                   size_t count, cudaMemcpyKind kind_ ) {

  checkCopyMode( MEM_TRANSFER_MEMCPY, kind_ );

  // check that this is a valid length. it must be a multiple of sizeof(int) because the kernel
  // copies everything as ints.
  if( count % sizeof(int) != 0 ) {
    fprintf(stderr, "ERROR: invalid size request for bufferedCudaMemcpy.\n" );
    abort_program();
  }

  // check for flush. for device to device transfer, we assume the meta buf will never overflow and never
  // automatically flush. this is done for simplicity of implementation. if necessary, this could be rewritten
  // to pack both data and metadata in the same comm buff, or specifically to pack DtoD metadata in the data buf.
  // this would be more complex.
  if( kind == cudaMemcpyHostToDevice or kind == cudaMemcpyDeviceToHost ) {

    if( h_data_cur + count > h_data_buf + data_buf_size ) {

      if( ncopies_buffered == 0 ) {
        fprintf(stderr, "ERROR: data_buf_size is too small for single data transfer.\n" );
        abort_program();
      }

      flushBufferedCudaMemcpy();

    }

  }

  // realloc hcopy_meta if not big enough
  if( (ncopies_buffered + 1) * sizeof(t_bufferedMemcpyHostMeta) > hcopy_meta_buf_size ) {
    hcopy_meta_buf_size *= 2;
    safeRealloc( (void **) &h_hcopy_meta_buf, hcopy_meta_buf_size , &hcopy_meta_buf_size, true );
  }

  // realloc dcopy_meta bufs if not big enough
  if( (ncopies_buffered + 1) * sizeof(t_bufferedMemcpyDeviceMeta) > dcopy_meta_buf_size ) {
    size_t dcopy_meta_buf_size_old = dcopy_meta_buf_size;
    dcopy_meta_buf_size *= 2;
    safeRealloc( (void **) &h_dcopy_meta_buf, dcopy_meta_buf_size, &dcopy_meta_buf_size_old, true );
    safeCudaRealloc( (void **) &d_dcopy_meta_buf, dcopy_meta_buf_size, &dcopy_meta_buf_size_old, true, true );
  }

  // no global meta used for bufferedMemcpy, so no need to reallocate

  // prepare to set metadata
  t_bufferedMemcpyDeviceMeta *d_meta_tmp = (t_bufferedMemcpyDeviceMeta *) ( h_dcopy_meta_buf ) + ncopies_buffered;
  t_bufferedMemcpyHostMeta   *h_meta_tmp = (t_bufferedMemcpyHostMeta *)   ( h_hcopy_meta_buf ) + ncopies_buffered;

  if( kind == cudaMemcpyHostToDevice ) {

    // pack h_data_buf with data.
    memcpy( h_data_cur, src, count );

    // meta for device copy (currently lives on host)
    d_meta_tmp->src = d_data_buf + ( h_data_cur - h_data_buf );
    d_meta_tmp->dst = dst;
    d_meta_tmp->count = count;

    //  no need to track host meta after copy.
  }

  else if( kind == cudaMemcpyDeviceToHost ) {

    // unable to pack comm buf in this mode

    // meta for device copy (currently lives on host)
    d_meta_tmp->src = src;
    d_meta_tmp->dst = d_data_buf + ( h_data_cur - h_data_buf );
    d_meta_tmp->count = count;

    // track host meta
    h_meta_tmp->src = h_data_cur;
    h_meta_tmp->dst = dst;
    h_meta_tmp->count = count;
  }

  else if( kind == cudaMemcpyDeviceToDevice ) {

    // unable to pack comm buf in this mode
    // meta for device copy (currently lives on host)
    d_meta_tmp->src = src;
    d_meta_tmp->dst = dst ;
    d_meta_tmp->count = count;

    // no host meta in this mode
  }

  h_data_cur += count;
  ncopies_buffered++;

  return;
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::flushBufferedCudaMemcpy( void ) {

  const int nblocks_ = ncopies_buffered;

  const size_t total_meta_size = ncopies_buffered * sizeof(t_bufferedMemcpyDeviceMeta);
  const size_t total_data_size = h_data_cur - h_data_buf;

  if( ncopies_buffered == 0 ) return;

  if( kind == cudaMemcpyHostToDevice ) {

    // com buf is already packed. send to device

    // copy H to D
    chkErr( cudaMemcpy( d_data_buf, h_data_buf, total_data_size, cudaMemcpyHostToDevice ));

    chkErr( cudaMemcpy( d_dcopy_meta_buf, h_dcopy_meta_buf, total_meta_size, cudaMemcpyHostToDevice ));

    bufferedCudaMemcpyKernel<<<nblocks_, block_size>>>( (int *) d_data_buf,
                                                        (t_bufferedMemcpyDeviceMeta *) d_dcopy_meta_buf );
    chkLastErr("bufferedCudaMemcpyKernel");
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif

  }

  else if( kind == cudaMemcpyDeviceToHost ) {

    chkErr( cudaMemcpy( d_dcopy_meta_buf, h_dcopy_meta_buf, total_meta_size, cudaMemcpyHostToDevice ));

    bufferedCudaMemcpyKernel<<<nblocks_, block_size>>>( (int *) d_data_buf,
                                                        (t_bufferedMemcpyDeviceMeta *) d_dcopy_meta_buf );
    chkLastErr("bufferedCudaMemcpyKernel");
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif

    // copy D to H
    chkErr( cudaMemcpy( h_data_buf, d_data_buf, total_data_size, cudaMemcpyDeviceToHost ));

    // finally copy out to H dsts
    t_bufferedMemcpyHostMeta *meta;
    for( int i=0; i<ncopies_buffered; i++ ) {
      meta = (t_bufferedMemcpyHostMeta *)( h_hcopy_meta_buf ) + i;
      memcpy( meta->dst, meta->src, meta->count );
    }

  }

  else if( kind == cudaMemcpyDeviceToDevice ) {

    chkErr( cudaMemcpy( d_dcopy_meta_buf, h_dcopy_meta_buf, total_meta_size, cudaMemcpyHostToDevice ));

    bufferedCudaMemcpyKernel<<<nblocks_, block_size>>>( NULL, (t_bufferedMemcpyDeviceMeta *) d_dcopy_meta_buf );
    chkLastErr("bufferedCudaMemcpyKernel");
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif

  }

  // reset state
  ncopies_buffered = 0;
  h_data_cur = h_data_buf;
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::bufferedChunkCopy( p_k_part *x, p_k_part *p, p_k_part *q,
                                                  int *ix, t_chunk **chunks, int npar,
                                                  int npar_present, int chunk_size,
                                                  cudaMemcpyKind kind_ ) {

  checkCopyMode( MEM_TRANSFER_CHUNKCOPY, kind_, npar_present );

  size_t mem_size_x = npar * sizeof(p_k_part) * p_x_dim;
  size_t mem_size_p = npar * sizeof(p_k_part) * p_p_dim;
  size_t mem_size_q = npar * sizeof(p_k_part) * 1;
  size_t mem_size_ix = npar * sizeof(int) * p_x_dim;

  size_t mem_size_total = mem_size_x + mem_size_p + mem_size_q + mem_size_ix;

  // Ensure mem_size_total is a multiple sizeof(p_k_part) so that alignment requirements
  // are satisfied on the device.  We only have to do this at the end because the ints
  // (possibly of a different size) are packed at the end.
  mem_size_total = ( (mem_size_total+sizeof(p_k_part)-1) / sizeof(p_k_part) ) * sizeof(p_k_part);

  // flush if this would overflow the buffer
  if( h_data_cur + mem_size_total > h_data_buf + data_buf_size ) {

    if( ncopies_buffered == 0 ) {
      char data_buf_size_str[64];
      strmem( data_buf_size_str, data_buf_size );
      fprintf(stderr, "ERROR: com_buf_size is too small for single data transfer. Try increasing it. Current size is %s\n", data_buf_size_str);
      abort_program();
    }

    flushBufferedChunkCopy();
  }

  // realloc hcopy_meta if not big enough
  if( (ncopies_buffered + 1) * sizeof(t_bufferedMemcpyHostMeta) > hcopy_meta_buf_size ) {
    hcopy_meta_buf_size *= 2;
    safeRealloc( (void **) &h_hcopy_meta_buf, hcopy_meta_buf_size, &hcopy_meta_buf_size, false );
  }

  // realloc dcopy_meta bufs if not big enough
  if( (ncopies_buffered + 1) * sizeof(t_bufferedMemcpyDeviceMeta) > dcopy_meta_buf_size ) {
    size_t dcopy_meta_buf_size_old = dcopy_meta_buf_size;
    dcopy_meta_buf_size *= 2;
    safeRealloc( (void **) &h_dcopy_meta_buf, dcopy_meta_buf_size, &dcopy_meta_buf_size_old, false );
    safeCudaRealloc( (void **) &d_dcopy_meta_buf, dcopy_meta_buf_size, &dcopy_meta_buf_size_old, true, true );
  }

  // realloc global meta
  if( sizeof(t_bufferedChunkCopyGlobalMeta) > global_meta_size ) {
    global_meta_size = sizeof(t_bufferedChunkCopyGlobalMeta);
    safeRealloc( (void **) &h_hcopy_meta_buf, global_meta_size, &global_meta_size, false );
  }

  // prepare to set metadata
  t_bufferedChunkCopyDeviceMeta *d_meta_tmp = (t_bufferedChunkCopyDeviceMeta *) ( h_dcopy_meta_buf ) + ncopies_buffered;
  t_bufferedChunkCopyHostMeta   *h_meta_tmp = (t_bufferedChunkCopyHostMeta *)   ( h_hcopy_meta_buf ) + ncopies_buffered;
  t_bufferedChunkCopyGlobalMeta *global_meta_tmp = (t_bufferedChunkCopyGlobalMeta *) ( global_meta );

  global_meta_tmp->chunk_size = chunk_size;

  if( kind == cudaMemcpyHostToDevice ) {

    char *dst_tmp = h_data_cur;

    // pack h_data_buf
    memcpy( dst_tmp, x, mem_size_x );

    dst_tmp += mem_size_x;
    memcpy( dst_tmp, p, mem_size_p );

    dst_tmp += mem_size_p;
    memcpy( dst_tmp, q, mem_size_q );

    dst_tmp += mem_size_q;
    memcpy( dst_tmp, ix, mem_size_ix );

    //  no need to track host meta after copy.

    // meta for device copy (currently lives on host)
    d_meta_tmp->packed_data = d_data_buf + ( h_data_cur - h_data_buf );
    d_meta_tmp->chunks = chunks;
    d_meta_tmp->npar = npar;
    d_meta_tmp->npar_present = npar_present;

  }

  else{

    // unable to pack comm buf in this mode

    // meta for device copy (currently lives on host)
    d_meta_tmp->packed_data = d_data_buf + ( h_data_cur - h_data_buf );
    d_meta_tmp->chunks = chunks;
    d_meta_tmp->npar = npar;
    d_meta_tmp->npar_present = npar_present; // never actually used

    // track host meta
    h_meta_tmp->x = x;
    h_meta_tmp->p = p;
    h_meta_tmp->q = q;
    h_meta_tmp->ix = ix;
    h_meta_tmp->packed_data = h_data_cur;
    h_meta_tmp->npar = npar;

  }

  h_data_cur += mem_size_total;
  ncopies_buffered++;

  return;
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_mem_transfer_manager::flushBufferedChunkCopy() {

  const int nblocks_ = ncopies_buffered;

  const size_t total_meta_size = ncopies_buffered * sizeof(t_bufferedChunkCopyDeviceMeta);
  const size_t total_data_size = h_data_cur - h_data_buf;

  t_bufferedChunkCopyGlobalMeta *global_meta_tmp = (t_bufferedChunkCopyGlobalMeta *) ( global_meta );

  if( ncopies_buffered == 0 ) return;

  if( kind == cudaMemcpyHostToDevice ) {

    // com buf is already packed. send to device

    // copy H to D
    chkErr( cudaMemcpy( d_data_buf, h_data_buf, total_data_size, cudaMemcpyHostToDevice ));

    chkErr( cudaMemcpy( d_dcopy_meta_buf, h_dcopy_meta_buf, total_meta_size,
                        cudaMemcpyHostToDevice ));

    // check if we are appending or copying
    if ( copy_not_append ) {

      bufferedChunkCopyKernel<<<nblocks_, block_size>>>( d_data_buf,
                                       (t_bufferedChunkCopyDeviceMeta *) d_dcopy_meta_buf,
                                       global_meta_tmp->chunk_size, kind );
      chkLastErr("bufferedChunkCopyKernel");
#ifdef __DEBUG__
      chkErr(cudaDeviceSynchronize());
#endif

    } else {

      bufferedChunkAppendKernel<<<nblocks_, block_size>>>( d_data_buf,
                                       (t_bufferedChunkCopyDeviceMeta *) d_dcopy_meta_buf,
                                       global_meta_tmp->chunk_size, kind );
      chkLastErr("bufferedChunkAppendKernel");
#ifdef __DEBUG__
      chkErr(cudaDeviceSynchronize());
#endif

    }
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif
  }

  else if( kind == cudaMemcpyDeviceToHost ) {

    chkErr( cudaMemcpy( d_dcopy_meta_buf, h_dcopy_meta_buf, total_meta_size,
                        cudaMemcpyHostToDevice ));

    bufferedChunkCopyKernel<<<nblocks_, block_size>>>( d_data_buf,
                                       (t_bufferedChunkCopyDeviceMeta *) d_dcopy_meta_buf,
                                       global_meta_tmp->chunk_size, kind );
    chkLastErr("bufferedChunkCopyKernel");
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif

    // copy D to H
    chkErr( cudaMemcpy( h_data_buf, d_data_buf, total_data_size, cudaMemcpyDeviceToHost ));

    // finally copy out to H dsts
    t_bufferedChunkCopyHostMeta *meta;

    for( int i=0; i<ncopies_buffered; i++ ) {

      meta = (t_bufferedChunkCopyHostMeta *)( h_hcopy_meta_buf ) + i;

      size_t mem_size_x = meta->npar * sizeof(p_k_part) * p_x_dim;
      size_t mem_size_p = meta->npar * sizeof(p_k_part) * p_p_dim;
      size_t mem_size_q = meta->npar * sizeof(p_k_part) * 1;
      size_t mem_size_ix = meta->npar * sizeof(int) * p_x_dim;

      memcpy( meta->x, meta->packed_data, mem_size_x );

      meta->packed_data += mem_size_x;
      memcpy( meta->p, meta->packed_data, mem_size_p );

      meta->packed_data += mem_size_p;
      memcpy( meta->q, meta->packed_data, mem_size_q );

      meta->packed_data += mem_size_q;
      memcpy( meta->ix, meta->packed_data, mem_size_ix );

    }

  }

  // reset state
  ncopies_buffered = 0;
  h_data_cur = h_data_buf;
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
__global__ void bufferedCudaMemcpyKernel( int *data_buf, t_bufferedMemcpyDeviceMeta *copy_meta_buf ) {

  int i = blockIdx.x;

  t_bufferedMemcpyDeviceMeta meta = copy_meta_buf[i];

  // cast these as int * so that they copy faster than char
  int *dst = (int *) ( meta.dst );
  const int *src =  (int *) ( meta.src );
  const size_t count = meta.count / sizeof(int);


  for( int j=threadIdx.x; j<count; j+=blockDim.x ) {
  dst[j] = src[j];
  }

}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
__global__ void bufferedChunkCopyKernel( char *data_buf, t_bufferedChunkCopyDeviceMeta *copy_meta_buf,
                                         const int chunk_size, cudaMemcpyKind kind ) {

  int i = blockIdx.x;

  t_bufferedChunkCopyDeviceMeta meta = copy_meta_buf[i];

  char *packed_data = meta.packed_data;
  t_chunk **chunks = meta.chunks;
  int npar = meta.npar;

  int nchunks = CEIL_DIV( npar, chunk_size );

  int num_par_copy;

  p_k_part *x_chunk, *x_buf;
  p_k_part *p_chunk, *p_buf;
  p_k_part *q_chunk, *q_buf;
  int *ix_chunk, *ix_buf;

  if( kind == cudaMemcpyHostToDevice ) {

    unpack_chunk( &x_buf, &p_buf, &q_buf, &ix_buf, (t_chunk **)(&packed_data), npar );

    const p_k_part *x_buf_ = x_buf;
    const p_k_part *p_buf_ = p_buf;
    const p_k_part *q_buf_ = q_buf;
    const int *ix_buf_ = ix_buf;

    for( int i=0; i<nchunks; i++ ) {

      num_par_copy = min( chunk_size, npar - chunk_size * i );

      unpack_chunk( &x_chunk, &p_chunk, &q_chunk, &ix_chunk, chunks + i, chunk_size );

      for( int j=threadIdx.x; j<num_par_copy; j += blockDim.x ) {

        for( int k=0; k<p_x_dim; k++ ) {
           x_chunk[ k * chunk_size + j ] =  x_buf_[ ( i * chunk_size + j ) * p_x_dim + k ];
          ix_chunk[ k * chunk_size + j ] = ix_buf_[ ( i * chunk_size + j ) * p_x_dim + k ];
        }

        for( int k=0; k<p_p_dim; k++ ) {
          p_chunk[ k * chunk_size + j ] = p_buf_[ ( i * chunk_size + j ) * p_p_dim + k ];
        }

        q_chunk[ j ] = q_buf_[ i * chunk_size + j ];

      } // j

    } // i

  // device to host
  } else {

    for( int i=0; i<nchunks; i++ ) {

      num_par_copy = min( chunk_size, npar - chunk_size * i );

      unpack_chunk( &x_chunk, &p_chunk, &q_chunk, &ix_chunk, chunks + i, chunk_size );
      unpack_chunk( &x_buf, &p_buf, &q_buf, &ix_buf, (t_chunk **)(&packed_data), npar );

      const p_k_part *x_chunk_ = x_chunk;
      const p_k_part *p_chunk_ = p_chunk;
      const p_k_part *q_chunk_ = q_chunk;
      const int *ix_chunk_ = ix_chunk;

      for( int j=threadIdx.x; j<num_par_copy; j += blockDim.x ) {

        for( int k=0; k<p_x_dim; k++ ) {
           x_buf[ ( i * chunk_size + j ) * p_x_dim + k ] =  x_chunk_[ k * chunk_size + j ];
          ix_buf[ ( i * chunk_size + j ) * p_x_dim + k ] = ix_chunk_[ k * chunk_size + j ];
        }

        for( int k=0; k<p_p_dim; k++ ) {
          p_buf[ ( i * chunk_size + j ) * p_p_dim + k ] = p_chunk_[ k * chunk_size + j ];
        }

        q_buf[ i * chunk_size + j ] = q_chunk_[ j ];

      } // j

    } // i

  }

}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
__global__ void bufferedChunkAppendKernel( char *data_buf, t_bufferedChunkCopyDeviceMeta *copy_meta_buf,
                                           const int chunk_size, cudaMemcpyKind kind ) {

  int i = blockIdx.x;

  t_bufferedChunkCopyDeviceMeta meta = copy_meta_buf[i];

  char *packed_data = meta.packed_data;
  t_chunk **chunks = meta.chunks;
  int npar = meta.npar;
  int npar_present = meta.npar_present;

  int start_chunk = npar_present / chunk_size;
  int start_idx = npar_present % chunk_size;
  int nchunks = 1 + CEIL_DIV( npar - chunk_size + start_idx, chunk_size );

  int num_par_copy;

  p_k_part *x_chunk, *x_buf;
  p_k_part *p_chunk, *p_buf;
  p_k_part *q_chunk, *q_buf;
  int *ix_chunk, *ix_buf;

  if( kind == cudaMemcpyHostToDevice ) {

    unpack_chunk( &x_buf, &p_buf, &q_buf, &ix_buf, (t_chunk **)(&packed_data), npar );

    const p_k_part *x_buf_ = x_buf;
    const p_k_part *p_buf_ = p_buf;
    const p_k_part *q_buf_ = q_buf;
    const int *ix_buf_ = ix_buf;

    for( int i=0; i<nchunks; i++ ) {

      num_par_copy = min( chunk_size - start_idx*(i==0), npar - chunk_size * i + start_idx*(i>0) );

      unpack_chunk( &x_chunk, &p_chunk, &q_chunk, &ix_chunk, chunks + start_chunk + i, chunk_size );

      for( int j=threadIdx.x+start_idx*(i==0); j<num_par_copy+start_idx*(i==0); j += blockDim.x ) {


        for( int k=0; k<p_x_dim; k++ ) {
           x_chunk[ k * chunk_size + j ] =  x_buf_[ ( i * chunk_size - start_idx + j ) * p_x_dim + k ];
          ix_chunk[ k * chunk_size + j ] = ix_buf_[ ( i * chunk_size - start_idx + j ) * p_x_dim + k ];
        }

        for( int k=0; k<p_p_dim; k++ ) {
          p_chunk[ k * chunk_size + j ] = p_buf_[ ( i * chunk_size - start_idx + j ) * p_p_dim + k ];
        }

        q_chunk[ j ] = q_buf_[ i * chunk_size - start_idx + j ];

      } // j

    } // i

  // device to host
  } else {

    // has not yet been implemented, and should raise an error
    assert(0);

  }

}



// void test_mem_transfer_manager( void ) {



// }

