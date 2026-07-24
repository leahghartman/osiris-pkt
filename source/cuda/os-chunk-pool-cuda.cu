#include <stdio.h>
#include <cuda_runtime.h>

#include "os-chunk-pool-cuda.h"
#include "os-sys-cuda.h"



// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
t_chunk_pool::t_chunk_pool( int chunk_pool_size_, int chunk_size_, int n_tils_, int my_aid ){
                            

  chunk_pool_size = chunk_pool_size_;
  chunk_size      = chunk_size_;
  n_tils          = n_tils_;

  if ( chunk_size%__BLOCK_SIZE__ != 0 ) {
    fprintf(stderr, "ERROR: chunk size, %d, must be divisible by the block size, %d.\n",
                    chunk_size, __BLOCK_SIZE__);
    abort_program();
  }


  // offsets for where data is stored within a chunk, total chunk memory
  x_offset = 0;
  p_offset = x_offset + p_x_dim * sizeof( p_k_part ) * chunk_size;
  q_offset = p_offset +  p_p_dim * sizeof( p_k_part ) * chunk_size;
  ix_offset = q_offset + sizeof( p_k_part ) * chunk_size;

  chunk_mem = ix_offset + p_x_dim * sizeof( int ) * chunk_size;

  size_t free_mem, total_mem;
  chkErr( cudaMemGetInfo( &free_mem, &total_mem ) );

  // ideal max chunk pool size (in reality, max size is slightly less than this due to 
  // alignment issues, etc, and the fact that we allocate a couple things after this)
  int max_chunk_pool_size = free_mem / chunk_mem;

  // set chunk_pool_size to default value, 75% of available memory, should be relatively safe
  if ( chunk_pool_size == -1 ) {
    chunk_pool_size = 0.75 * max_chunk_pool_size;
  }

  size_t chunk_pool_mem = chunk_mem * chunk_pool_size;

  // subtle point -- this is one chunk larger than necessary so that we can
  // leave one space empty when it's filled at initialization. this makes it so that the queue
  // can never become full (if new chunks aren't added to it). this allows us to distinguish
  // between the states start = end as being empty or full: start = end is empty.
  chkErr(cudaMalloc(reinterpret_cast<void **>(&chunk_queue), sizeof(t_chunk *) * ( chunk_pool_size )  ));

  char chunk_mem_str[64];
  char chunk_pool_mem_str[64];
  char free_mem_str[64];
  char total_mem_str[64];

  strmem( chunk_mem_str, chunk_mem );
  strmem( chunk_pool_mem_str, chunk_pool_mem );
  strmem( free_mem_str, free_mem );
  strmem( total_mem_str, total_mem );

  if (my_aid==1) {
    printf( "   GPU memory allocation summary:\n"
            "   --- total_mem per GPU              = %s\n"
            "   --- free GPU mem before alloc      = %s\n"
            "   --- chunk size                     = %s (%d particles)\n"
            "   --- chunk_pool size                = %s (%d chunks)\n"
            "   --- max chunk_pool_size            = %d chunks\n",
            total_mem_str, 
            free_mem_str,
            chunk_mem_str, chunk_size, 
            chunk_pool_mem_str, chunk_pool_size,
            max_chunk_pool_size );
  }

  chkErrCstm(cudaMalloc(reinterpret_cast<void **>(&chunks), chunk_pool_mem ),
             "ERROR: chunk_pool_size exceeds available memory.");

  // now check how much mem we used
  chkErr( cudaMemGetInfo( &free_mem, &total_mem ) );

  strmem( free_mem_str, free_mem );

  if (my_aid==1) {
    printf( "   --- mem remaining after allocation = %s\n", free_mem_str );
  }

  // allocate host state pointers
  state_size = 2 * sizeof(int);
  h_state = (int *) malloc( state_size );
  h_start = h_state;
  h_end = h_state + 1;

  // allocate device state pointers
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_state), 3 * sizeof(int) ));
  d_start = d_state;
  d_end = d_state + 1;
  d_end_tmp = d_state + 2;

  *h_start = 0;
  *h_end = chunk_pool_size - 1;

  copy_state_to_device();

  // populate all chunks
  init_chunks <<< 1, 1024 >>> ( chunk_queue, chunks, chunk_pool_size, chunk_mem );
  chkLastErr("init_chunks");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

  // alloc a copy of ourself and put it on the device. we can then effectively feed ourself into kernels.
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_chunk_pool), sizeof(t_chunk_pool) ));
  chkErr(cudaMemcpy( d_chunk_pool, this, sizeof( t_chunk_pool ), cudaMemcpyHostToDevice));

#ifdef DEBUG_CHUNKS
  PRINT( "Testing chunk pool...");
 
  // uncomment to test get_chunk_pool_device
  // test_chunk_pool_device<<< 32, 32 >>> ( d_chunk_pool );
  // chkLastErr("test_chunk_pool_device");
  // chkErr(cudaDeviceSynchronize());

  test_chunk_pool_host( this );

  PRINT( "Success" );
  exit(-1);
#endif

}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void t_chunk_pool::cleanup( void ) {

  free( h_state );

  chkErr( cudaFree( d_state ) );

  chkErr( cudaFree( chunks ) );
  chkErr( cudaFree( chunk_queue) );

  chkErr( cudaFree( d_chunk_pool ) );

}




// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
__host__ void t_chunk_pool::copy_state_to_device( void ) {
    chkErr(cudaMemcpy( d_state, &h_state, state_size, cudaMemcpyHostToDevice));
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
__host__ void t_chunk_pool::copy_state_from_device( void ) {
    chkErr(cudaMemcpy( &h_state, d_state, state_size, cudaMemcpyDeviceToHost));
}


// --------------------------------------------------------------------------------------
// initialize vdfs with the write dimensions for referencing the data.
// --------------------------------------------------------------------------------------
__host__ __device__ void t_chunk_pool::init_vdfs( t_f1<p_k_part> *x, t_f1<p_k_part> *p,
                                                t_f0<p_k_part> *q, t_f1<int> *ix ) {

    x->init( chunk_size, p_x_dim, NULL, false );
    p->init( chunk_size, p_p_dim, NULL, false );
    q->init( chunk_size, NULL, false );
    ix->init( chunk_size, p_x_dim, NULL, false );

}


// --------------------------------------------------------------------------------------
// todo: think of better name for this function
// --------------------------------------------------------------------------------------
__host__ __device__ void t_chunk_pool::set_data_pointers( p_k_part **x, p_k_part **p,
                                                          p_k_part **q, int **ix, t_chunk **chunk ) {

  *x = (p_k_part *) ( *chunk + x_offset );
  *p = (p_k_part *) ( *chunk + p_offset );
  *q = (p_k_part *) ( *chunk + q_offset );
  *ix = (int *)     ( *chunk + ix_offset );

}

// --------------------------------------------------------------------------------------
// todo: think of better name for this function
// --------------------------------------------------------------------------------------
__host__ __device__ void unpack_chunk( p_k_part **x, p_k_part **p,
                                     p_k_part **q, int **ix, t_chunk **chunk, const int chunk_size  ) {

  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * sizeof( p_k_part ) * chunk_size;
  const size_t q_offset = p_offset +  p_p_dim * sizeof( p_k_part ) * chunk_size;
  const size_t ix_offset = q_offset + sizeof( p_k_part ) * chunk_size;

  *x = (p_k_part *) ( *chunk + x_offset );
  *p = (p_k_part *) ( *chunk + p_offset );
  *q = (p_k_part *) ( *chunk + q_offset );
  *ix = (int *)     ( *chunk + ix_offset );

}


// --------------------------------------------------------------------------------------
// todo: think of better name for this function
// --------------------------------------------------------------------------------------
__host__ __device__ void t_chunk_pool::set_vdf_pointers( t_f1<p_k_part> *x, t_f1<p_k_part> *p,
                                                t_f0<p_k_part> *q, t_f1<int> *ix, t_chunk **chunk ) {

  x->g = (p_k_part *) ( *chunk + x_offset );
  p->g = (p_k_part *) ( *chunk + p_offset );
  q->g = (p_k_part *) ( *chunk + q_offset );
  ix->g = (int *)     ( *chunk + ix_offset );

}


// __host__ t_chunk_pool::copy_chunk( t_chunk *dst, t_chunk *src,
//                                    cudaMemcpyKind kind, t_mem_transfer_manager *mem_transfer_manager ) {

//   mem_transfer_manager->bufferedCudaMemcpy( dst, src, chunk_mem, kind );

// }




// --------------------------------------------------------------------------------------
// debug functions
// --------------------------------------------------------------------------------------
#ifdef DEBUG_CHUNKS

// --------------------------------------------------------------------------------------
// test contention: many blocks and threads grab chunks from the pool simultaneously. should never fail. 
// --------------------------------------------------------------------------------------
__global__ void test_chunk_pool_device( t_chunk_pool *chunk_pool ) {

  t_chunk *chunks[2];

  int ret;

  while(1) {
    ret = chunk_pool->get_chunks_device( 2, chunks );
    // PRINT_INT( ret ) ;
    if( ! ret ) break;
  }

}

// --------------------------------------------------------------------------------------
// test that queue works for reading and writing past boundary 
// set chunk_pool_size to 8 in input deck for this test 
// --------------------------------------------------------------------------------------
void test_chunk_pool_host( t_chunk_pool *chunk_pool ) {

  const int test_buf_size = 8; 

  t_chunk **d_chunks;  // = (t_chunk **) malloc( sizeof(t_chunk *) * 8 );
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_chunks), sizeof(t_chunk *) * test_buf_size ) );

  t_chunk **h_chunks = (t_chunk **) malloc( sizeof(t_chunk *) * test_buf_size );

  int ret;

  // get 8 chunks 
  PRINT("Getting initial chunks...");
  ret = chunk_pool->get_chunks_host( 8, d_chunks, NULL );
  chkErr(cudaDeviceSynchronize());
  if(!ret) PRINT( "Serious problem." );

  // chkErr(cudaMemcpy( h_chunks, d_chunks, test_buf_size * sizeof( t_chunk *), cudaMemcpyDeviceToHost));
  // for( int i=0; i<test_buf_size; i++ ) PRINT_PTR( h_chunks[i] );

  // return 4 -- this does not write past end of queue, but appends to end. 
  PRINT("Returning initial chunks...");
  ret = chunk_pool->return_chunks_host( 4, d_chunks + 4 );
  chkErr(cudaDeviceSynchronize());
  if(!ret) PRINT( "Serious problem." );

  // chkErr(cudaMemcpy( h_chunks, d_chunks, test_buf_size * sizeof( t_chunk *), cudaMemcpyDeviceToHost));
  // for( int i=0; i<test_buf_size; i++ ) PRINT_PTR( h_chunks[i] );

  // now read 4 chunks. this will read past the end, since the start of the queue is at 8.
  PRINT("Getting chunks past end...");
  ret = chunk_pool->get_chunks_host( 4, d_chunks + 4, NULL );
  chkErr(cudaDeviceSynchronize());
  if(!ret) PRINT( "ERROR: failed to wrap around queue for get_chunks_host.");

  // chkErr(cudaMemcpy( h_chunks, d_chunks, test_buf_size * sizeof( t_chunk *), cudaMemcpyDeviceToHost));
  // for( int i=0; i<test_buf_size; i++ ) PRINT_PTR( h_chunks[i] );

  // now give all chunks back to host. this will wrap around.
  PRINT("Returning chunks past end...");
  ret = chunk_pool->return_chunks_host( 8, d_chunks );
  chkErr(cudaDeviceSynchronize());
  if(!ret) PRINT( "ERROR: failed to wrap around queue for return_chunks_host.");

  PRINT("Testing return_chunks error handling...");
  ret = chunk_pool->return_chunks_host( 8, d_chunks );
  chkErr(cudaDeviceSynchronize());
  if(ret) PRINT( "ERROR: error should have been returned.");

  // chkErr(cudaMemcpy( h_chunks, d_chunks, test_buf_size * sizeof( t_chunk *), cudaMemcpyDeviceToHost));
  // for( int i=0; i<test_buf_size; i++ ) PRINT_PTR( h_chunks[i] );

  PRINT("Getting all chunks again...");
  ret = chunk_pool->get_chunks_host( 8, d_chunks, NULL );
  chkErr(cudaDeviceSynchronize());
  if(!ret) PRINT( "ERROR: failed to get chunks again.");

  PRINT("Testing get_chunks error handling...");
  ret = chunk_pool->get_chunks_host( 8, d_chunks, NULL );
  chkErr(cudaDeviceSynchronize());
  if(ret) PRINT( "ERROR: error should have been returned.");

  // chkErr(cudaMemcpy( h_chunks, d_chunks, test_buf_size * sizeof( t_chunk *), cudaMemcpyDeviceToHost));
  // for( int i=0; i<test_buf_size; i++ ) PRINT_PTR( h_chunks[i] );

  chkErr( cudaFree( d_chunks ) );
  free( h_chunks );

}

#endif




// --------------------------------------------------------------------------------------
// serial version
// --------------------------------------------------------------------------------------
__host__ int t_chunk_pool::get_chunks_host( int nchunks, t_chunk **dest,
                                            t_mem_transfer_manager *mem_transfer_manager ) {

  int nchunks_available;

  if( nchunks == 0 ) return 1;

  // compute num available chunks in the chunk pool
  if( *h_start <= *h_end ) {
    nchunks_available = *h_end - *h_start;
  } else {
    nchunks_available = chunk_pool_size + *h_end - *h_start;
  }

  // fail if we are requesting more chunks than available
  if( nchunks_available <= nchunks ){
    fprintf( stderr, "ERROR: unable to get chunks -- not enough available in the pool.\n" );
    return 0;
  }

  int start_old = *h_start;

  *h_start = (*h_start + nchunks) % chunk_pool_size;

  // use buffered mem copy if desired by calling function
  if( mem_transfer_manager ) {

   if( start_old < *h_start ) {
      mem_transfer_manager->bufferedCudaMemcpy( dest, chunk_queue + start_old, nchunks * sizeof( t_chunk *), cudaMemcpyDeviceToDevice);
    } else {
      int nchunks1 = chunk_pool_size - start_old;
      int nchunks2 = nchunks - nchunks1;
      mem_transfer_manager->bufferedCudaMemcpy( dest, chunk_queue + start_old, nchunks1 * sizeof( t_chunk *), cudaMemcpyDeviceToDevice);
      mem_transfer_manager->bufferedCudaMemcpy( dest + nchunks1, chunk_queue, nchunks2 * sizeof( t_chunk *), cudaMemcpyDeviceToDevice);
    }

  }

  // else use regular memcpy
  else {

   if( start_old < *h_start ) {
      chkErr(cudaMemcpy( dest, chunk_queue + start_old, nchunks * sizeof( t_chunk *),
                         cudaMemcpyDeviceToDevice));
    } else {
      int nchunks1 = chunk_pool_size - start_old;
      int nchunks2 = nchunks - nchunks1;
      chkErr(cudaMemcpy( dest, chunk_queue + start_old, nchunks1 * sizeof( t_chunk *),
                         cudaMemcpyDeviceToDevice));
      chkErr(cudaMemcpy( dest + nchunks1, chunk_queue, nchunks2 * sizeof( t_chunk *),
                         cudaMemcpyDeviceToDevice));
    }

  }

  // // debug
  // check_chunks( dest, nchunks );

  return 1;

}


// --------------------------------------------------------------------------------------
// serial version
// --------------------------------------------------------------------------------------
__host__ int t_chunk_pool::return_chunks_host( int nchunks, t_chunk **source ) {

  int nchunks_available;

  if( nchunks == 0 ) return 1;

  if( *h_start <= *h_end ) {
    nchunks_available = *h_end - *h_start;
  } else {
    nchunks_available = chunk_pool_size + *h_end - *h_start;
  }

  // Sanity check: this condition should be impossible to be reached under normal use
  if( nchunks >= chunk_pool_size - nchunks_available ){
    fprintf( stderr, "ERROR: unable to return chunks -- not enough space available in pool.\n" );
    return 0;
  }

  int end_old = *h_end;

  *h_end = ( *h_end + nchunks ) % ( chunk_pool_size );

  if( end_old < *h_end ) {
    chkErr(cudaMemcpy( chunk_queue + end_old, source, nchunks * sizeof( t_chunk * ),
                       cudaMemcpyDeviceToDevice));
  } else {
    int nchunks1 = chunk_pool_size - end_old;
    int nchunks2 = nchunks - nchunks1;
    chkErr(cudaMemcpy( chunk_queue + end_old, source, nchunks1 * sizeof( t_chunk * ),
                       cudaMemcpyDeviceToDevice));
    chkErr(cudaMemcpy( chunk_queue, source + nchunks1, nchunks2 * sizeof( t_chunk * ),
                       cudaMemcpyDeviceToDevice));
  }

  return 1;

}


// --------------------------------------------------------------------------------------
// get chunks in parallel from the device
// --------------------------------------------------------------------------------------
__device__ int t_chunk_pool::get_chunks_device( int nchunks, t_chunk **dest )
{
  int nchunks_available;
  int start_old_tmp, start_new, end_old;

  if( nchunks == 0 ) return 1;

  int start_old = * d_start;

  while( 1 ) {

    // double check that enough chunks are available
    end_old = *d_end;

    if( start_old <= end_old ) {
     nchunks_available = end_old - start_old;
    } else {
     nchunks_available = chunk_pool_size + end_old - start_old;
    }

    if( nchunks_available <= nchunks ) return 0;

    start_old_tmp = start_old;

    start_new = (start_old + nchunks) % chunk_pool_size;

    start_old = atomicCAS( d_start, start_old, start_new );

    if( start_old == start_old_tmp ) break;
  }

  // printf( "start_old=%d, start_new=%d\n", start_old, start_new );

  int i;

  if( start_old < start_new ) {

    for( i=start_old; i<start_new; i++ ) {
      dest[i] = chunk_queue[i];
    }

  }
  else {

    for( i=start_old; i<chunk_pool_size; i++ ) {
      dest[i] = chunk_queue[i];
    }

    for( i=0; i<start_new; i++ ) {
      dest[i] = chunk_queue[i];
    }

  }

 return 1;
}




// // this function will not be implemented unless necessary

// __device__ int t_chunk_pool::return_chunks_device( int nchunks, t_chunk *source ) {

//  int end_tmp_cur_cur, nchunks_available;

//  while( 1 ) {
//    end_tmp_cur_cur = end_tmp_cur;
//    if( __sync_val_compare_and_swap( &end_tmp_cur, end_tmp_cur_cur, (end_tmp_cur_cur + nchunks)%chunk_pool_size ) != end_tmp_cur ) break;
//  }

//  int start_old = h_start;

//  if( start_old < end_tmp_cur_cur ) {
//    nchunks_available = h_end - h_start;
//  } else {
//    nchunks_available = chunk_pool_size + h_end - h_start;
//  }

//  // this condition should be impossible to be reached under normal use
//  if( nchunks > chunk_pool_size - nchunks_available ) return 0;

//  for( int i=0; i<nchunks; i++ ) {
//    chunk_idx_queue[ (end_tmp_cur_cur + i)%chunk_pool_size ] = source[i].idx;
//  }

//  // only one guy updates the h_end at a time
//  // do not update until the actual queue h_end has caught up to where we inserted data
//  while( 1 ) {
//    // int h_end = h_end;
//    if( __sync_val_compare_and_swap( &h_end, end_tmp_cur_cur, (end_tmp_cur_cur + nchunks)%chunk_pool_size ) != h_end ) break;
//  }

//  return 1;
// }




// #ifdef DEBUG_CHUNKS
// __host__ void t_chunk_pool::check_chunks( t_chunk **chunks, int nchunks ) {

//  size_t mem = nchunks * sizeof( t_chunk *);

//  t_chunk **h_chunks = (t_chunk **) malloc( mem );

//     chkErr( cudaMemcpy( h_chunks, chunks, mem, cudaMemcpyDeviceToHost ) );

//     int status = 1;

//     for( int i=0; i<nchunks; i++ ) {
//      if( h_chunks[i]->debug_num != CHUNK_POOL_DEBUG_NUM ) {
//        status = 0;
//      }
//      // PRINT_INT( h_chunks[i].debug_num );
//     }

//     if( ! status ) {
//      PRINT( "ERROR: chunks are invalid." );
//     } else {
//      PRINT("Chunks are valid.");
//     }

//     free( h_chunks );

// }
// #endif




// --------------------------------------------------------------------------------------
// init chunk queue
// --------------------------------------------------------------------------------------
__global__ void init_chunks( t_chunk **chunk_queue, t_chunk *chunks,
                            int chunk_pool_size, size_t chunk_mem ) {

  for( int i=threadIdx.x; i<chunk_pool_size; i += blockDim.x ) {

    chunk_queue[i] = chunks + i * chunk_mem;

// #ifdef DEBUG_CHUNKS
//     chunk_queue[i].debug_num = CHUNK_POOL_DEBUG_NUM;
// #endif
  }

}




