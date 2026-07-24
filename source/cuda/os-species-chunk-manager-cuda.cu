#include "os-species-chunk-manager-cuda.h"


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
t_species_chunk_manager::t_species_chunk_manager( void ) {

  x = NULL;
  p = NULL;
  q = NULL;
  ix = NULL;

  h_chunks = NULL;
  d_chunks = NULL;

  h_nchunks_max = NULL;
  d_nchunks_max = NULL;

  h_nchunks = NULL;
  d_nchunks = NULL;

  h_state_changed = false;

  nchunks_pad = 0;

  chunk_growth_factor = 0;

  chunk_pool = NULL;
  mem_transfer_manager = NULL;

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::init_group( t_f1<p_k_part> *x_, t_f1<p_k_part> *p_,
                                          t_f0<p_k_part> *q_, t_f1<int> *ix_,
                                          int n_tils, int nchunks_pad_, int chunk_growth_factor_,
                                          t_chunk_pool *chunk_pool_,
                                          t_mem_transfer_manager *mem_transfer_manager_ ){

  x = x_;
  p = p_;
  q = q_;
  ix = ix_;

  nchunks_pad = nchunks_pad_;
  chunk_growth_factor = chunk_growth_factor_;
  chunk_pool = chunk_pool_;
  mem_transfer_manager = mem_transfer_manager_;

  // prepare the allocations
  h_chunks     = (t_chunk ***) malloc( n_tils * sizeof( t_chunk **) );
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_chunks), n_tils * sizeof(t_chunk **) ));

  h_nchunks_max     = (int *) malloc( n_tils * sizeof(int) );
  h_nchunks         = (int *) malloc( n_tils * sizeof(int) );
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_nchunks_max), n_tils * sizeof(int) ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_nchunks), n_tils * sizeof(int) ));

  h_state_changed = false;

}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::init_tile( int til, int nchunks_max_init, int n_tils, 
                                         bool if_last_til ){


  // h_chunks is a host array storing device pointers
  chkErr( cudaMalloc( reinterpret_cast<void **>(&(h_chunks[til])),
                                        nchunks_max_init * sizeof(t_chunk *) ) );
  h_nchunks_max[til] = nchunks_max_init;
  h_nchunks[til]     = 0;

}


// ---------------------------------------------------------------------------------------
// copy to device once we've initialized the value for all tiles
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::init_tile_last( int n_tils ){

  // d_chunks is a device pointer storing device pointers
  chkErr( cudaMemcpy(d_chunks, h_chunks, n_tils * sizeof(t_chunk **), cudaMemcpyHostToDevice) );

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::init_tile_dlb( int til, t_chunk **d_chunk, int nchunks_init,
                                             int nchunks_max_init ){


  h_chunks[til]      = d_chunk;
  h_nchunks_max[til] = nchunks_max_init;
  h_nchunks[til]     = nchunks_init;

  // Note: copying of h_chunks to device is done later once we've recd our new tiles

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::cleanup( int n_tils ) {

  for(int i=0; i<n_tils; i++) {
    cudaFree( h_chunks[i] );
  }

  free( h_chunks );

  free( h_nchunks_max );
  free( h_nchunks );

  cudaFree( d_chunks );
  cudaFree( d_nchunks_max );
  cudaFree( d_nchunks );

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::cleanup_dlb( int *tils_send_old, int n_send ) {

  int til;

  for(int i=0; i<n_send; i++) {
    til = tils_send_old[i];
    cudaFree( h_chunks[til] );
  }

  free( h_chunks );

  free( h_nchunks_max );
  free( h_nchunks );

  cudaFree( d_chunks );
  cudaFree( d_nchunks_max );
  cudaFree( d_nchunks );

}


// ---------------------------------------------------------------------------------------
// snag some mo chunks from the queueueue
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::mo_chunks( int til, int nchunks_mo ) {

  // new num chunks this tile will hold
  int nchunks_new = nchunks_mo + h_nchunks[til];

  // first grow the chunk buffer if we need more space
  if( nchunks_new > h_nchunks_max[til] ) {

    // new size of tchunk ** for this tile
    int nchunks_max_new = int( chunk_growth_factor * nchunks_new ) + 1;

    // t_chunk **h_chunks_til_old = hh_chunks[ til ];
    t_chunk **d_chunks_til_old = h_chunks[ til ];

    // allocate new chunk array of appropriate size
    // t_chunk **h_chunks_til_new;
    t_chunk **d_chunks_til_new;

    // h_chunks_til_new = (t_chunk **) malloc( nchunks_max_new * sizeof(t_chunk *));
    chkErr( cudaMalloc(reinterpret_cast<void **>(&d_chunks_til_new), nchunks_max_new * sizeof(t_chunk *)) );

    // copy values of old chunk array to new chunk array (both are on the device)
    // this memcpy can't be buffered unless the subsequent cudaFree is also buffered. this could
    // be done but would require nontrivial modification of the memoryTransferManager

    chkErr( cudaMemcpy( d_chunks_til_new, d_chunks_til_old, h_nchunks[til] * sizeof(t_chunk *), cudaMemcpyDeviceToDevice ) );

    // mem_transfer_manager->bufferedCudaMemcpy( d_chunks_til_new, d_chunks_til_old, h_nchunks[til] * sizeof(t_chunk *), cudaMemcpyDeviceToDevice );

    chkErr(cudaFree( d_chunks_til_old ));

    h_chunks[til] = d_chunks_til_new;

    h_nchunks_max[til] = nchunks_max_new;

  }

  // now there's enough space, grab the chunks

  // get reference to chunk array for this tile on host
  t_chunk **d_chunks_til = h_chunks[ til ];

  // add new chunks onto end of array
  int status = chunk_pool->get_chunks_host( nchunks_mo, d_chunks_til + h_nchunks[til],
                                            mem_transfer_manager );

  h_nchunks[til] = nchunks_new;

  // set flag so that we know to write to the device
  h_state_changed = true;

  return status;

}


// ---------------------------------------------------------------------------------------
// return nchunks_blow chunks to the queue
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::blow_chunks( int til, int nchunks_blow ) {

  // new num chunks this tile will hold
  int nchunks_new = h_nchunks[til] - nchunks_blow;

  // get reference to chunk array for this tile on host
  t_chunk **d_chunks_til = h_chunks[ til ];

  // return chunks from the end of the array to the pool
  int status = chunk_pool->return_chunks_host( nchunks_blow, d_chunks_til + nchunks_new );

  h_nchunks[til] = nchunks_new;

  // set flag so that we know to write to the device
  h_state_changed = true;

  return status;

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_chunk_manager::sync_state( int n_tils, cudaMemcpyKind kind ) {

  int *nchunks_src, *nchunks_dst;
  int *nchunks_max_src, *nchunks_max_dst;

  if( kind == cudaMemcpyHostToDevice ) {

    if( ! h_state_changed ) return;

    nchunks_src = h_nchunks;
    nchunks_dst = d_nchunks;

    nchunks_max_src = h_nchunks_max;
    nchunks_max_dst = d_nchunks_max;

  } else {

    if( h_state_changed ) {
      fprintf( stderr, "ERROR: attempted to sync chunk_manager state from device to host when host state had changed." );
      exit(-1);
    }

    nchunks_src = d_nchunks;
    nchunks_dst = h_nchunks;

    nchunks_max_src = d_nchunks_max;
    nchunks_max_dst = h_nchunks_max;

  }

  // this can only happen in one direction. only applies when cudaMalloc is called to
  // reallocate elements of h_chunks on the host.
  if( kind == cudaMemcpyHostToDevice )
    mem_transfer_manager->bufferedCudaMemcpy( d_chunks, h_chunks, n_tils * sizeof(t_chunk **),
                                              kind );

  mem_transfer_manager->bufferedCudaMemcpy( nchunks_dst, nchunks_src, n_tils * sizeof(int),
                                            kind );
  mem_transfer_manager->bufferedCudaMemcpy( nchunks_max_dst, nchunks_max_src, n_tils * sizeof(int),
                                            kind );

  mem_transfer_manager->flush();

  h_state_changed = false;

}


// ---------------------------------------------------------------------------------------
// Calculate how many chunks **would** be necessary to make sure each tile has exactly 
// nchunks_pad number of empty chunks. Don't actually get the chunks here though.
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::prune_and_pad_chunks_dry( int n_tils, int *num_par ) {

  int chunk_size = chunk_pool->chunk_size;
  int chunk_pool_size = chunk_pool->chunk_pool_size;

  // compute num chunks...
  int nchunks_held = 0; // currently held
  int nchunks_full_req = 0; // that will be required to hold all incoming particles
  int nchunks_pad_req = n_tils * nchunks_pad; // needed to store pad chunks

  for (int til=0; til<n_tils; ++til) {
    // calculate num chunks neeeded to hold all particles
    // by full here we really just mean not empty
    int nchunks_full = num_par[til] / chunk_size;
    if( num_par[til] % chunk_size > 0 ) nchunks_full ++;
    nchunks_full_req += nchunks_full;

    nchunks_held += h_nchunks[til];
  }

  // compute num available chunks in the chunk pool
  int nchunks_available;
  if( *(chunk_pool->h_start) <= *(chunk_pool->h_end) ) {
    nchunks_available = *(chunk_pool->h_end) - *(chunk_pool->h_start);
  } else {
    nchunks_available = chunk_pool->chunk_pool_size + *(chunk_pool->h_end) - *(chunk_pool->h_start);
  }

  if ( nchunks_full_req+nchunks_pad_req > chunk_pool_size ){
    fprintf(stderr, "ERROR: Not enough chunks available in the pool. "
           "num tiles: %d. "
           "num chunks required to store all particles: %d. "
           "n_chunks_pad: %d. "
           "chunk_pool_size: %d. "
           "nchunks available: %d. "
           "nchunks held: %d\n",
           n_tils, nchunks_full_req, nchunks_pad, chunk_pool_size, nchunks_available,
           nchunks_held);
    return 0;
  }

  return 1;

}


// ---------------------------------------------------------------------------------------
// Make sure each tile has exactly nchunks_pad number of empty chunks
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::prune_and_pad_chunks( int n_tils, int *num_par ) {

  int status;

  // Check to see if there are enough chunks to give every tile it's share
  status = prune_and_pad_chunks_dry( n_tils, num_par );
  if (!status) return status;

  for (int til=0; til<n_tils; ++til) {
    status = prune_and_pad_chunks_til( til, num_par[til] );
    if (!status) return status;
  }

  // we need to flush the mem transfer manager after calling mo_chunks()
  mem_transfer_manager->flush();

  sync_state( n_tils, cudaMemcpyHostToDevice );

  return 1;

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::prune_and_pad_chunks_til( int til, int num_par ){

  int status = 1;
  int chunk_size = chunk_pool->chunk_size;

  // calculate num chunks neeeded to hold all particles
  // by full here we really just mean not empty
  int nchunks_full = num_par / chunk_size;
  if( num_par % chunk_size > 0 ) nchunks_full ++;

  int nchunks_pad_actual = h_nchunks[til] - nchunks_full;

  if(  nchunks_pad_actual < nchunks_pad ) {
    // get more empty chunks if we don't have enough
    int nchunks_mo = nchunks_pad - nchunks_pad_actual;
    if( nchunks_mo <= 0 ) return 1;
    status = mo_chunks( til, nchunks_mo );
  } else if ( nchunks_pad_actual > nchunks_pad ) {
    // return empty chunks if we had too many
    int nchunks_blow = nchunks_pad_actual - nchunks_pad;
    if( nchunks_blow <= 0 ) return 1;
    status = blow_chunks( til, nchunks_blow );
  }

  return status;

}


// ---------------------------------------------------------------------------------------
// this function requires that enough chunks are already available for copying
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::copy_particles_to_device( int til, int num_par_til ) {

  int chunk_size = chunk_pool->chunk_size;

  if (num_par_til < 1) return 1;

  int nchunks_copy = CEIL_DIV(num_par_til,chunk_size);

  // grab more chunks since we didn't have space before
  if( nchunks_copy > h_nchunks[til] ) {
    fprintf( stderr, "ERROR: attempted to copy particles to device when not enough "
                     "chunks were available.\n");
    return 0;
  }


  t_chunk **d_chunks_til = h_chunks[til];

  p_k_part *h_x_tmp = x[til].g;
  p_k_part *h_p_tmp = p[til].g;
  p_k_part *h_q_tmp = q[til].g;
  int *h_ix_tmp     = ix[til].g;

  mem_transfer_manager->bufferedChunkCopy( h_x_tmp, h_p_tmp, h_q_tmp, h_ix_tmp,
                                           d_chunks_til, num_par_til, p_none_present,
                                           chunk_pool->chunk_size, cudaMemcpyHostToDevice );

  return 1;

}


// ---------------------------------------------------------------------------------------
// this function requires that enough chunks are already available for copying
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::append_particles_to_device( int til, int num_par_til,
                                                            int num_par_present ) {

  int chunk_size = chunk_pool->chunk_size;

  if (num_par_til+num_par_present < 1) return 1;

  int nchunks_copy = CEIL_DIV( num_par_til+num_par_present, chunk_size );

  // grab more chunks since we didn't have space before
  if( nchunks_copy > h_nchunks[til] ) {
    fprintf( stderr, "ERROR: attempted to append particles to device when not enough chunks were available.\n");
    return 0;
  }


  t_chunk **d_chunks_til = h_chunks[til];

  p_k_part *h_x_tmp = x[til].g;
  p_k_part *h_p_tmp = p[til].g;
  p_k_part *h_q_tmp = q[til].g;
  int *h_ix_tmp     = ix[til].g;

  mem_transfer_manager->bufferedChunkCopy( h_x_tmp, h_p_tmp, h_q_tmp, h_ix_tmp,
                                           d_chunks_til, num_par_til, num_par_present,
                                           chunk_pool->chunk_size, cudaMemcpyHostToDevice );

  return 1;

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
int t_species_chunk_manager::copy_particles_from_device( int til, int num_par_til,
                                                          int num_par_max_til ) {

  if ( num_par_til > num_par_max_til ) {
    fprintf(stderr,
        "ERROR: Attempted to copy particles from device when CPU buffers were too small: "
        "num_par = %d greater than num_par_max = %d\n",num_par_til, num_par_max_til );
    return 0;
  }

  if (num_par_til < 1) return 1;

  t_chunk **d_chunks_til = h_chunks[til];

  p_k_part *h_x_tmp = x[til].g;
  p_k_part *h_p_tmp = p[til].g;
  p_k_part *h_q_tmp = q[til].g;
  int *h_ix_tmp =     ix[til].g;

  mem_transfer_manager->bufferedChunkCopy( h_x_tmp, h_p_tmp, h_q_tmp, h_ix_tmp,
                                           d_chunks_til, num_par_til, p_none_present,
                                           chunk_pool->chunk_size, cudaMemcpyDeviceToHost );

  return 1;

}
