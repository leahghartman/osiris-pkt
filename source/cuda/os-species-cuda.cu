#include <stdio.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <cstring>
using namespace std;

#include "os-species-cuda.h"
#include "os-spec-push-cuda.h"
#include "os-sort-cuda.h"
#include "os-spec-validate-cuda.h"


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
t_species_c::t_species_c( void ) {
  n_tils = 0;
  rqm = 0;
  free_stream = false;
  for (int i=0; i<p_max_dim; ++i) {
    dx[i] = 0;
    move_num[i] = 0;
  }
  interpolation = 0;
  grid_center = false;
  dt = 0;
  d_my_nx_p = NULL;

  num_par_max = NULL;
  num_par = NULL;
  d_num_par = NULL;
  energy_fortran = NULL;
  energy = NULL;
  d_energy = NULL;

  x = NULL;
  p = NULL;
  q = NULL;
  ix = NULL;

  d_n_hole = NULL;
  d_n_recv_loc = NULL;
  d_n_comm_ext = NULL;

  next = NULL;

  chunk_pool = NULL;

  nblocks = 0;

  d_block_start_tiles = NULL;
  d_block_start_chunks = NULL;
  d_block_num_par = NULL;

  // populated in the input deck
  ihole_size = 0;
  nchunks_pad = 0;
  nblocks_requested = 0;
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::init_group( int n_tils_, t_chunk_pool *chunk_pool_,
                              t_mem_transfer_manager *mem_transfer_manager_,
                              p_k_part rqm_, bool free_stream_, p_double *dx_,
                              int interpolation_, bool grid_center_, p_double dt_,
                              int ihole_size_, p_double ihole_size_frac, int nchunks_pad_,
                              p_double nchunks_pad_frac, p_double chunk_growth_factor_,
                              int nblocks_requested_, int tiles_per_node ) {

  n_tils      = n_tils_;
  rqm         = rqm_;
  free_stream = free_stream_;

  for (int i=0; i<p_max_dim; ++i)
    dx[i] = dx_[i];

  interpolation        = interpolation_;
  grid_center          = grid_center_;
  dt                   = dt_;
  chunk_pool           = chunk_pool_;
  mem_transfer_manager = mem_transfer_manager_;
  nblocks_requested    = nblocks_requested_;

  // Can set ihole_size as a fraction of (a proxy for) num particles/tile.
  // In practice this number will be much larger (~2-3x) than the actual number of particles
  // per tiles, so ihole_size_frac should correspondingly be lower than you'd normally think.
  if ( ihole_size_ == -1 ) {
    int num_par_per_tile = chunk_pool->chunk_pool_size * chunk_pool->chunk_size / tiles_per_node;
    ihole_size = static_cast<int>(ihole_size_frac * num_par_per_tile);
  } else {
    ihole_size = ihole_size_;
  }

  if ( nchunks_pad_ == -1 ) {
    // Can set nchunks_pad as a fraction of the number of chunks/tile
    nchunks_pad = static_cast<int>(nchunks_pad_frac * chunk_pool->chunk_pool_size / tiles_per_node);
  } else {
    nchunks_pad = nchunks_pad_;
  }

  // Set up n_tils arrays of specified lengths that contain hole indices
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_ihole), n_tils * sizeof(t_f0<int>) ));
  t_f0<int> h_single_ihole = t_f0<int>( ihole_size, false );
  for( int i=0; i<n_tils; i++ ) {
    h_single_ihole.malloc_copy_to_device( d_ihole + i );
  }
  h_single_ihole.cleanup();

  my_nx_p = t_f1<int>( p_x_dim, n_tils, false );

  num_par_max    = new int[n_tils];
  num_par        = new int[n_tils];
  energy_fortran = new p_double*[n_tils];
  energy         = new p_double[n_tils];

  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_energy), n_tils * sizeof(p_double) ));

  x = (t_f1<p_k_part>*) malloc(n_tils * sizeof (t_f1<p_k_part>));
  p = (t_f1<p_k_part>*) malloc(n_tils * sizeof (t_f1<p_k_part>));
  q = (t_f0<p_k_part>*) malloc(n_tils * sizeof (t_f0<p_k_part>));
  ix= (t_f1<int>*)      malloc(n_tils * sizeof (t_f1<int>));

  // initialize chunk managers: chunk_manager for regular chunks (holds most of data)
  // and chunk_manager for external communication (typically just a few chunks per tile)
  chunk_manager = new t_species_chunk_manager;
  chunk_manager->init_group( x, p, q, ix, n_tils, nchunks_pad, chunk_growth_factor_,
                             chunk_pool, mem_transfer_manager );

  chunk_manager_ext = new t_species_chunk_manager;
  chunk_manager_ext->init_group( x, p, q, ix, n_tils, nchunks_pad, chunk_growth_factor_,
                                 chunk_pool, mem_transfer_manager  );

  // Set up d_n_hole, d_n_recv_loc, d_n_comm_ext
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_n_hole), n_tils * sizeof(int) ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_n_recv_loc), n_tils * sizeof(int) ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_n_comm_ext), n_tils * sizeof(int) ));

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::init_tile( int til, int num_par_max_, int num_par_, p_double *energy_,
                             int *my_nx_p_, p_k_part *x_, p_k_part *p_, p_k_part *q_,
                             int *ix_, int nchunks_max_init_, bool if_last_til ) {

  num_par_max[til]    = num_par_max_;
  num_par[til]        = num_par_;
  energy_fortran[til] = energy_;

  for( int i=0; i<p_x_dim; i++ ) {
    // only need the number of cells
    my_nx_p.fh(i,til) = my_nx_p_[3*i+2];
  }

  new (x+til) t_f1<p_k_part>( p_x_dim, num_par_max_, x_, false );
  new (p+til) t_f1<p_k_part>( p_p_dim, num_par_max_, p_, false );
  new (q+til) t_f0<p_k_part>( num_par_max_, q_, false );
  new (ix+til)t_f1<int>     ( p_x_dim, num_par_max_, ix_, false );

  // init chunk manager, but we don't grab chunks until later once we know the num particles
  chunk_manager->init_tile( til, nchunks_max_init_, n_tils, if_last_til );

  // init ext chunk manager and grab chunks since we already know how many we want
  chunk_manager_ext->init_tile( til, nchunks_max_init_, n_tils, if_last_til  );
  int status = chunk_manager_ext->prune_and_pad_chunks_til(til, 0);
  if (!status) abort_program();

  if( til == n_tils-1 or if_last_til ) {
    // flush after calling prune and pad chunks
    mem_transfer_manager->flush();

    init_tile_last();
  }

}


// ---------------------------------------------------------------------------------------
// Complete initialization for tiles
// ---------------------------------------------------------------------------------------
void t_species_c::init_tile_last(){

  // sync state after calling prune and pad chunks
  chunk_manager_ext->sync_state( n_tils, cudaMemcpyHostToDevice );

  chkErr(cudaMalloc( reinterpret_cast<void **>(&d_my_nx_p), my_nx_p.size()*sizeof(int) ));
  chkErr(cudaMemcpy( d_my_nx_p, my_nx_p.g, my_nx_p.size()*sizeof(int), cudaMemcpyHostToDevice ));

  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_num_par), sizeof(int)*n_tils));
  chkErr(cudaMemcpy(d_num_par, num_par, sizeof(int)*n_tils, cudaMemcpyHostToDevice));

  chunk_manager -> init_tile_last( n_tils );
  chunk_manager_ext -> init_tile_last( n_tils );

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::init_group_dlb( int n_tils_new, t_species_c *spec,
				  int ihole_size_, p_double ihole_size_frac,
				  int nchunks_pad_, p_double nchunks_pad_frac ) {

  n_tils         = n_tils_new;
  rqm            = spec->rqm;
  free_stream    = spec->free_stream;

  for (int i=0; i<p_max_dim; ++i)
    dx[i] = spec->dx[i];

  interpolation        = spec->interpolation;
  grid_center          = spec->grid_center;
  dt                   = spec->dt;
  chunk_pool           = spec->chunk_pool;
  mem_transfer_manager = spec->mem_transfer_manager;
  nblocks_requested    = spec->nblocks_requested;

  if ( ihole_size_ == -1 ) {
    int num_par_per_tile = chunk_pool->chunk_pool_size * chunk_pool->chunk_size / n_tils;
    ihole_size = static_cast<int>(ihole_size_frac * num_par_per_tile);
  } else {
    ihole_size = ihole_size_;
  }

  // compute new chunks_pad if necessary
  if ( nchunks_pad_ == -1 ) {
    // Can set nchunks_pad as a fraction of the number of chunks/tile
    nchunks_pad = static_cast<int>(nchunks_pad_frac * chunk_pool->chunk_pool_size / n_tils );
  } else {
    nchunks_pad = nchunks_pad_;
  }
  
  // Set up n_tils_new arrays of specified lengths that contain hole indices
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_ihole), n_tils_new * sizeof(t_f0<int>) ));
  t_f0<int> h_single_ihole = t_f0<int>( ihole_size, false );
  for( int i=0; i<n_tils_new; i++ ) {
    h_single_ihole.malloc_copy_to_device( d_ihole + i );
  }
  h_single_ihole.cleanup();

  my_nx_p = t_f1<int>( p_x_dim, n_tils_new, false );

  num_par_max    = new int[n_tils_new];
  num_par        = new int[n_tils_new];
  energy_fortran = new p_double*[n_tils_new];
  energy         = new p_double[n_tils_new];

  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_energy), n_tils_new * sizeof(p_double) ));

  x = (t_f1<p_k_part>*) malloc(n_tils_new * sizeof (t_f1<p_k_part>));
  p = (t_f1<p_k_part>*) malloc(n_tils_new * sizeof (t_f1<p_k_part>));
  q = (t_f0<p_k_part>*) malloc(n_tils_new * sizeof (t_f0<p_k_part>));
  ix= (t_f1<int>*)      malloc(n_tils_new * sizeof (t_f1<int>));

  // initialize chunk managers: chunk_manager for regular chunks (holds most of data)
  // and chunk_manager for external communication (typically just a few chunks per tile)
  
  p_double chunk_growth_factor_ = spec->chunk_manager->chunk_growth_factor;
  t_chunk_pool *chunk_pool_     = spec->chunk_manager->chunk_pool;

  t_mem_transfer_manager *mem_transfer_manager_ = spec->chunk_manager->mem_transfer_manager;

  chunk_manager = new t_species_chunk_manager;
  chunk_manager->init_group( x, p, q, ix, n_tils_new, nchunks_pad, chunk_growth_factor_,
                             chunk_pool_, mem_transfer_manager_ );

  chunk_manager_ext = new t_species_chunk_manager;
  chunk_manager_ext->init_group( x, p, q, ix, n_tils_new, nchunks_pad, chunk_growth_factor_,
                                 chunk_pool_, mem_transfer_manager_  );

  // we are going to have to re-sync with device after we've recieved any tiles from dlb
  chunk_manager -> h_state_changed = true;
  chunk_manager_ext -> h_state_changed = true;

  // Set up d_n_hole, d_n_recv_loc, d_n_comm_ext
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_n_hole), n_tils_new * sizeof(int) ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_n_recv_loc), n_tils_new * sizeof(int) ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_n_comm_ext), n_tils_new * sizeof(int) ));

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::init_tile_dlb( int til, int til_old, t_species_c *spec ) {

  int nchunks_init;
  int nchunks_max_init;
  t_chunk **d_chunk;

  num_par_max[til]    = spec->num_par_max[til_old];
  num_par[til]        = spec->num_par[til_old];
  energy_fortran[til] = spec->energy_fortran[til_old];

  for( int i=0; i<p_x_dim; i++ ) {
    // only need the number of cells
    my_nx_p.fh(i,til) = spec->my_nx_p.fh(i,til_old);
  }

  x[til] = spec -> x[til_old];
  p[til] = spec -> p[til_old];
  q[til] = spec -> q[til_old];
  ix[til] = spec -> ix[til_old];

  // init chunk manager
  nchunks_init     = spec->chunk_manager->h_nchunks[til_old];
  nchunks_max_init = spec->chunk_manager->h_nchunks_max[til_old];
  d_chunk          = spec->chunk_manager->h_chunks[til_old];
  chunk_manager->init_tile_dlb( til, d_chunk, nchunks_init, nchunks_max_init );

  // init ext chunk manager
  nchunks_init     = spec->chunk_manager_ext->h_nchunks[til_old];
  nchunks_max_init = spec->chunk_manager_ext->h_nchunks_max[til_old];
  d_chunk          = spec->chunk_manager_ext->h_chunks[til_old];
  chunk_manager_ext->init_tile_dlb( til, d_chunk, nchunks_init, nchunks_max_init );
  int status = chunk_manager_ext->prune_and_pad_chunks_til(til, 0);
  if (!status) abort_program();

  // Note: final initialization steps are done (ie init_tile_last is called) either after
  // recieving tiles, if any, otherwise at the very end of dlb

}


// ---------------------------------------------------------------------------------------
// Update C as to the locations of the new particle buffers
// Called everytime particle buffers are resized (read: reallocated) in fortran
// ---------------------------------------------------------------------------------------
void t_species_c::sync_buffer( int til, int num_par_max_, p_k_part *x_, p_k_part *p_,
                               p_k_part *q_, int *ix_ ) {

  num_par_max[til] = num_par_max_;

  (x+til)->init( p_x_dim, num_par_max_, x_, false );
  (p+til)->init( p_p_dim, num_par_max_, p_, false );
  (q+til)->init( num_par_max_, q_, false );
  (ix+til)->init( p_x_dim, num_par_max_, ix_, false );

}


#if 1 
// The old get_block_partition prioritized launching an ideal number of blocks
// regardless of the number of tiles on that node, specified by the input deck
// parameter n_blocks_requested. 
// In this version, we just want one block per tile, which is ideal for shared memory.


// ---------------------------------------------------------------------------------------
// Computes parameters for particle processing kernels (nblocks, d_block_num_par,
// d_block_start_chunks, d_block_start_tiles)
// ---------------------------------------------------------------------------------------
void t_species_c::get_block_partition( int chunk_size ){

  nblocks = n_tils;

  t_f0<int> h_block_start_tiles = t_f0<int>( nblocks, false );
  t_f0<int> h_block_start_chunks = t_f0<int>( nblocks, false );
  t_f0<int> h_block_num_par = t_f0<int>( nblocks, false );

  // set starting tile, and chunk on that tile, for each block
  for( int til=0; til<n_tils; ++til) {
    h_block_start_tiles.fh(til) = til;
    h_block_start_chunks.fh(til) = 0;
    h_block_num_par.fh(til) = num_par[til];
  }

// print block partition
#if 0
  for( int i=0; i<nblocks; i++ ) {
    printf( "block %d: %d %d %d \n", i, h_block_start_tiles.fh(i), h_block_start_chunks.fh(i), h_block_num_par.fh(i) );
  }
#endif

  // copy to device
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_block_start_tiles), sizeof(int)*nblocks ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_block_start_chunks), sizeof(int)*nblocks ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_block_num_par), sizeof(int)*nblocks ));

  chkErr(cudaMemcpy(d_block_start_tiles, h_block_start_tiles.f, sizeof(int)*nblocks, cudaMemcpyHostToDevice));
  chkErr(cudaMemcpy(d_block_start_chunks, h_block_start_chunks.f, sizeof(int)*nblocks, cudaMemcpyHostToDevice));
  chkErr(cudaMemcpy(d_block_num_par, h_block_num_par.f, sizeof(int)*nblocks, cudaMemcpyHostToDevice));

  // Clean up memory
  h_block_start_tiles.cleanup();
  h_block_start_chunks.cleanup();
  h_block_num_par.cleanup();


}

#else

// ---------------------------------------------------------------------------------------
// Computes parameters for particle processing kernels (nblocks, d_block_num_par,
// d_block_start_chunks, d_block_start_tiles)
// ---------------------------------------------------------------------------------------
void t_species_c::get_block_partition( int chunk_size ){

  int tot_nchunks = 0;
  for( int til=0; til<n_tils; til++) {
    tot_nchunks += (num_par[til]+chunk_size-1) / chunk_size;
  }

  // use one block per chunk if there aren't enough chunks (there should be
  // more chunks than blocks in ordinary usage)
  int nchunks_per_block = max( tot_nchunks / nblocks_requested, 1 );

  // ensure there's at least one block per tile
  int nblocks_max = 0;
  for( int til=0; til<n_tils; til++) {
    int nchunks_on_this_tile = (num_par[til]+chunk_size-1) / chunk_size;
    nblocks_max += max( nchunks_on_this_tile / nchunks_per_block, 1 );
  }

  t_f0<int> h_block_start_tiles = t_f0<int>( nblocks_max, false );
  t_f0<int> h_block_start_chunks = t_f0<int>( nblocks_max, false );
  t_f0<int> h_block_num_par = t_f0<int>( nblocks_max, false );

  nblocks = 0;

  // TODO: Make this a kernel
  for( int til=0; til<n_tils; ++til) {

    // calculate number of chunks on this tile
    int nchunks_on_this_tile = (num_par[til]+chunk_size-1) / chunk_size;

    // set the number of blocks on this tile
    // make sure there's at least one block per tile
    int nblocks_on_this_tile = max( nchunks_on_this_tile / nchunks_per_block, 1 );

    for( int j = 0; j<nblocks_on_this_tile; j++ ) {

      // set starting tile, and chunk on that tile, for each block
      h_block_start_tiles.fh(nblocks) = til;
      h_block_start_chunks.fh(nblocks) = nchunks_per_block * j;

      // set the number of particles to be processed by each block, h_block_num_par
      // these blocks have only completely full chunks
      if( j < nblocks_on_this_tile - 1 ) {

        h_block_num_par.fh(nblocks) = chunk_size * nchunks_per_block;

      // the last chunk on the last block might not be completely full of particles
      } else {

        int num_par_remaining = num_par[til] - ( nblocks_on_this_tile - 1 ) * chunk_size * nchunks_per_block;

        h_block_num_par.fh(nblocks) = num_par_remaining;

      }

      nblocks++;

    }
  }

// print block partition
#if 0
  for( int i=0; i<nblocks; i++ ) {
    printf( "block %d: %d %d %d \n", i, h_block_start_tiles.fh(i), h_block_start_chunks.fh(i), h_block_num_par.fh(i) );
  }

  printf( "nblocks actual / requested / max: %d / %d / %d\n", nblocks, nblocks_requested, nblocks_max );
#endif

  // copy to device
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_block_start_tiles), sizeof(int)*nblocks_max ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_block_start_chunks), sizeof(int)*nblocks_max ));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_block_num_par), sizeof(int)*nblocks_max ));

  chkErr(cudaMemcpy(d_block_start_tiles, h_block_start_tiles.f, sizeof(int)*nblocks_max, cudaMemcpyHostToDevice));
  chkErr(cudaMemcpy(d_block_start_chunks, h_block_start_chunks.f, sizeof(int)*nblocks_max, cudaMemcpyHostToDevice));
  chkErr(cudaMemcpy(d_block_num_par, h_block_num_par.f, sizeof(int)*nblocks_max, cudaMemcpyHostToDevice));

  // Clean up memory
  h_block_start_tiles.cleanup();
  h_block_start_chunks.cleanup();
  h_block_num_par.cleanup();


}

#endif


// ---------------------------------------------------------------------------------------
// Move window for particles, then receive any additional particles from the CPU
// ---------------------------------------------------------------------------------------
void t_species_c::move_window( int *nmove, int *num_par_, t_node_conf_c *no_co ){

  move_num[0] = nmove[0];
  move_num[1] = nmove[1];
  move_num[2] = nmove[2];

  // Setup execution parameters
  int chunk_size = chunk_pool->chunk_size;

  get_block_partition( chunk_size );

// This is redundant if particles are validated after update_boundary_2. But I'll keep it here
// just in case you want to use it for an extra sanity check
// #ifdef __DEBUG__
//   validate( no_co, false, "before move window" );
// #endif

  // shift particle positions
  shift_positions_wrapper( chunk_manager->d_chunks, d_block_start_tiles,
                           d_block_start_chunks, d_block_num_par,
                           nmove[0], nmove[1], nmove[2], nblocks, chunk_size);

  // copy particles from host to device from moving window injection
  for (int til=0; til<n_tils; ++til) {
    int status = chunk_manager->append_particles_to_device( til, num_par_[til], num_par[til] );
    if (!status) abort_program();

    num_par[til] += num_par_[til]; // update num_par for this species
    num_par_[til] = num_par[til];  // update array that CPU has a handle on with new num_par
  }
  mem_transfer_manager->flush();

  // update device num_par buffer
  chkErr(cudaMemcpy(d_num_par, num_par_, sizeof(int)*n_tils, cudaMemcpyHostToDevice));

// This is redundant if particles are validated before advance_deposit. But I'll keep it here
// just in case you want to use it for an extra sanity check
// #ifdef __DEBUG__
//   // first, need to recalculate block partition since number of particles has changed
//   chkErr(cudaFree(d_block_start_tiles));
//   chkErr(cudaFree(d_block_start_chunks));
//   chkErr(cudaFree(d_block_num_par));

//   get_block_partition( chunk_size );

//   validate( no_co, false, "after move window" );
// #endif


  // Clean up memory for block partition
  chkErr(cudaFree(d_block_start_tiles));
  chkErr(cudaFree(d_block_start_chunks));
  chkErr(cudaFree(d_block_num_par));

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::pre_advance_deposit( t_node_conf_c *no_co, t_current_c *jay, t_emf_c *emf ){

  int chunk_size = chunk_pool->chunk_size;

  get_block_partition( chunk_size );

  for( int i=0; i<n_tils; i++) {
    energy[i] = 0;
  }
  if (!free_stream) {
    chkErr(cudaMemcpy( d_energy, energy, sizeof(p_double)*n_tils, cudaMemcpyHostToDevice));
  }

  // Wait for all streams to finish
  chkErr(cudaStreamSynchronize( no_co->s_j ));

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::advance_deposit( t_node_conf_c *no_co, t_current_c *jay, t_emf_c *emf,
                                   bool report_energy ){

  int chunk_size = chunk_pool->chunk_size;

  // calculate amount of shared memory to allocate in advance_deposit/dudt
  // this is equal to the memory size of the largest field grid for a tile
  size_t shmem_size_adv_dep, shmem_size_dudt;
  shmem_size_adv_dep = sizeof(p_k_fld) * jay->nx_tot_max;
  shmem_size_dudt    = sizeof(p_k_fld) * emf->nx_tot_max * 2; // half for e, half for b

  // // ***************************
  // // To disable shared memory in dudt, uncomment this, and comment/uncomment
  // // the appropriate DUDT_BORIS_*D in os-spec-push-cuda.cu
  // // ***************************
  // shmem_size_dudt = 0;
  // // ***************************

  // check if tiles are too big to fit into shared memory for this device
  if ( (shmem_size_dudt > no_co->maxSharedMemoryPerBlockOptin) or 
             shmem_size_adv_dep > no_co->maxSharedMemoryPerBlockOptin){
    // printf("   Shmem per block requested for dudt:    %zu bytes\n"
    //        "   Shmem per block requested for adv dep: %zu bytes\n",
    //        shmem_size_dudt, shmem_size_adv_dep);
    // printf("   Max shared memory per block:            %zu bytes\n"
    //        "   Max shared memory per block (Opt in):   %zu bytes\n",
    //         no_co->maxSharedMemoryPerBlock,
    //         no_co->maxSharedMemoryPerBlockOptin);
    fprintf(stderr,"ERROR: Too much shared memory per block requested. "
                   "Try decreasing tile size.\n");
    abort_program();
  }

#ifdef __DEBUG__
  validate( no_co, false, "before advance deposit" );
#endif

  // Call kernel stub - Update momenta
  // energy is currently broken (should be by tile)
  if ( ! free_stream ){
    dudt_boris_wrapper( chunk_manager->d_chunks, d_block_start_tiles, d_block_start_chunks,
                        d_block_num_par, rqm, emf, interpolation, grid_center, d_energy,
                        report_energy, dt, nblocks, chunk_size,
                        no_co->s_adv_dep, shmem_size_dudt );
  }

  // Call kernel stub - advance position of particles and deposit current
  advance_deposit_wrapper( chunk_manager->d_chunks, d_block_start_tiles, d_block_start_chunks,
                           d_block_num_par, dx, jay, interpolation, dt, nblocks,
                           chunk_size, d_n_hole, d_n_comm_ext,
                           d_n_recv_loc, no_co->s_adv_dep, shmem_size_adv_dep );

#ifdef __DEBUG__
  validate( no_co, true, "after advance deposit" );
#endif

}


// ---------------------------------------------------------------------------------------
// Sort particles and gather information on holes in particle buffer
// ---------------------------------------------------------------------------------------
void t_species_c::update_boundary_1( t_node_conf_c *no_co, int *num_par_ ){

  int chunk_size = chunk_pool->chunk_size;

  // make sure each tile has some empty chunks
  int status = chunk_manager->prune_and_pad_chunks( n_tils, num_par );
  if (!status) abort_program();

  // Call kernel stub - sort particles
  sort_wrapper( chunk_manager->d_chunks, chunk_manager_ext->d_chunks,
                chunk_manager->d_nchunks_max, chunk_manager->d_nchunks,
                chunk_manager_ext->d_nchunks_max, chunk_manager_ext->d_nchunks,
                d_my_nx_p, d_block_start_tiles, d_block_start_chunks,
                d_block_num_par, chunk_size, d_ihole, ihole_size, d_n_hole,
                no_co->d_neighbor_til_id, no_co->d_shift, d_n_comm_ext, d_n_recv_loc,
                d_num_par, nblocks, no_co->d_err );

  // check error code from kernel
  check_sort_error( no_co, "sort" );

  // update buffer containing how many particles each tile needs to send
  chkErr(cudaMemcpy(num_par_, d_n_comm_ext, sizeof(int)*n_tils, cudaMemcpyDeviceToHost));

}


// ---------------------------------------------------------------------------------------
// copy particles from device to host for MPI comms and physical boundaries
// ---------------------------------------------------------------------------------------
void t_species_c::copy_particles_ext_from_device( int *num_par_ ){

  for (int til=0; til<n_tils; ++til) {
    int status = chunk_manager_ext->copy_particles_from_device( til, num_par_[til], num_par_max[til] );
    if (!status) abort_program();
  }
  mem_transfer_manager->flush();

}


// ---------------------------------------------------------------------------------------
// Fill holes in each tile's particle buffer
// ---------------------------------------------------------------------------------------
void t_species_c::update_boundary_2( t_node_conf_c *no_co, int *num_par_ ){

  int chunk_size = chunk_pool->chunk_size;
  int nblocks = n_tils;

  // here, num_par_ contains the number of particles received from other nodes
  chkErr(cudaMemcpy(d_n_comm_ext, num_par_, sizeof(int)*n_tils, cudaMemcpyHostToDevice));

  // copy particles from host to device from MPI comms and physical boundaries
  copy_particles_ext_to_device( num_par_ );

  // make sure there are enough chunks to store new number of particles
  pad_chunks_ub2( num_par_ );

  // Call kernel stub - fill iholes in each tile's particle buffer
  compactify_wrapper( chunk_manager->d_chunks, chunk_manager_ext->d_chunks,
                      chunk_manager->d_nchunks_max, chunk_manager->d_nchunks,
                      chunk_manager_ext->d_nchunks_max, chunk_manager_ext->d_nchunks,
                      chunk_size, d_ihole, d_n_hole,
                      d_n_comm_ext, d_num_par, d_n_recv_loc, nblocks, no_co->d_err );

  // check error code from kernel
  check_sort_error( no_co, "compactify" );

  // update num_par to buffer that fortran has a handle on
  // now, num_par_ contains the number of particles for each tile again, as normal
  chkErr(cudaMemcpy(num_par_, d_num_par, sizeof(int)*n_tils, cudaMemcpyDeviceToHost));

  // update species member data num_par
  for (int til=0; til<n_tils; ++til)
    num_par[til] = num_par_[til];

#ifdef __DEBUG__
  // first, need to recalculate block partition since number of particles has changed
  chkErr(cudaFree(d_block_start_tiles));
  chkErr(cudaFree(d_block_start_chunks));
  chkErr(cudaFree(d_block_num_par));

  get_block_partition( chunk_size );

  // now, validate particles
  validate( no_co, false, "after update boundary" );
#endif

  // Clean up memory
  chkErr(cudaFree(d_block_start_tiles));
  chkErr(cudaFree(d_block_start_chunks));
  chkErr(cudaFree(d_block_num_par));

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::copy_particles_ext_to_device( int *num_par_ ){

  int status = chunk_manager_ext->prune_and_pad_chunks( n_tils, num_par_ );
  if (!status) abort_program();

  for (int til=0; til<n_tils; ++til) {
    int status = chunk_manager_ext->copy_particles_to_device( til, num_par_[til] );
    if (!status) abort_program();
  }

  mem_transfer_manager->flush();

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::pad_chunks_ub2(int *n_recv_ext){

  // n_recv_ext is number of particles received from other nodes during host update bound

  // new number of particles we'll have after compactify
  int *num_par_new = (int *) malloc( sizeof(int)*n_tils );

  // number of holes that were left in particle array after sort
  int *n_hole = (int *) malloc( sizeof(int)*n_tils );
  chkErr(cudaMemcpy(n_hole, d_n_hole, sizeof(int)*n_tils, cudaMemcpyDeviceToHost));

  // number of particles recieved via local communicaitions in sort
  int *n_recv_loc = (int *) malloc( sizeof(int)*n_tils );
  chkErr(cudaMemcpy(n_recv_loc, d_n_recv_loc, sizeof(int)*n_tils, cudaMemcpyDeviceToHost));

  // calculate the new number of particles we'll have after the compactify
  // num_par is the old (pre-sort) number of particles
  for (int til=0; til<n_tils; ++til) {
    num_par_new[til]  = num_par[til] + n_recv_loc[til] + n_recv_ext[til] - n_hole[til];

    // only prune and pad if we gained particles (ie don't prune, there could be particles
    // in there that need to get compactified)
    if ( num_par_new[til] > num_par[til] ){
      int status = chunk_manager->prune_and_pad_chunks_til( til, num_par_new[til] );
      if (!status) abort_program();
    }
  }

  mem_transfer_manager->flush();

  chunk_manager->sync_state( n_tils, cudaMemcpyHostToDevice );

  free(num_par_new);
  free(n_hole);
  free(n_recv_loc);

}


// ---------------------------------------------------------------------------------------
// Check error code from sort/compactify kernels
// ---------------------------------------------------------------------------------------
void t_species_c::check_sort_error( t_node_conf_c *no_co, const char *msg ){

  // copy error code from host to device
  chkErr( cudaMemcpy( no_co->h_err, no_co->d_err, sizeof(int), cudaMemcpyDeviceToHost ) );

  if ( *(no_co->h_err) == p_err_nchunks_max ) {
    fprintf(stderr,"ERROR: tried to access an invalid element of chunks in %s. "
                   "Try increasing nchunks_pad\n", msg);
    abort_program();
  } else if ( *(no_co->h_err) == p_err_nchunks_max_ext ) {
    fprintf(stderr,"ERROR: tried to access an invalid element of chunks_ext in %s. "
                   "Try increasing nchunks_pad\n", msg);
    abort_program();
  } else if ( *(no_co->h_err) == p_err_pierre_ihole ) {
    fprintf(stderr,"ERROR: ran out of space in ihole in %s. "
                   "Try increasing ihole_size. Current ihole_size is %d.\n",
                   msg, ihole_size);
    abort_program();
  } else if ( *(no_co->h_err) != 0 ) {
    fprintf(stderr,"ERROR: mysterious error triggered in %s\n", msg);
    abort_program();
  }

}


// ---------------------------------------------------------------------------------------
// Check if all particle values are ok
// ---------------------------------------------------------------------------------------
void t_species_c::validate( t_node_conf_c *no_co, bool over, const char *msg ){

  validate_wrapper( chunk_manager->d_chunks, chunk_manager->d_nchunks, d_my_nx_p,
                    d_block_start_tiles, d_block_start_chunks,
                    d_block_num_par, chunk_pool->chunk_size, over,
                    move_num[0], move_num[1], move_num[2],
                    no_co->d_err, no_co->my_aid, nblocks );

  // check error code from kernel
  chkErr( cudaMemcpy( no_co->h_err, no_co->d_err, sizeof(int), cudaMemcpyDeviceToHost ) );

  if ( *(no_co->h_err) != 0 ) {
    fprintf(stderr,"ERROR: Validate species %s failed\n", msg);
    abort_program();
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::copy_report_from_device( bool copy_spec, bool report_energy ){

  // if necessary for diagnostics, copy particles from device to host
  if (copy_spec) {
    for (int til=0; til<n_tils; ++til) {
      int status = chunk_manager->copy_particles_from_device( til, num_par[til], num_par_max[til] );
      if (!status) abort_program();
    }
    mem_transfer_manager->flush();
  }

  // if necessary for diagnostics, copy energy from device to host
  if (report_energy) {
    chkErr(cudaMemcpy( energy, d_energy, sizeof(p_double)*n_tils, cudaMemcpyDeviceToHost ));
    for (int til=0; til<n_tils; ++til) {
      *energy_fortran[til] = energy[til];
    }
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::cleanup( int n_tils ) {

  my_nx_p.cleanup();
  chkErr(cudaFree(d_my_nx_p));
  chkErr(cudaFree(d_num_par));

  for( int i=0; i<n_tils; i++ ) {
    (d_ihole+i)->cleanup( true );
  }
  chkErr(cudaFree(d_ihole));

  chunk_manager->cleanup( n_tils );
  chunk_manager_ext->cleanup( n_tils );

  delete chunk_manager;
  delete chunk_manager_ext;

  delete[] num_par_max;
  delete[] num_par;
  delete[] energy_fortran;
  delete[] energy;
  chkErr(cudaFree(d_energy));

  for (int til=0; til<n_tils; ++til) {
    x[til].cleanup();
    p[til].cleanup();
    q[til].cleanup();
    ix[til].cleanup();
    x[til].~t_f1();
    p[til].~t_f1();
    q[til].~t_f0();
    ix[til].~t_f1();
  }

  free(x);
  free(p);
  free(q);
  free(ix);

  chkErr(cudaFree( d_n_hole ));
  chkErr(cudaFree( d_n_recv_loc ));
  chkErr(cudaFree( d_n_comm_ext ));
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::cleanup_dlb( int *tils_send_old, int n_send ) {

  int til;

  my_nx_p.cleanup();
  chkErr(cudaFree(d_my_nx_p));
  chkErr(cudaFree(d_num_par));

  // for now I just allocated a complete new on one the device upon dlb (see init_group_dlb)
  // so clean up everything
  // Could revise this choice if it is slow
  for( int i=0; i<n_tils; i++ ) {
    (d_ihole+i)->cleanup( true );
  }
  chkErr(cudaFree(d_ihole));

  chunk_manager->cleanup_dlb( tils_send_old, n_send );
  chunk_manager_ext->cleanup_dlb( tils_send_old, n_send );

  delete chunk_manager;
  delete chunk_manager_ext;

  delete[] num_par_max;
  delete[] num_par;
  delete[] energy_fortran;
  delete[] energy;
  chkErr(cudaFree(d_energy));

  for (int i=0; i<n_send; ++i) {
    til = tils_send_old[i];
    x[til].cleanup();
    p[til].cleanup();
    q[til].cleanup();
    ix[til].cleanup();
    x[til].~t_f1();
    p[til].~t_f1();
    q[til].~t_f0();
    ix[til].~t_f1();
  }

  free(x);
  free(p);
  free(q);
  free(ix);

  chkErr(cudaFree( d_n_hole ));
  chkErr(cudaFree( d_n_recv_loc ));
  chkErr(cudaFree( d_n_comm_ext ));
}


// ---------------------------------------------------------------------------------------
// return all chunks to the pool for this species/til
// ---------------------------------------------------------------------------------------
void t_species_c::delete_particles_from_device( int til ){

  int n_chunks_blow;

  // return chunks
  n_chunks_blow = chunk_manager->h_nchunks[til];
  chunk_manager->blow_chunks( til, n_chunks_blow );

  // return chunks_ext
  n_chunks_blow = chunk_manager_ext->h_nchunks[til];
  chunk_manager_ext->blow_chunks( til, n_chunks_blow );

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_species_c::test_particle_copy( void ) {

  PRINT( "Testing particle copy..." );

  int status;

  int til = 0;

  // duplicate buffers for test
  int num_par_cur = num_par[til];

  size_t til_mem_size_x      = sizeof(p_k_part) * p_x_dim * num_par_cur;
  size_t til_mem_size_p      = sizeof(p_k_part) * p_p_dim * num_par_cur;
  size_t til_mem_size_q      = sizeof(p_k_part) * num_par_cur;
  size_t til_mem_size_ix     = sizeof(int) * p_x_dim * num_par_cur;

  p_k_part *xtmp = (p_k_part *) malloc( til_mem_size_x );
  p_k_part *ptmp = (p_k_part *) malloc( til_mem_size_p );
  p_k_part *qtmp = (p_k_part *) malloc( til_mem_size_q );
  int *ixtmp = (int *) malloc( til_mem_size_ix );

  memcpy( xtmp, x[til].g, til_mem_size_x );
  memcpy( ptmp, p[til].g, til_mem_size_p );
  memcpy( qtmp, q[til].g, til_mem_size_q );
  memcpy( ixtmp, ix[til].g, til_mem_size_ix );

  PRINT( "copying particles to device..." );
  chunk_manager->copy_particles_to_device( til, num_par_cur );
  mem_transfer_manager->flush();

  // zero out arrays to confirm they're actually copied back
  memset( x[til].g, 0, til_mem_size_x );
  memset( p[til].g, 0, til_mem_size_p );
  memset( q[til].g, 0, til_mem_size_q );
  memset( ix[til].g, 0, til_mem_size_ix );

  PRINT( "copying particles back..." );
  status = chunk_manager->copy_particles_from_device( til, num_par_cur, num_par_max[til] );
  if (!status) abort_program();

  mem_transfer_manager->flush();

  status = 1;

  // int idx = 5120 * 4;
  // PRINT_FLOAT( x[til].g[ idx ] );
  // PRINT_FLOAT( xtmp[ idx ] );

  if( memcmp( xtmp, x[til].g, til_mem_size_x ) ) {
    PRINT( "ERROR: x copy failed.");
    status = 0;
  }

  if( memcmp( ptmp, p[til].g, til_mem_size_p ) ) {
    PRINT( "ERROR: p copy failed.");
    status = 0;
  }

  if( memcmp( qtmp, q[til].g, til_mem_size_q ) ) {
    PRINT( "ERROR: q copy failed.");
    status = 0;
  }

  if( memcmp( ixtmp, ix[til].g, til_mem_size_ix ) ) {
    PRINT( "ERROR: ix copy failed.");
    status = 0;
  }

  if( status ) {
    PRINT( "---> All tests succeeded." );
  }

  free( xtmp );
  free( ptmp );
  free( qtmp );
  free( ixtmp );

  exit(-1);

}
