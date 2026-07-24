#include <stdio.h>
#include <cstdlib>

#include "os-sim-cuda.h"
#include "os-spec-push-cuda.h"
#include "os-dlb-cuda.h"
#include "os-sys-cuda.h"
#include "os-param-cuda.h"

// define SIM_CUDA
t_simulation_c SIM_CUDA;

// Use this function to refer to SIM_CUDA
t_simulation_c &can_i_haz_simulation_c(void) {
  return SIM_CUDA;
}

// Object definition
// Assume all til_id and spec_id come in as 1-indexed, go out as 0-indexed
t_simulation_c::t_simulation_c( void ) {

  no_co = NULL;
  jay = NULL;
  emf = NULL;
  species = NULL;
  chunk_pool = NULL;
  mem_transfer_manager = NULL;

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_node_conf( int n_tils, int my_aid, int no_num,
                                     int *neighbor_til_id, int *shift ) {

#ifdef __DEBUG__
  if ( my_aid==1 )
    printf("\n****** WARNING: Running with __DEBUG__ defined - performance will be limited ******\n\n");
#endif

  // Allocate objects
  if (!no_co) no_co = new t_node_conf_c;

  // Initialize object
  no_co->init( n_tils, my_aid, no_num, neighbor_til_id, shift );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_chunk_pool( int chunk_pool_size, int chunk_size, int n_tils ) {
  if (!chunk_pool) chunk_pool = new t_chunk_pool( chunk_pool_size, chunk_size, n_tils,
                                                  no_co->my_aid );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_mem_transfer_manager( int com_buf_size ) {
  if (!mem_transfer_manager) mem_transfer_manager = new t_mem_transfer_manager( com_buf_size );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::flush_mem_transfer_manager( void ) {
  mem_transfer_manager->flush();
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_flds_group( void ) {
  // Allocate objects
  if (!jay) jay = new t_current_c;
  if (!emf) emf = new t_emf_c;

  // Initialize arrays
  jay->init_group( no_co->n_tils );
  emf->init_group( no_co->n_tils );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_flds_tile( int til_id, int *jay_nx, int *jay_gc_num,
                                     p_k_fld *jay_, int *emf_nx, int *emf_gc_num,
                                     p_k_fld *e, p_k_fld *b, bool if_last_til ) {
  // Initialize data for each tile
  jay->init_tile( til_id-1, jay_nx, jay_gc_num, jay_, no_co->n_tils,
                  mem_transfer_manager, if_last_til );
  emf->init_tile( til_id-1, emf_nx, emf_gc_num, e, b, mem_transfer_manager );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_spec_group( int spec_id, p_k_part rqm, bool free_stream, p_double *dx,
                                      int interpolation, bool grid_center, p_double dt,
                                      int ihole_size, p_double ihole_size_frac, int nchunks_pad,
                                      p_double nchunks_pad_frac, p_double chunk_growth_factor,
                                      int nblocks_requested, int tiles_per_node ) {

  t_species_c *spec = species;
  t_species_c *tail = NULL;
  for (int i=0; i<spec_id-1; i++) {
    tail = spec;
    spec = spec->next;
  }
  if (spec) {
    fprintf(stderr, "ERROR: Species pointer is somehow already associated.");
    abort_program();
  }

  spec = new t_species_c;
  spec->init_group( no_co->n_tils, chunk_pool, mem_transfer_manager,
                    rqm, free_stream, dx, interpolation,
                    grid_center, dt, ihole_size, ihole_size_frac, nchunks_pad,
                    nchunks_pad_frac, chunk_growth_factor, nblocks_requested,
                    tiles_per_node );

  if (spec_id > 1) {
    tail->next = spec;
  } else {
    // is this line necessary?
    species = spec;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_spec_tile( int til_id, int spec_id, int num_par_max,
                                     int num_par, p_double *energy, int *my_nx_p,
                                     p_k_part *x, p_k_part *p, p_k_part *q, int *ix,
                                     int nchunks_max_init, bool if_last_til ) {

  t_species_c *spec = species;
  for (int i=0; i<spec_id-1; i++)
    spec = spec->next;

  spec->init_tile( til_id-1, num_par_max, num_par, energy, my_nx_p, x, p, q, ix,
                   nchunks_max_init, if_last_til );

}


// ---------------------------------------------------------------------------------------
// Final initialization step after dynamic load balancing
// Only gets called if this MPI rank didn't receive any tiles
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_tile_last_dlb(){

  jay->init_tile_last( no_co->n_tils, mem_transfer_manager );

  t_species_c *spec = species;
  while( spec ) {
    spec->init_tile_last();

    spec->chunk_manager->sync_state( no_co->n_tils, cudaMemcpyHostToDevice );

    spec = spec->next;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::sync_buffer_spec( int til_id, int spec_id, int num_par_max,
                                       p_k_part *x, p_k_part *p, p_k_part *q, int *ix ) {

  t_species_c *spec = species;
  for (int i=0; i<spec_id-1; i++)
    spec = spec->next;

  spec->sync_buffer( til_id-1, num_par_max, x, p, q, ix );

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::copy_particles_to_device( void ) {

  int status;

  // debug: make sure copying to and from device works correctly.
#if 0
  species->test_particle_copy();
#endif

  // copy particles in, one tile at a time
  t_species_c *spec = species;
  while( spec ){
    // Get chunks from the pool to make room for the particles
    status = spec->chunk_manager->prune_and_pad_chunks( no_co->n_tils, spec->num_par );
    if (!status) abort_program();

    for( int til=0; til<no_co->n_tils; til++ ) {
      status = spec->chunk_manager->copy_particles_to_device( til, spec->num_par[til] );
      if (!status) abort_program();
    }

    mem_transfer_manager->flush();

    spec = spec->next;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::pre_advance_deposit() {

  // copy host memory to device, one tile at a time
  for (int til=0; til < no_co->n_tils; ++til) {
    emf->buffered_copy_to_device( til, mem_transfer_manager );
  }
  mem_transfer_manager->flush();

  // zero out current array
  zero_jay_wrapper( jay->d_prefix_size, jay->prefix_size[no_co->n_tils+1], jay,
                    jay->prefix_size[no_co->n_tils], no_co->s_j );

  t_species_c *spec = species;
  int j = 0;
  while(spec) {
    spec->pre_advance_deposit( no_co, jay, emf );
    spec = spec->next;
    j += 1;
  }

}


// ---------------------------------------------------------------------------------------
// we separate this from pre_advance_deposit and post_advance_deposit so that we can 
// isolate the timings for just the advance_deposit kernels in fortran
// ---------------------------------------------------------------------------------------
void t_simulation_c::advance_deposit( bool *report_energy ) {

  t_species_c *spec = species;
  int j = 0;
  while(spec) {
    spec->advance_deposit( no_co, jay, emf, report_energy[j] );
    spec = spec->next;
    j += 1;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::post_advance_deposit() {

  // copy result from device to host, one tile at a time
  chkErr(cudaStreamSynchronize( no_co->s_adv_dep ));
  for (int til=0; til < no_co->n_tils; ++til) {
    jay->buffered_copy_from_device( til, mem_transfer_manager );
  }
  mem_transfer_manager->flush();

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::move_window_part( int *nmove, int *num_par ) {

  int i = 0;
  int j = 0;
  int num_par_new[no_co->n_tils];

  t_species_c *spec = species;

  while(spec) {
    // Get chunks from the pool to make room for the particles
    for (int til=0; til<no_co->n_tils; ++til) {
      num_par_new[til] = spec->num_par[til] + (num_par+i)[til];
    }
    int status = spec->chunk_manager->prune_and_pad_chunks( no_co->n_tils, num_par_new );
    if (!status) abort_program();

    spec->move_window( nmove, num_par+i, no_co );

    spec = spec->next;
    i += no_co->n_tils;
    j += 1;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::update_boundary_part_1( int *num_par ) {
  t_species_c *spec = species;
  int i = 0;
  while(spec) {
    spec->update_boundary_1( no_co, num_par+i );
    spec = spec->next;
    i += no_co->n_tils;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::copy_particles_ext_from_device( int *num_par ) {
  t_species_c *spec = species;
  int i = 0;
  while( spec ){
    spec->copy_particles_ext_from_device( num_par+i );
    spec = spec->next;
    i += no_co->n_tils;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::update_boundary_part_2( int *num_par ) {
  t_species_c *spec = species;
  int i = 0;
  while(spec) {
    spec->update_boundary_2( no_co, num_par+i );
    spec = spec->next;
    i += no_co->n_tils;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::copy_report_from_device( bool *copy_spec, bool *report_energy ) {
  t_species_c *spec = species;
  int j = 0;
  while(spec) {
    spec->copy_report_from_device( copy_spec[j], report_energy[j] );
    spec = spec->next;
    j += 1;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_node_conf_dlb( int n_tils, int my_aid, int no_num,
                                         int *neighbor_til_id, int *shift ) {
  // Initialize object
  no_co->init_dlb( n_tils, my_aid, no_num, neighbor_til_id, shift );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::copy_particles_from_device_dlb( int n_tils_send, int *tils_send ) {

  int til;
  t_species_c *spec;

  for (int i=0; i<n_tils_send; i++) {
    til = tils_send[i];
    spec = species;
    while(spec) {
      int status = spec->chunk_manager->copy_particles_from_device(
                                                                 til, spec->num_par[til], 
                                                                 spec->num_par_max[til] );
      if (!status) abort_program();

      spec = spec->next;
    }
  }

  // flush any remaining data from the buffer
  mem_transfer_manager->flush();

  // return all chunks belonging to these tiles to the chunk pool
  for (int i=0; i<n_tils_send; i++) {
    til = tils_send[i];
    spec = species;
    while(spec) {
      spec->delete_particles_from_device( til );
      spec = spec->next;
    }
  }

  // flush device to device mem transfers from returning chunks to the pool
  mem_transfer_manager->flush();

}


// ---------------------------------------------------------------------------------------
// reallocate simulation objects to accomadate new, post dynamic load balance, number 
// of tiles, and copy pointers to stationary tiles to the new structures
// ---------------------------------------------------------------------------------------
void t_simulation_c::init_dlb_cuda( int n_tils_new, int *h_tils_stnry, int *h_tils_stnry_old,
				    int n_stnry, int *tils_send_old, int n_send,
				    int ihole_size, p_double ihole_size_frac,
				    int nchunks_pad, p_double nchunks_pad_frac  ){

  // post and pre load balance tile id's, respectively
  // these are only different when using Viktor topology
  int til, til_old;

  // reinitialize flds group
  t_current_c *jay_temp;
  jay_temp = new t_current_c;
  jay_temp->init_group( n_tils_new );

  t_emf_c *emf_temp;
  emf_temp = new t_emf_c;
  emf_temp->init_group( n_tils_new );

  // reinitialize species group
  t_species_c *species_temp;

  t_species_c *spec_temp = NULL;
  t_species_c *tail_temp = NULL;

  t_species_c *spec = species;
  int spec_id = 0;
  while( spec ) {
    spec_temp = new t_species_c;
    spec_temp->init_group_dlb( n_tils_new, spec, ihole_size, ihole_size_frac,
			       nchunks_pad, nchunks_pad_frac );

    if (spec_id > 0) {
      tail_temp->next = spec_temp;
    } else {
      species_temp = spec_temp;
    }

    tail_temp = spec_temp;
    spec_temp = spec_temp->next;

    spec = spec->next;
    spec_id += 1;
  }

  // copy host pointers to stationary tiles into temp arrays
  for (int i=0; i<n_stnry; i++) {
    til = h_tils_stnry[i];
    til_old = h_tils_stnry_old[i];

    jay_temp->init_tile_dlb( til, til_old, jay );

    emf_temp->init_tile_dlb( til, til_old, emf );

    t_species_c *spec      = species;
    t_species_c *spec_temp = species_temp;
    while( spec ) {
      spec_temp->init_tile_dlb( til, til_old, spec );

      spec_temp = species_temp->next;
      spec      = spec->next;
    }
  }
  // flush after calling prune_and_pad_chunks (sync state happens later)
  mem_transfer_manager->flush();

  // Now we deal with the device arrays:

  // copy list of stationary tiles to device
  if ( n_stnry>0 ){
    int *d_tils_stnry;
    int *d_tils_stnry_old;

    chkErr(cudaMalloc(reinterpret_cast<void **>(&d_tils_stnry), sizeof(int)*n_stnry ));
    chkErr(cudaMalloc(reinterpret_cast<void **>(&d_tils_stnry_old), sizeof(int)*n_stnry ));

    chkErr(cudaMemcpy(d_tils_stnry, h_tils_stnry, sizeof(int)*n_stnry,
                      cudaMemcpyHostToDevice));
    chkErr(cudaMemcpy(d_tils_stnry_old, h_tils_stnry_old, sizeof(int)*n_stnry,
                      cudaMemcpyHostToDevice));

    // block partition: one thread per tile
    int nblocks = CEIL_DIV( n_stnry, __BLOCK_SIZE__ );

    // call kernels to copy device pointers to stationary tiles into temp arrays
    // Note: we only need to handle fields here, device arrays for species are just handled
    // with a few memcopies
    copy_tils_stnry_wrapper( d_tils_stnry, d_tils_stnry_old, n_stnry, jay, jay_temp,
                             emf, emf_temp, nblocks );

    cudaFree( d_tils_stnry );
    cudaFree( d_tils_stnry_old );
  }

  // cleanup no longer needed objects, temp objects become the new objects
  jay->cleanup_dlb( tils_send_old, n_send );
  delete(jay);
  jay = jay_temp;

  emf->cleanup_dlb( tils_send_old, n_send );
  delete(emf);
  emf = emf_temp;

  spec = species;
  t_species_c *tail = NULL;
  while(spec) {
    spec->cleanup_dlb( tils_send_old, n_send ); 
    tail = spec;
    spec = spec->next;
    delete tail;
  }
  species = species_temp;

  no_co->cleanup_dlb();

  // some other housekeeping items as long as we're at it
  chunk_pool->n_tils = n_tils_new;

}


// ---------------------------------------------------------------------------------------
// Here we copy particles to device from tiles we received during load balancing
// ---------------------------------------------------------------------------------------
void t_simulation_c::copy_particles_to_device_dlb( int n_tils_recv, int *tils_recv ) {

  int til;
  int status;

  // debug: make sure copying to and from device works correctly.
#if 0
  species->test_particle_copy();
#endif

  // copy particles in, one tile at a time
  t_species_c *spec = species;
  while( spec ){
    // Get chunks from the pool to make room for the particles
    status = spec->chunk_manager->prune_and_pad_chunks( no_co->n_tils, spec->num_par );
    if (!status) abort_program();

    for( int i=0; i<n_tils_recv; i++ ) {
      til = tils_recv[i];
      status = spec->chunk_manager->copy_particles_to_device( til, spec->num_par[til] );
      if (!status) abort_program();
    }

    spec = spec->next;
  }

  mem_transfer_manager->flush();

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_simulation_c::cleanup( void ) {

  jay->cleanup( no_co->n_tils );
  emf->cleanup( no_co->n_tils );
  delete jay;
  delete emf;

  t_species_c *spec = species;
  t_species_c *tail = NULL;
  while(spec) {
    spec->cleanup( no_co->n_tils );
    tail = spec;
    spec = spec->next;
    delete tail;
  }

  mem_transfer_manager->cleanup();
  delete mem_transfer_manager;

  chunk_pool->cleanup();
  delete chunk_pool;

  no_co->cleanup();
  delete no_co;
}
