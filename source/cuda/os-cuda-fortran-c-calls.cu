// ---------------------------------------------------------------------------------------
// This file contains all functions that are called from Fortran
// ---------------------------------------------------------------------------------------

#include <stdio.h>
#include <algorithm>
using namespace std;

#include "os-param-cuda.h"
#include "os-sys-cuda.h"
#include "os-sim-cuda.h"

extern "C" void init_node_conf( int n_tils, int my_aid, int no_num, int *neighbor_til_id,
                                int *shift ){

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_node_conf( n_tils, my_aid, no_num, neighbor_til_id, shift );

}

extern "C" void init_chunk_pool( int chunk_pool_size, int chunk_size, int n_tils ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_chunk_pool( chunk_pool_size, chunk_size, n_tils );

}

extern "C" void init_mem_transfer_manager( int com_buf_size ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_mem_transfer_manager( com_buf_size );

}

extern "C" void flush_mem_transfer_manager( void ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.flush_mem_transfer_manager();

}

extern "C" void init_flds_group( void ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_flds_group();

}

extern "C" void init_flds_tile( int til_id, int *jay_nx, int *jay_gc_num, p_k_fld *jay,
                                int *emf_nx, int *emf_gc_num, p_k_fld *e, p_k_fld *b,
                                bool if_last_til ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_flds_tile( til_id, jay_nx, jay_gc_num, jay, emf_nx, emf_gc_num, e, b, if_last_til );

}

extern "C" void init_spec_group( int spec_id, p_k_part rqm, bool free_stream, p_double *dx,
                                 int interpolation, bool grid_center, p_double dt,
                                 int ihole_size, p_double ihole_size_frac, int nchunks_pad,
                                 p_double nchunks_pad_frac, p_double chunk_growth_factor,
                                 int nblocks_requested, int tiles_per_node ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_spec_group( spec_id, rqm, free_stream, dx, interpolation, grid_center, dt,
                       ihole_size, ihole_size_frac, nchunks_pad, nchunks_pad_frac,
                       chunk_growth_factor, nblocks_requested, tiles_per_node );

}

extern "C" void init_spec_tile( int til_id, int spec_id, int num_par_max, int num_par,
                                p_double *energy, int *my_nx_p, p_k_part *x, p_k_part *p,
                                p_k_part *q, int *ix, int nchunks_max_init, bool if_last_til ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_spec_tile( til_id, spec_id, num_par_max, num_par, energy, my_nx_p, x, p, q, 
                      ix, nchunks_max_init, if_last_til );

}

extern "C" void init_tile_last_dlb( void ){

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_tile_last_dlb();

}

extern "C" void sync_buffer_spec( int til_id, int spec_id, int num_par_max,
                                  p_k_part *x, p_k_part *p, p_k_part *q, int *ix ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.sync_buffer_spec( til_id, spec_id, num_par_max, x, p, q, ix );

}

extern "C" void copy_particles_to_device( void ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.copy_particles_to_device();

}

extern "C" void copy_particles_ext_from_device( int *num_par ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.copy_particles_ext_from_device(num_par);

}

extern "C" void pre_advance_deposit_cuda( void ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.pre_advance_deposit();

}

extern "C" void advance_deposit_cuda( bool *report_energy ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.advance_deposit( report_energy );

}

extern "C" void post_advance_deposit_cuda( void ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.post_advance_deposit();

}

extern "C" void move_window_part_cuda( int *nmove, int *num_par ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.move_window_part( nmove, num_par );

}

extern "C" void update_boundary_part_cuda_1( int *num_par ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.update_boundary_part_1( num_par );

}

extern "C" void update_boundary_part_cuda_2( int *num_par ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.update_boundary_part_2( num_par );

}

extern "C" void init_dlb_cuda_( int n_tils_new, int *tils_stnry, int *tils_stnry_old,
				int n_stnry, int *tils_send_old, int n_send,
				int ihole_size, p_double ihole_size_frac, 
				int nchunks_pad, p_double nchunks_pad_frac ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_dlb_cuda( n_tils_new, tils_stnry, tils_stnry_old, n_stnry, tils_send_old,
		     n_send, ihole_size, ihole_size_frac, nchunks_pad, nchunks_pad_frac );
 
}

extern "C" void copy_particles_to_device_dlb( int n_tils_recv, int *tils_recv ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.copy_particles_to_device_dlb( n_tils_recv, tils_recv);
 
}

extern "C" void copy_particles_from_device_dlb( int n_tils_send, int *tils_send ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.copy_particles_from_device_dlb( n_tils_send, tils_send );
 
}

extern "C" void init_node_conf_dlb( int n_tils, int my_aid, int no_num, int *neighbor_til_id,
                                    int *shift ){

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.init_node_conf_dlb( n_tils, my_aid, no_num, neighbor_til_id, shift );

}

extern "C" void copy_report_from_device( bool *copy_spec, bool *report_energy ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.copy_report_from_device( copy_spec, report_energy );

}

extern "C" void cleanup_cuda( void ) {

  t_simulation_c &sim = can_i_haz_simulation_c();

  sim.cleanup();

}
