#ifndef OS_SIM_CUDA_H
#define OS_SIM_CUDA_H

#include "os-nconf-cuda.h"
#include "os-current-cuda.h"
#include "os-emf-cuda.h"
#include "os-species-cuda.h"
#include "os-chunk-pool-cuda.h"
#include "os-mem-transfer-manager-cuda.h"

// to store global cuda state
class t_simulation_c {
public:
  t_node_conf_c *no_co;
  t_current_c *jay;
  t_emf_c *emf;
  t_species_c *species;
  t_chunk_pool *chunk_pool;
  t_mem_transfer_manager *mem_transfer_manager;

  t_simulation_c( void );

  // Member functions serving purpose analagous to init_sim_group() in fortran.
  // Need a bunch of different initializatoin routines(can't just have one main init like
  // in fortran) due to limitations on what can be passed from fortran to C
  virtual void init_node_conf( int n_tils, int my_aid, int no_num, int *neighbor_til_id,
                               int *shift );
  virtual void init_chunk_pool( int chunk_pool_size, int chunk_size, int n_tils );
  virtual void init_mem_transfer_manager( int com_buf_size );
  virtual void flush_mem_transfer_manager( void );
  virtual void init_flds_group( void );
  virtual void init_flds_tile( int til_id, int *jay_nx, int *jay_gc_num, p_k_fld *jay_,
                               int *emf_nx, int *emf_gc_num, p_k_fld *e, p_k_fld *b,
                               bool if_last_til );
  virtual void init_spec_group( int spec_id, p_k_part rqm, bool free_stream, p_double *dx,
                                int interpolation, bool grid_center, p_double dt,
                                int ihole_size, p_double ihole_size_frac, int nchunks_pad,
                                p_double nchunks_pad_frac, p_double chunk_growth_factor,
                                int nblocks_requested, int tiles_per_node );
  virtual void init_spec_tile( int til_id, int spec_id, int num_par_max, int num_par,
                               p_double *energy, int *my_nx_p, p_k_part *x, p_k_part *p,
                               p_k_part *q, int *ix, int nchunks_max_init, bool if_last_til );
  virtual void init_tile_last_dlb();
  virtual void sync_buffer_spec( int til_id, int spec_id, int num_par_max,
                                 p_k_part *x, p_k_part *p, p_k_part *q, int *ix );
  virtual void copy_particles_to_device( void );
  virtual void pre_advance_deposit();
  virtual void advance_deposit( bool *report_energy );
  virtual void post_advance_deposit();
  virtual void move_window_part( int *nmove, int *num_par );
  virtual void update_boundary_part_1( int *num_par );
  virtual void copy_particles_ext_from_device( int *num_par );
  virtual void update_boundary_part_2( int *num_par );
  virtual void copy_report_from_device( bool *copy_spec, bool *report_energy );

  // Load balancing procedures
  virtual void copy_particles_from_device_dlb(int n_send, int *tils );
  virtual void init_dlb_cuda( int n_tils_new, int *h_tils_stnry, int *h_tils_stnry_old,
                              int n_stnry, int *tils_send_old, int n_send,
			      int ihole_size, p_double ihole_size_frac,
			      int nchunks_pad, p_double nchunks_pad_frac );
  virtual void copy_particles_to_device_dlb( int n_tils_recv, int *tils_recv );
  virtual void init_node_conf_dlb( int n_tils, int my_aid, int no_num, int *neighbor_til_id,
                                   int *shift );

  virtual void cleanup( void );
};

// Use this function to refer to SIM_CUDA
t_simulation_c& can_i_haz_simulation_c(void);

#endif
