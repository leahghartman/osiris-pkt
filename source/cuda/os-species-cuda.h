#ifndef OS_SPEC_CUDA_H
#define OS_SPEC_CUDA_H

#include "os-param-cuda.h"
#include "os-sys-cuda.h"
#include "os-vdf-cuda.h"
#include "os-nconf-cuda.h"
#include "os-current-cuda.h"
#include "os-emf-cuda.h"
#include "os-chunk-pool-cuda.h"
#include "os-mem-transfer-manager-cuda.h"
#include "os-species-chunk-manager-cuda.h"

class t_species_c {
public:

  int n_tils;

  // Parameters identical for all tiles
  p_k_part rqm;
  bool free_stream;
  p_double dx[p_max_dim];
  int interpolation;
  bool grid_center;
  p_double dt;

  int move_num[p_max_dim];

  // number of cells per dimension per tile (third element of my_nx_p from Fortran)
  // dimension(p_x_dim,n_tils)
  t_f1<int> my_nx_p;
  // int array of above to avoid complication
  int *d_my_nx_p;

  t_f0<int> *d_ihole; // n_tils array of t_f0 containing arrays for the holes in all species
  int ihole_size;     // initial size of each buffer in d_ihole

  // number of particles associated with chunks
  int *num_par_max, *num_par, *d_num_par;
  p_double **energy_fortran, *energy, *d_energy;

  // One for each tile, these point to CPU particle buffers
  t_f1<p_k_part> *x, *p;
  t_f0<p_k_part> *q;
  t_f1<int> *ix;

  int *d_n_hole;     // num particles leaving each tile (num holes in particle array)
  int *d_n_recv_loc; // num particles each tile is receiving from tiles on same node
  int *d_n_comm_ext; // num particles each tile is receiving from/sending to tiles on
                     // other nodes

  t_species_c *next;

  // chunk handling stuff
  t_chunk_pool *chunk_pool; // pointer to global chunk_pool
  t_mem_transfer_manager *mem_transfer_manager; // pointer to global mem_transfer_manager 

  t_species_chunk_manager *chunk_manager;
  t_species_chunk_manager *chunk_manager_ext;

  int nchunks_pad;                        // number of empty chunks to have before sort

  // kernel invocation parameters
  int nblocks_requested;

  // block partition buffers
  int nblocks;

  int *d_block_start_tiles;
  int *d_block_start_chunks;
  int *d_block_num_par;

  t_species_c( void );

  virtual void init_group( int n_tils_, t_chunk_pool *chunk_pool_,
                           t_mem_transfer_manager *mem_transfer_manager_,
                           p_k_part rqm_, bool free_stream_, p_double *dx_,
                           int interpolation_, bool grid_center_, p_double dt_,
                           int ihole_size_, p_double ihole_size_frac, int nchunks_pad_,
                           p_double nchunks_pad_frac, p_double chunk_growth_factor_,
                           int nblocks_requested_, int tiles_per_node );

  
  virtual void init_tile( int til, int num_par_max_, int num_par_, p_double *energy_,
                          int *my_nx_p_, p_k_part *x_, p_k_part *p_, p_k_part *q_,
                          int *ix_, int nchunks_max_init_, bool if_last_til );

  virtual void init_tile_last();

  virtual void init_group_dlb( int n_tils_new, t_species_c *spec,
			       int ihole_size_, p_double ihole_size_frac,
			       int nchunks_pad, p_double nchunks_pad_frac );
  
  virtual void init_tile_dlb( int til, int til_old, t_species_c *spec );

  virtual void sync_buffer( int til, int num_par_max_, p_k_part *x_, p_k_part *p_,
                            p_k_part *q_, int *ix_ );

  virtual void get_block_partition( int chunk_size );

  virtual void pre_advance_deposit( t_node_conf_c *no_co, t_current_c *jay, t_emf_c *emf );

  virtual void advance_deposit( t_node_conf_c *no_co, t_current_c *jay, t_emf_c *emf,
                                bool report_energy );

  virtual void move_window( int *nmove, int *num_par_inject, t_node_conf_c *no_co );

  virtual void update_boundary_1( t_node_conf_c *no_co, int *num_par_ );

  virtual void copy_particles_ext_from_device( int *num_par_ );
  
  virtual void update_boundary_2( t_node_conf_c *no_co, int *num_par_ );

  virtual void copy_particles_ext_to_device( int *num_par_ );

  virtual void pad_chunks_ub2(int *n_recv_ext);

  virtual void check_sort_error( t_node_conf_c *no_co, const char *msg );

  virtual void validate( t_node_conf_c *no_co, bool over, const char *msg );

  virtual void copy_report_from_device( bool copy_spec, bool report_energy );
  
  virtual void cleanup( int n_tils );
  virtual void cleanup_dlb( int *tils_send_old, int n_send );

  virtual void delete_particles_from_device( int til );
  
  virtual ~t_species_c( void ) {};

  void test_particle_copy( void );

};



#endif
