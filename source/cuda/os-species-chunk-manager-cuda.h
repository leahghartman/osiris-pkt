#ifndef OS_SPECIES_CHUNK_MANAGER_CUDA_H
#define OS_SPECIES_CHUNK_MANAGER_CUDA_H

#include "os-param-cuda.h"
#include "os-sys-cuda.h"
#include "os-chunk-pool-cuda.h"
#include "os-mem-transfer-manager-cuda.h"




class t_species_chunk_manager{
public:

  // associated particle buffers for transfer to and from CPU
  t_f1<p_k_part> *x, *p;
  t_f0<p_k_part> *q;
  t_f1<int> *ix;

  // these sets of variables are typically the same but may be updated from either the host
  // or device. synchronize using the sync_state function.

  // array of array of device chunk pointers for each tile
  // d_chunks is a triple device pointer
  // h_chunks is a host pointer to double device pointers
  t_chunk ***h_chunks;
  t_chunk ***d_chunks;

  // max number of chunks each tile can store
  int *h_nchunks_max;
  int *d_nchunks_max;

  // number of chunks currently held by each tile
  int *h_nchunks;
  int *d_nchunks;

  // track whether sync is needed between host and device
  bool h_state_changed;
  // int *h_chunks_state_changed;

  // number of empty chunks to maintain on end of h_chunks and d_chunks
  int nchunks_pad;

  p_double chunk_growth_factor;   // specifies how much to grow **d_chunks buffer if we run out of space

  t_chunk_pool *chunk_pool; // pointer to global chunk_pool
  t_mem_transfer_manager *mem_transfer_manager; // pointer to global mem_transfer_manager

  t_species_chunk_manager( void );

  void init_group( t_f1<p_k_part> *x, t_f1<p_k_part> *p, t_f0<p_k_part> *q, t_f1<int> *ix,
                   int n_tils, int nchunks_pad, int chunk_growth_factor,
                   t_chunk_pool *chunk_pool_, t_mem_transfer_manager *mem_transfer_manager_ );

  void init_tile( int til, int nchunks_max, int n_tils, bool if_last_til );

  virtual void init_tile_last( int n_tils );

  void init_tile_dlb( int til, t_chunk **d_chunk, int nchunks_init, int nchunks_max_init );

  void cleanup( int n_tils );
  void cleanup_dlb( int *tils_send_old, int n_send );

  void sync_state( int n_tils, cudaMemcpyKind kind );

  // chunk management
  int mo_chunks( int til, int nchunks_recv );
  int blow_chunks( int til, int nchunks_blow );
  int prune_and_pad_chunks( int n_tils, int *num_par );
  int prune_and_pad_chunks_dry( int n_tils, int *num_par );
  int prune_and_pad_chunks_til( int til, int num_par );

  // data transfer functions
  int copy_particles_to_device( int til, int num_par_til );
  int append_particles_to_device( int til, int num_par_til, int num_par_present );
  int copy_particles_from_device( int til, int num_par_til, int num_par_max_til );

};






#endif // OS_SPECIES_CHUNK_MANAGER_H
