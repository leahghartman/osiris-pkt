#ifndef OS_NCONF_CUDA_H
#define OS_NCONF_CUDA_H

#include <cuda_runtime.h>
#include "os-vdf-cuda.h"

class t_node_conf_c {
public:
  int n_tils;
  int my_aid;

  size_t maxSharedMemoryPerBlock;
  size_t maxSharedMemoryPerBlockOptin;

  // dimensions are ( 3^p_x_dim, n_tils )
  t_f1<int> *d_neighbor_til_id;
  t_f2<int> *d_shift;

  // for kernel error handling
  int *h_err, *d_err;

  cudaStream_t s_j, s_adv_dep;

  t_node_conf_c( void );
  virtual void init( int n_tils_, int my_aid, int no_num, int *neighbor_til_id,
                     int *shift );
  virtual void init_dlb( int n_tils_, int my_aid_, int no_num, int *neighbor_til_id,
                              int *shift );
  virtual void cleanup( void );
  virtual void cleanup_dlb( void );
  virtual ~t_node_conf_c( void ) {};
};

#endif
