#ifndef OS_CURRENT_CUDA_H
#define OS_CURRENT_CUDA_H

#include <cuda_runtime.h>
#include "os-param-cuda.h"
#include "os-vdf-cuda.h"
#include "os-mem-transfer-manager-cuda.h"

class t_current_c {
public:
  // These arrays are *,n_tils or *,*,n_tils in size, containing information about
  // all tiles
  t_f1<int> *nx;
  t_f2<int> *gc_num;

  t_f1<int> *dims;

  // number of cells in largest tile on this node (including guard cells)
  // used to allocate shmem in dudt/adv dep
  int nx_tot_max; 

  // These will be arrays n_tils in size containing field data for each tile
  t_f1<p_k_fld> *h_jd1, *d_jd1;
  t_f2<p_k_fld> *h_jd2, *d_jd2;
  t_f3<p_k_fld> *h_jd3, *d_jd3;

  // host array of pointers to device field data arrays (ie they point directly
  // to device field data arrays, rather than to an object from which you can
  // get the field data pointer)
  // This enables us to copy fields directly from device to host, without having
  // to do an extra memcopy so we can get this info out of d_jd1/d_jd2/d_jd3
  p_k_fld **h_j_data_arr;

  // used to determine kernel configuration parameters when zeroing current array
  int *prefix_size, *d_prefix_size;

  t_current_c( void );
  virtual void init_group( int n_tils );
  virtual void init_tile( int til, int *nx_, int *gc_num_, p_k_fld *jay,
                          int n_tils, t_mem_transfer_manager *mem_transfer_manager,
                          bool if_last_til );
  virtual void init_tile_last( int n_tils, t_mem_transfer_manager *mem_transfer_manager );
  virtual void init_tile_dlb( int til, int til_old, t_current_c *jay );
  virtual void buffered_copy_from_device( int til,
                                          t_mem_transfer_manager *mem_transfer_manager );
  virtual void cleanup( int n_tils );
  virtual void cleanup_dlb( int *tils_send_old, int n_send );
  virtual ~t_current_c( void ) {};
};

#endif
