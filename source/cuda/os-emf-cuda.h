#ifndef OS_EMF_CUDA_H
#define OS_EMF_CUDA_H

#include <cuda_runtime.h>
#include "os-param-cuda.h"
#include "os-vdf-cuda.h"
#include "os-mem-transfer-manager-cuda.h"

class t_emf_c {
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
  t_f1<p_k_fld> *h_ed1, *h_bd1, *d_ed1, *d_bd1;
  t_f2<p_k_fld> *h_ed2, *h_bd2, *d_ed2, *d_bd2;
  t_f3<p_k_fld> *h_ed3, *h_bd3, *d_ed3, *d_bd3;

  // host array of pointers to device field data arrays (ie they point directly
  // to device field data arrays, rather than to an object from which you can
  // get the field data pointer)
  // This enables us to copy fields directly from device to host, without having
  // to do an extra memcopy so we can get this info out of d_jd1/d_jd2/d_jd3
  p_k_fld **h_e_data_arr, **h_b_data_arr;

  t_emf_c( void );
  virtual void init_group( int n_tils );
  virtual void init_tile( int til, int *nx_, int *gc_num_, p_k_fld *e, p_k_fld *b,
                          t_mem_transfer_manager *mem_transfer_manager );
  virtual void init_tile_dlb( int til, int til_old, t_emf_c *emf );
  virtual void buffered_copy_to_device( int til, t_mem_transfer_manager *mem_transfer_manager );
  virtual void buffered_copy_from_device( int til, t_mem_transfer_manager *mem_transfer_manager );
  virtual void cleanup( int n_tils );
  virtual void cleanup_dlb( int *tils_send_old, int n_send );
  virtual ~t_emf_c( void ) {};

};

#endif
