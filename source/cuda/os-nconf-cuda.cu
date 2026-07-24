#include <stdio.h>
#include <cuda_runtime.h>
#include <algorithm>
using namespace std;

#include "os-nconf-cuda.h"
#include "os-param-cuda.h"
#include "os-sys-cuda.h"

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
t_node_conf_c::t_node_conf_c( void ) {

  n_tils = 0;
  my_aid = 0;
  d_neighbor_til_id = NULL;
  d_shift = NULL;

}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_node_conf_c::init( int n_tils_, int my_aid_, int no_num, int *neighbor_til_id,
                          int *shift ) {

  n_tils     = n_tils_;
  my_aid     = my_aid_;

  // set device for this mpi rank
  // Currently we are assuming that mpi ranks whose my_aid's are close are also on the
  // same node
  int device_count = 0;
  cudaGetDeviceCount( &device_count );

  if (my_aid==1) {
    printf("   Found %d GPUs per physical node\n", device_count);
    printf("   Running with block size %d\n", __BLOCK_SIZE__);
  }

  int device = my_aid % device_count;

  chkErrStDv(cudaSetDevice(device), device_count);

  // get the default max shared memory per block, and the max shared mem per block
  // which can be opted into using cudaFuncSetAttribute
  int value;

  cudaDeviceGetAttribute( &value, cudaDevAttrMaxSharedMemoryPerBlock, device );
  maxSharedMemoryPerBlock = size_t(value);

  cudaDeviceGetAttribute( &value, cudaDevAttrMaxSharedMemoryPerBlockOptin, device );
  maxSharedMemoryPerBlockOptin = size_t(value);

  // printf("   Max shared memory per block:            %zu bytes\n"
  //        "   Max shared memory per block (Opt in):   %zu bytes\n",
  //         maxSharedMemoryPerBlock,
  //         maxSharedMemoryPerBlockOptin);

  // set up device no_co arrays
  int dim1 = 1;
  for( int i=0; i<p_x_dim; i++) {
    dim1 *= 3;
  }

  t_f1<int> h_neighbor_til_id = t_f1<int>( dim1, n_tils, neighbor_til_id, false );
  chkErr(cudaMalloc( reinterpret_cast<void **>(&d_neighbor_til_id), sizeof(t_f1<int>) ));
  h_neighbor_til_id.malloc_copy_to_device( d_neighbor_til_id );

  t_f2<int> h_shift = t_f2<int>( p_x_dim, dim1, n_tils, shift, false );
  chkErr(cudaMalloc( reinterpret_cast<void **>(&d_shift), sizeof(t_f2<int>) ));
  h_shift.malloc_copy_to_device( d_shift );


  // initialize kernel error handling member data
  h_err = (int *) malloc(sizeof(int));
  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_err), sizeof(int)));

  *h_err = 0;
  chkErr(cudaMemcpy( d_err, h_err, sizeof(int), cudaMemcpyHostToDevice));

  // initialize streams
  // streams with NonBlocking flag can execute concurrently with the default stream
  chkErr(cudaStreamCreateWithFlags( &s_j, cudaStreamNonBlocking ));
  chkErr(cudaStreamCreateWithFlags( &s_adv_dep, cudaStreamNonBlocking ));

}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_node_conf_c::init_dlb( int n_tils_, int my_aid_, int no_num, int *neighbor_til_id,
                              int *shift ) {

  n_tils = n_tils_;

  // set up device no_co arrays
  int dim1 = 1;
  for( int i=0; i<p_x_dim; i++) {
    dim1 *= 3;
  }

  t_f1<int> h_neighbor_til_id = t_f1<int>( dim1, n_tils, neighbor_til_id, false );
  chkErr(cudaMalloc( reinterpret_cast<void **>(&d_neighbor_til_id), sizeof(t_f1<int>) ));
  h_neighbor_til_id.malloc_copy_to_device( d_neighbor_til_id );

  t_f2<int> h_shift = t_f2<int>( p_x_dim, dim1, n_tils, shift, false );
  chkErr(cudaMalloc( reinterpret_cast<void **>(&d_shift), sizeof(t_f2<int>) ));
  h_shift.malloc_copy_to_device( d_shift );

}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_node_conf_c::cleanup( void ) {
  d_neighbor_til_id->cleanup( true );
  chkErr(cudaFree(d_neighbor_til_id));
  d_shift->cleanup( true );
  chkErr(cudaFree(d_shift));
  free(h_err);
  chkErr(cudaFree(d_err));
  chkErr(cudaStreamDestroy(s_j));
  chkErr(cudaStreamDestroy(s_adv_dep));
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_node_conf_c::cleanup_dlb( void ) {
  d_neighbor_til_id->cleanup( true );
  chkErr(cudaFree(d_neighbor_til_id));
  d_shift->cleanup( true );
  chkErr(cudaFree(d_shift));
}
