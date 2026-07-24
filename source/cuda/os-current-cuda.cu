#include <stdio.h>
#include <algorithm>
using namespace std;

#include "os-current-cuda.h"
#include "os-sys-cuda.h"

t_current_c::t_current_c( void ) {
  nx = NULL;
  gc_num = NULL;
  h_jd1 = NULL;
  h_jd2 = NULL;
  h_jd3 = NULL;
  d_jd1 = NULL;
  d_jd2 = NULL;
  d_jd3 = NULL;
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_current_c::init_group( int n_tils ) {
  nx = new t_f1<int>( p_max_dim, n_tils, false );
  dims = new t_f1<int>( p_x_dim+1, n_tils, false );
  nx_tot_max = -1;
  gc_num = new t_f2<int>( p_upper+1, p_max_dim, n_tils, false );
  h_j_data_arr = new p_k_fld*[n_tils];
  prefix_size = new int[n_tils+2];

  switch(p_x_dim) {
    case 1:
      h_jd1 = (t_f1<p_k_fld>*) malloc(n_tils * sizeof (t_f1<p_k_fld>));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_jd1), sizeof(t_f1<p_k_fld>)*n_tils));
      break;
    case 2:
      h_jd2 = (t_f2<p_k_fld>*) malloc(n_tils * sizeof (t_f2<p_k_fld>));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_jd2), sizeof(t_f2<p_k_fld>)*n_tils));
      break;
    case 3:
      h_jd3 = (t_f3<p_k_fld>*) malloc(n_tils * sizeof (t_f3<p_k_fld>));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_jd3), sizeof(t_f3<p_k_fld>)*n_tils));
      break;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_current_c::init_tile( int til, int *nx_, int *gc_num_, p_k_fld *jay,
                             int n_tils, t_mem_transfer_manager *mem_transfer_manager,
                             bool if_last_til ) {

  // Fill in variables for tile

  // tot number of cells in this tile
  int nx_tot = p_f_dim;

  dims->fh(0,til) = p_f_dim;
  for (int i=0; i<p_x_dim; ++i) {
    nx->fh(i,til) = nx_[i];
    gc_num->fh(p_lower,i,til) = gc_num_[i*(p_upper+1) + p_lower];
    gc_num->fh(p_upper,i,til) = gc_num_[i*(p_upper+1) + p_upper];
    dims->fh(i+1,til) = nx_[i] + gc_num_[i*(p_upper+1) + p_lower]
                              + gc_num_[i*(p_upper+1) + p_upper];

    nx_tot *= dims->fh(i+1,til);
  }

  if ( nx_tot > nx_tot_max ) {
    nx_tot_max = nx_tot;
  }

  // Allocate field array objects for this tile
  // Point buffers to the first cell that is not a guard cell
  // Malloc and copy device buffers
  switch(p_x_dim) {
    case 1:
      new (h_jd1+til) t_f1<p_k_fld>( p_f_dim,
                          nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                          jay, false );

      prefix_size[til+1] = h_jd1[til].size();

      (h_jd1+til)->offset( 0, gc_num->fh(p_lower,0,til) );

      h_j_data_arr[til] = h_jd1[til].buffered_malloc_copy_to_device( d_jd1+til, mem_transfer_manager, false );
      break;
    case 2:
      new (h_jd2+til) t_f2<p_k_fld>( p_f_dim,
                          nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                          nx->fh(1,til)+gc_num->fh(p_lower,1,til)+gc_num->fh(p_upper,1,til),
                          jay, false );

      prefix_size[til+1] = h_jd2[til].size();

      (h_jd2+til)->offset( 0, gc_num->fh(p_lower,0,til), gc_num->fh(p_lower,1,til) );

      h_j_data_arr[til] = h_jd2[til].buffered_malloc_copy_to_device( d_jd2+til, mem_transfer_manager, false );
      break;
    case 3:
      new (h_jd3+til) t_f3<p_k_fld>( p_f_dim,
                          nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                          nx->fh(1,til)+gc_num->fh(p_lower,1,til)+gc_num->fh(p_upper,1,til),
                          nx->fh(2,til)+gc_num->fh(p_lower,2,til)+gc_num->fh(p_upper,2,til),
                          jay, false );

      prefix_size[til+1] = h_jd3[til].size();

      (h_jd3+til)->offset( 0, gc_num->fh(p_lower,0,til), gc_num->fh(p_lower,1,til),
                         gc_num->fh(p_lower,2,til) );

      h_j_data_arr[til] = h_jd3[til].buffered_malloc_copy_to_device( d_jd3+til, mem_transfer_manager, false );
      break;
  }

  if ( til == n_tils-1 or if_last_til ) {
    init_tile_last( n_tils, mem_transfer_manager );
  }
}


// ---------------------------------------------------------------------------------------
// Complete initialization for tiles
// ---------------------------------------------------------------------------------------
void t_current_c::init_tile_last( int n_tils,
                                  t_mem_transfer_manager *mem_transfer_manager ) {

  // Compute prefix_scan of number of blocks required for each tile
  prefix_size[0] = 0;
  for (int i=1; i<n_tils+1; i++) {
    prefix_size[i] = (prefix_size[i] + __BLOCK_SIZE__ - 1) / __BLOCK_SIZE__ + prefix_size[i-1];
  }

  // Store average (ceiling) in last slot
  prefix_size[n_tils+1] = (prefix_size[n_tils] + n_tils - 1) / n_tils;

  chkErr(cudaMalloc(reinterpret_cast<void **>(&d_prefix_size), (n_tils+2) * sizeof(int) ));
  mem_transfer_manager->bufferedCudaMemcpy( d_prefix_size, prefix_size, (n_tils+2) * sizeof(int),
                                            cudaMemcpyHostToDevice );

}


// ---------------------------------------------------------------------------------------
// Copy pointers to host data for stationary tiles to the temp current object, which will
// ultimately replace the old current object
// Copying of pointers to device data is done later
// ---------------------------------------------------------------------------------------
void t_current_c::init_tile_dlb( int til, int til_old, t_current_c *jay ){

  int nx_tot = p_f_dim; // tot number of cells in this tile

  dims->fh(0,til) = jay->dims->fh(0,til_old);
  for (int i=0; i<p_x_dim; ++i) {
    nx->fh(i,til)             = jay->nx->fh(i,til_old);
    gc_num->fh(p_lower,i,til) = jay->gc_num->fh(p_lower,i,til_old);
    gc_num->fh(p_upper,i,til) = jay->gc_num->fh(p_upper,i,til_old);
    dims->fh(i+1,til)         = jay->dims->fh(i+1,til_old);

    nx_tot *= dims->fh(i+1,til);
  }

  if ( nx_tot > nx_tot_max ) {
    nx_tot_max = nx_tot;
  }

  switch(p_x_dim) {
    case 1:
      h_jd1[til] = jay->h_jd1[til_old];
      prefix_size[til+1] = h_jd1[til].size();
      break;
    case 2:
      h_jd2[til] = jay->h_jd2[til_old];
      prefix_size[til+1] = h_jd2[til].size();
      break;
    case 3:
      h_jd3[til] = jay->h_jd3[til_old];
      prefix_size[til+1] = h_jd3[til].size();
      break;
  }

  h_j_data_arr[til]          = jay->h_j_data_arr[til_old];

  // Note: computing prefix scan is done later once we've recvd our new tils

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_current_c::buffered_copy_from_device( int til,
                                             t_mem_transfer_manager *mem_transfer_manager ){

  // copy host memory to device
  switch(p_x_dim) {
    case 1:
      h_jd1[til].buffered_copy_from_device( h_j_data_arr[til], mem_transfer_manager );
      break;
    case 2:
      h_jd2[til].buffered_copy_from_device( h_j_data_arr[til], mem_transfer_manager );
      break;
    case 3:
      h_jd3[til].buffered_copy_from_device( h_j_data_arr[til], mem_transfer_manager );
      break;
  }

}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_current_c::cleanup( int n_tils ) {

  nx->cleanup();
  dims->cleanup();
  gc_num->cleanup();

  delete nx;
  delete dims;
  delete gc_num;
  delete[] h_j_data_arr;
  delete[] prefix_size;
  chkErr(cudaFree(d_prefix_size));

  switch(p_x_dim) {
    case 1:
      for (int til=0; til<n_tils; ++til) {
        h_jd1[til].cleanup();
        h_jd1[til].~t_f1();
        (d_jd1+til)->cleanup( true );
      }
      free(h_jd1);
      chkErr(cudaFree(d_jd1));
      break;
    case 2:
      for (int til=0; til<n_tils; ++til) {
        h_jd2[til].cleanup();
        h_jd2[til].~t_f2();
        (d_jd2+til)->cleanup( true );
      }
      free(h_jd2);
      chkErr(cudaFree(d_jd2));
      break;
    case 3:
      for (int til=0; til<n_tils; ++til) {
        h_jd3[til].cleanup();
        h_jd3[til].~t_f3();
        (d_jd3+til)->cleanup( true );
      }
      free(h_jd3);
      chkErr(cudaFree(d_jd3));
      break;
  }

}


// ---------------------------------------------------------------------------------------
// cleanup objects that will no longer be needed post load balance
// ---------------------------------------------------------------------------------------
void t_current_c::cleanup_dlb( int *tils_send_old, int n_send ) {

  int til;

  nx->cleanup();
  dims->cleanup();
  gc_num->cleanup();

  delete nx;
  delete dims;
  delete gc_num;
  delete[] h_j_data_arr;
  delete[] prefix_size;
  chkErr(cudaFree(d_prefix_size));

  switch(p_x_dim) {
    case 1:
      for (int i=0; i<n_send; ++i) {
        til = tils_send_old[i];
        h_jd1[til].cleanup();
        h_jd1[til].~t_f1();
        (d_jd1+til)->cleanup( true );
      }
      free(h_jd1);
      chkErr(cudaFree(d_jd1));
      break;
    case 2:
      for (int i=0; i<n_send; ++i) {
        til = tils_send_old[i];
        h_jd2[til].cleanup();
        h_jd2[til].~t_f2();
        (d_jd2+til)->cleanup( true );
      }
      free(h_jd2);
      chkErr(cudaFree(d_jd2));
      break;
    case 3:
      for (int i=0; i<n_send; ++i) {
        til = tils_send_old[i];
        h_jd3[til].cleanup();
        h_jd3[til].~t_f3();
        (d_jd3+til)->cleanup( true );
      }
      free(h_jd3);
      chkErr(cudaFree(d_jd3));
      break;
  }

}
