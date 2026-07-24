#include <stdio.h>
#include <algorithm>
using namespace std;

#include "os-emf-cuda.h"
#include "os-sys-cuda.h"

t_emf_c::t_emf_c( void ) {
  nx = NULL;
  gc_num = NULL;
  h_ed1 = NULL; h_bd1 = NULL;
  h_ed2 = NULL; h_bd2 = NULL;
  h_ed3 = NULL; h_bd3 = NULL;
  d_ed1 = NULL; d_bd1 = NULL;
  d_ed2 = NULL; d_bd2 = NULL;
  d_ed3 = NULL; d_bd3 = NULL;
  h_e_data_arr = NULL; h_b_data_arr = NULL;
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_emf_c::init_group( int n_tils ) {
  nx = new t_f1<int>( p_max_dim, n_tils, false );
  dims = new t_f1<int>( p_x_dim+1, n_tils, false );
  nx_tot_max = -1;
  gc_num = new t_f2<int>( p_upper+1, p_max_dim, n_tils, false );
  h_e_data_arr = new p_k_fld*[n_tils];
  h_b_data_arr = new p_k_fld*[n_tils];

  switch(p_x_dim) {
    case 1:
      h_ed1 = (t_f1<p_k_fld>*) malloc(n_tils * sizeof (t_f1<p_k_fld>));
      h_bd1 = (t_f1<p_k_fld>*) malloc(n_tils * sizeof (t_f1<p_k_fld>));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_ed1), sizeof(t_f1<p_k_fld>)*n_tils));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_bd1), sizeof(t_f1<p_k_fld>)*n_tils));
      break;
    case 2:
      h_ed2 = (t_f2<p_k_fld>*) malloc(n_tils * sizeof (t_f2<p_k_fld>));
      h_bd2 = (t_f2<p_k_fld>*) malloc(n_tils * sizeof (t_f2<p_k_fld>));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_ed2), sizeof(t_f2<p_k_fld>)*n_tils));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_bd2), sizeof(t_f2<p_k_fld>)*n_tils));
      break;
    case 3:
      h_ed3 = (t_f3<p_k_fld>*) malloc(n_tils * sizeof (t_f3<p_k_fld>));
      h_bd3 = (t_f3<p_k_fld>*) malloc(n_tils * sizeof (t_f3<p_k_fld>));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_ed3), sizeof(t_f3<p_k_fld>)*n_tils));
      chkErr(cudaMalloc(reinterpret_cast<void **>(&d_bd3), sizeof(t_f3<p_k_fld>)*n_tils));
      break;
  }
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_emf_c::init_tile( int til, int *nx_, int *gc_num_, p_k_fld *e, p_k_fld *b,
                         t_mem_transfer_manager *mem_transfer_manager ) {


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
      new (h_ed1+til) t_f1<p_k_fld>( p_f_dim,
                           nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                           e, false );
      new (h_bd1+til) t_f1<p_k_fld>( p_f_dim,
                           nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                           b, false );

      h_ed1[til].offset( 0, gc_num->fh(p_lower,0,til) );
      h_bd1[til].offset( 0, gc_num->fh(p_lower,0,til) );

      h_e_data_arr[til] = h_ed1[til].buffered_malloc_copy_to_device( d_ed1+til, mem_transfer_manager, false );
      h_b_data_arr[til] = h_bd1[til].buffered_malloc_copy_to_device( d_bd1+til, mem_transfer_manager, false );
      break;
    case 2:
      new (h_ed2+til) t_f2<p_k_fld>( p_f_dim,
                           nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                           nx->fh(1,til)+gc_num->fh(p_lower,1,til)+gc_num->fh(p_upper,1,til),
                           e, false );
      new (h_bd2+til) t_f2<p_k_fld>( p_f_dim,
                           nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                           nx->fh(1,til)+gc_num->fh(p_lower,1,til)+gc_num->fh(p_upper,1,til),
                           b, false );

      h_ed2[til].offset( 0, gc_num->fh(p_lower,0,til), gc_num->fh(p_lower,1,til) );
      h_bd2[til].offset( 0, gc_num->fh(p_lower,0,til), gc_num->fh(p_lower,1,til) );

      h_e_data_arr[til] = h_ed2[til].buffered_malloc_copy_to_device( d_ed2+til, mem_transfer_manager, false );
      h_b_data_arr[til] = h_bd2[til].buffered_malloc_copy_to_device( d_bd2+til, mem_transfer_manager, false );
      break;
    case 3:
      new (h_ed3+til) t_f3<p_k_fld>( p_f_dim,
                           nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                           nx->fh(1,til)+gc_num->fh(p_lower,1,til)+gc_num->fh(p_upper,1,til),
                           nx->fh(2,til)+gc_num->fh(p_lower,2,til)+gc_num->fh(p_upper,2,til),
                           e, false );
      new (h_bd3+til) t_f3<p_k_fld>( p_f_dim,
                           nx->fh(0,til)+gc_num->fh(p_lower,0,til)+gc_num->fh(p_upper,0,til),
                           nx->fh(1,til)+gc_num->fh(p_lower,1,til)+gc_num->fh(p_upper,1,til),
                           nx->fh(2,til)+gc_num->fh(p_lower,2,til)+gc_num->fh(p_upper,2,til),
                           b, false );

      h_ed3[til].offset( 0, gc_num->fh(p_lower,0,til), gc_num->fh(p_lower,1,til),
                       gc_num->fh(p_lower,2,til) );
      h_bd3[til].offset( 0, gc_num->fh(p_lower,0,til), gc_num->fh(p_lower,1,til),
                       gc_num->fh(p_lower,2,til) );

      h_e_data_arr[til] = h_ed3[til].buffered_malloc_copy_to_device( d_ed3+til, mem_transfer_manager, false );
      h_b_data_arr[til] = h_bd3[til].buffered_malloc_copy_to_device( d_bd3+til, mem_transfer_manager, false );
      break;
  }
}


// ---------------------------------------------------------------------------------------
// Copy pointers to host data for stationary tiles to the temp emf object, which will
// ultimately replace the old emf object
// Copying of pointers to device data is done later
// ---------------------------------------------------------------------------------------
void t_emf_c::init_tile_dlb( int til, int til_old, t_emf_c *emf ) {

  int nx_tot = p_f_dim; // tot number of cells in this tile

  dims->fh(0,til) = emf->dims->fh(0,til_old);
  for (int i=0; i<p_x_dim; ++i) {
    nx->fh(i,til)             = emf->nx->fh(i,til_old);
    gc_num->fh(p_lower,i,til) = emf->gc_num->fh(p_lower,i,til_old);
    gc_num->fh(p_upper,i,til) = emf->gc_num->fh(p_upper,i,til_old);
    dims->fh(i+1,til)         = emf->dims->fh(i+1,til_old);

    nx_tot *= dims->fh(i+1,til);
  }

  if ( nx_tot > nx_tot_max ) {
    nx_tot_max = nx_tot;
  }

  switch(p_x_dim) {
    case 1:
      h_ed1[til] = emf->h_ed1[til_old];
      h_bd1[til] = emf->h_bd1[til_old];
      break;
    case 2:
      h_ed2[til] = emf->h_ed2[til_old];
      h_bd2[til] = emf->h_bd2[til_old];
      break;
    case 3:
      h_ed3[til] = emf->h_ed3[til_old];
      h_bd3[til] = emf->h_bd3[til_old];
      break;
  }

  h_e_data_arr[til] = emf->h_e_data_arr[til_old];
  h_b_data_arr[til] = emf->h_b_data_arr[til_old];

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_emf_c::buffered_copy_to_device( int til, t_mem_transfer_manager *mem_transfer_manager ) {

  // copy host memory to device
  switch(p_x_dim) {
    case 1:
      h_bd1[til].buffered_copy_to_device( h_b_data_arr[til], mem_transfer_manager );
      h_ed1[til].buffered_copy_to_device( h_e_data_arr[til], mem_transfer_manager );
      break;
    case 2:
      h_bd2[til].buffered_copy_to_device( h_b_data_arr[til], mem_transfer_manager );
      h_ed2[til].buffered_copy_to_device( h_e_data_arr[til], mem_transfer_manager );
      break;
    case 3:
      h_bd3[til].buffered_copy_to_device( h_b_data_arr[til], mem_transfer_manager );
      h_ed3[til].buffered_copy_to_device( h_e_data_arr[til], mem_transfer_manager );
      break;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_emf_c::buffered_copy_from_device( int til, t_mem_transfer_manager *mem_transfer_manager ) {

  // copy host memory to device
  switch(p_x_dim) {
    case 1:
      h_bd1[til].buffered_copy_from_device( h_b_data_arr[til], mem_transfer_manager );
      h_ed1[til].buffered_copy_from_device( h_e_data_arr[til], mem_transfer_manager );
      break;
    case 2:
      h_bd2[til].buffered_copy_from_device( h_b_data_arr[til], mem_transfer_manager );
      h_ed2[til].buffered_copy_from_device( h_e_data_arr[til], mem_transfer_manager );
      break;
    case 3:
      h_bd3[til].buffered_copy_from_device( h_b_data_arr[til], mem_transfer_manager );
      h_ed3[til].buffered_copy_from_device( h_e_data_arr[til], mem_transfer_manager );
      break;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_emf_c::cleanup( int n_tils ) {

  nx->cleanup();
  dims->cleanup();
  gc_num->cleanup();

  delete nx;
  delete dims;
  delete gc_num;
  delete[] h_e_data_arr;
  delete[] h_b_data_arr;

  switch(p_x_dim) {
    case 1:
      for (int til=0; til<n_tils; ++til) {
        h_ed1[til].cleanup();
        h_bd1[til].cleanup();
        h_ed1[til].~t_f1();
        h_bd1[til].~t_f1();
        (d_ed1+til)->cleanup( true );
        (d_bd1+til)->cleanup( true );
      }
      free(h_ed1);
      free(h_bd1);
      chkErr(cudaFree(d_ed1));
      chkErr(cudaFree(d_bd1));
      break;
    case 2:
      for (int til=0; til<n_tils; ++til) {
        h_ed2[til].cleanup();
        h_bd2[til].cleanup();
        h_ed2[til].~t_f2();
        h_bd2[til].~t_f2();
        (d_ed2+til)->cleanup( true );
        (d_bd2+til)->cleanup( true );
      }
      free(h_ed2);
      free(h_bd2);
      chkErr(cudaFree(d_ed2));
      chkErr(cudaFree(d_bd2));
      break;
    case 3:
      for (int til=0; til<n_tils; ++til) {
        h_ed3[til].cleanup();
        h_bd3[til].cleanup();
        h_ed3[til].~t_f3();
        h_bd3[til].~t_f3();
        (d_ed3+til)->cleanup( true );
        (d_bd3+til)->cleanup( true );
      }
      free(h_ed3);
      free(h_bd3);
      chkErr(cudaFree(d_ed3));
      chkErr(cudaFree(d_bd3));
      break;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void t_emf_c::cleanup_dlb( int *tils_send_old, int n_send ) {

  int til;

  nx->cleanup();
  dims->cleanup();
  gc_num->cleanup();

  delete nx;
  delete dims;
  delete gc_num;
  delete[] h_e_data_arr;
  delete[] h_b_data_arr;

  switch(p_x_dim) {
    case 1:
      for (int i=0; i<n_send; ++i) {
        til = tils_send_old[i];
        h_ed1[til].cleanup();
        h_bd1[til].cleanup();
        h_ed1[til].~t_f1();
        h_bd1[til].~t_f1();
        (d_ed1+til)->cleanup( true );
        (d_bd1+til)->cleanup( true );
      }
      free(h_ed1);
      free(h_bd1);
      chkErr(cudaFree(d_ed1));
      chkErr(cudaFree(d_bd1));
      break;
    case 2:
      for (int i=0; i<n_send; ++i) {
        til = tils_send_old[i];
        h_ed2[til].cleanup();
        h_bd2[til].cleanup();
        h_ed2[til].~t_f2();
        h_bd2[til].~t_f2();
        (d_ed2+til)->cleanup( true );
        (d_bd2+til)->cleanup( true );
      }
      free(h_ed2);
      free(h_bd2);
      chkErr(cudaFree(d_ed2));
      chkErr(cudaFree(d_bd2));
      break;
    case 3:
      for (int i=0; i<n_send; ++i) {
        til = tils_send_old[i];
        h_ed3[til].cleanup();
        h_bd3[til].cleanup();
        h_ed3[til].~t_f3();
        h_bd3[til].~t_f3();
        (d_ed3+til)->cleanup( true );
        (d_bd3+til)->cleanup( true );
      }
      free(h_ed3);
      free(h_bd3);
      chkErr(cudaFree(d_ed3));
      chkErr(cudaFree(d_bd3));
      break;
  }

}
