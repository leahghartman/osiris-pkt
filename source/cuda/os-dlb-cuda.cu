#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "os-dlb-cuda.h"

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void copy_tils_stnry_1d( int *tils_stnry, int *tils_stnry_old, int n_tils_stnry,
                                    t_f1<p_k_fld> *e_f1, t_f1<p_k_fld> *b_f1,
                                    t_f1<p_k_fld> *e_f1_temp, t_f1<p_k_fld> *b_f1_temp,
                                    t_f1<p_k_fld> *jay_f1, t_f1<p_k_fld> *jay_f1_temp ){

  int tx = threadIdx.x;
  int bx = blockIdx.x;

  int idx = bx*BLOCK_SIZE + tx;

  int til, til_old;

  if ( idx < n_tils_stnry ) {
    til     = tils_stnry[idx];
    til_old = tils_stnry_old[idx];

    jay_f1_temp[til] = jay_f1[til_old];
    e_f1_temp[til]   = e_f1[til_old];
    b_f1_temp[til]   = b_f1[til_old];
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void copy_tils_stnry_2d( int *tils_stnry, int *tils_stnry_old, int n_tils_stnry,
                                    t_f2<p_k_fld> *e_f2, t_f2<p_k_fld> *b_f2,
                                    t_f2<p_k_fld> *e_f2_temp, t_f2<p_k_fld> *b_f2_temp,
                                    t_f2<p_k_fld> *jay_f2, t_f2<p_k_fld> *jay_f2_temp ){

  int tx = threadIdx.x;
  int bx = blockIdx.x;

  int idx = bx*BLOCK_SIZE + tx;

  int til, til_old;

  if ( idx < n_tils_stnry ) {
    til     = tils_stnry[idx];
    til_old = tils_stnry_old[idx];

    jay_f2_temp[til] = jay_f2[til_old];
    e_f2_temp[til]   = e_f2[til_old];
    b_f2_temp[til]   = b_f2[til_old];
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void copy_tils_stnry_3d( int *tils_stnry, int *tils_stnry_old, int n_tils_stnry,
                                    t_f3<p_k_fld> *e_f3, t_f3<p_k_fld> *b_f3,
                                    t_f3<p_k_fld> *e_f3_temp, t_f3<p_k_fld> *b_f3_temp,
                                    t_f3<p_k_fld> *jay_f3, t_f3<p_k_fld> *jay_f3_temp ){

  int tx = threadIdx.x;
  int bx = blockIdx.x;

  int idx = bx*BLOCK_SIZE + tx;

  int til, til_old;

  if ( idx < n_tils_stnry ) {
    til     = tils_stnry[idx];
    til_old = tils_stnry_old[idx];

    jay_f3_temp[til] = jay_f3[til_old];
    e_f3_temp[til]   = e_f3[til_old];
    b_f3_temp[til]   = b_f3[til_old];
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void copy_tils_stnry_wrapper( int *tils_stnry, int *tils_stnry_old, int n_tils_stnry,
                              t_current_c *jay, t_current_c *jay_temp,
                              t_emf_c *emf, t_emf_c *emf_temp,
                              int nblocks ){

  switch( p_x_dim ) {
    case 1:
      copy_tils_stnry_1d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>> ( tils_stnry, tils_stnry_old,
                                                        n_tils_stnry, emf->d_ed1, emf->d_bd1,
                                                        emf_temp->d_ed1, emf_temp->d_bd1,
                                                        jay->d_jd1, jay_temp->d_jd1 );
      break;
    case 2:
      copy_tils_stnry_2d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>> ( tils_stnry, tils_stnry_old,
                                                         n_tils_stnry, emf->d_ed2, emf->d_bd2,
                                                         emf_temp->d_ed2, emf_temp->d_bd2,
                                                         jay->d_jd2, jay_temp->d_jd2 );
      break;
    case 3:
      copy_tils_stnry_3d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>> ( tils_stnry, tils_stnry_old,
                                                         n_tils_stnry, emf->d_ed3, emf->d_bd3,
                                                         emf_temp->d_ed3, emf_temp->d_bd3,
                                                         jay->d_jd3, jay_temp->d_jd3 );
      break;
  }

  chkLastErr("copy_tils_stnry");

}
