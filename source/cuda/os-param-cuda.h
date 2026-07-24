#ifndef OS_PARAM_CUDA_H
#define OS_PARAM_CUDA_H

// ----------------------
// Mirrors os-param.f90
// ----------------------
#include <stddef.h>
#include <stdint.h>

/* #define __DEBUG__ */

typedef char t_chunk;

#define __BLOCK_SIZE__ 512

// num spatial components
#ifndef P_X_DIM
  #error You must define P_X_DIM = 1, 2 or 3. Right now it is not even defined (or at least C cannot see it)
#endif

#if P_X_DIM == 1
  #define p_x_dim 1
#elif P_X_DIM == 2
  #define p_x_dim 2
#elif P_X_DIM == 3
  #define p_x_dim 3
#else
  #error You must define P_X_DIM = 1, 2 or 3. It is set to an invalid value.
#endif

// Field and particle quantities
#ifdef PRECISION_SINGLE
  #define p_k_fld  float
  #define p_k_part float
#elif PRECISION_DOUBLE
  #define p_k_fld  double
  #define p_k_part double
#else
  // Make sure this choice is the same as in os-param.f90
  #define p_k_fld  double
  #define p_k_part double
#endif

#define _0_5_p_k_fld  static_cast<p_k_fld >(0.5)
#define _1_0_p_k_fld  static_cast<p_k_fld >(1.0)
#define _2_0_p_k_fld  static_cast<p_k_fld >(2.0)
#define _0_5_p_k_part static_cast<p_k_part>(0.5)
#define _1_0_p_k_part static_cast<p_k_part>(1.0)
#define _2_0_p_k_part static_cast<p_k_part>(2.0)
#define _0_5_p_double static_cast<p_double>(0.5)
#define _1_0_p_double static_cast<p_double>(1.0)
#define _2_0_p_double static_cast<p_double>(2.0)

// Must be the same as those defined in Fortran!
#define p_cell_low  1 // Probably not needed unless in cylindrical geometry
#define p_cell_near 2 // Probably not needed unless in cylindrical geometry
#define p_linear    1
#define p_quadratic 2
#define p_cubic     3
#define p_quartic   4

// parameters describing coordinates/dimensions
#define p_p_dim 3 // num velocity/momentum components
#define p_f_dim 3 // no. of comp. of a vectorfield

// parameters to define boundary operations at internal boundaries
// boundary indexes
// NOTE: differ from those in os-param.f90 due to 0 indexing in C
#define p_lower 0
#define p_upper 1
#define p_size  2

// parameter used to tag holes in particle array
#define p_hole -600

// parameter used by t_mem_transfer_manager to determine copy/append. Has to be negative
#define p_none_present -1

#define p_err_nchunks_max_ext -10 // we ran out of space in ***chunks_ext_buffer
#define p_err_nchunks_max     -11 // we ran out of space in ***chunks buffer
#define p_err_nchunks_avail   -12 // we ran out of chunks in the pool
#define p_err_pierre_ihole    -13 // we ran out of space in ihole

#endif
