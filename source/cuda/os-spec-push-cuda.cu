// if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "os-spec-push-cuda.h"

// ---------------------------------------------------------------------------------------
// Atomic add for architectures which don't support atomic adds for floats
// ---------------------------------------------------------------------------------------
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
__device__ double atomicAdd(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

// ---------------------------------------------------------------------------------------
// Perform exclusive parallel prefix sum on p_double array idata which is of size BLOCK_SIZE
// Note: Could implement avoiding bank conficts if this kernel is heavy
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__device__ void prescan(p_double *const idata) {

  int tx = threadIdx.x;
  int offset = 1;

  // build sum in place up the tree
  for (int d = BLOCK_SIZE>>1; d > 0; d >>= 1) {
    __syncthreads();
    if (tx < d) {
      int ai = offset*(2*tx+1)-1;
      int bi = offset*(2*tx+2)-1;
      idata[bi] += idata[ai];
    }
    offset *= 2;
  }

  // clear the last element
  if (tx == 0) {
    idata[BLOCK_SIZE - 1] = 0;
  }

  // traverse down tree & build scan
  for (int d = 1; d < BLOCK_SIZE; d *= 2) {
    offset >>= 1;
    __syncthreads();
    if (tx < d) {
      int ai = offset*(2*tx+1)-1;
      int bi = offset*(2*tx+2)-1;
      p_double t = idata[ai];
      idata[ai] = idata[bi];
      idata[bi] += t;
    }
  }

}

//----------------------------------------------------------------------------------------
// Splines
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Linear
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s1( const p_k_fld x, p_k_fld *const s ) {

  s[0] = _0_5_p_k_fld - x;
  s[1] = _0_5_p_k_fld + x;

}

__forceinline__ __device__ void splineh_s1( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  s[0] = ( 1 - h ) - x;
  s[1] = (     h ) + x;

}

//----------------------------------------------------------------------------------------
// Quadratic
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s2( const p_k_fld x, p_k_fld *const s ) {

  const p_k_fld t0 = _0_5_p_k_fld - x;
  const p_k_fld t1 = _0_5_p_k_fld + x;

  s[-1] = _0_5_p_k_fld * t0*t0;
  s[ 0] = _0_5_p_k_fld + t0*t1;
  s[ 1] = _0_5_p_k_fld * t1*t1;

}

__forceinline__ __device__ void splineh_s2( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  const p_k_fld t0 = ( 1 - h ) - x;
  const p_k_fld t1 = (     h ) + x;

  s[-1] = _0_5_p_k_fld * t0*t0;
  s[ 0] = _0_5_p_k_fld + t0*t1;
  s[ 1] = _0_5_p_k_fld * t1*t1;

}

//----------------------------------------------------------------------------------------
// Cubic
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s3( const p_k_fld x, p_k_fld *const s ) {

  const p_k_fld c1_6 = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c2_3 = _2_0_p_k_fld / p_k_fld(3.0);

  p_k_fld t0 = _0_5_p_k_fld - x;
  p_k_fld t1 = _0_5_p_k_fld + x;

  const p_k_fld t2 = t0 * t0;
  const p_k_fld t3 = t1 * t1;

  t0 *= t2;
  t1 *= t3;

  s[-1] = c1_6 * t0;
  s[ 0] = ( c2_3 - t3 ) + _0_5_p_k_fld*t1;
  s[ 1] = ( c2_3 - t2 ) + _0_5_p_k_fld*t0;
  s[ 2] = c1_6 * t1;

}

__forceinline__ __device__ void splineh_s3( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  const p_k_fld c1_6 = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c2_3 = _2_0_p_k_fld / p_k_fld(3.0);

  p_k_fld t0 = ( 1 - h ) - x;
  p_k_fld t1 = (     h ) + x;

  const p_k_fld t2 = t0 * t0;
  const p_k_fld t3 = t1 * t1;

  t0 *= t2;
  t1 *= t3;

  s[-1] = c1_6 * t0;
  s[ 0] = ( c2_3 - t3 ) + _0_5_p_k_fld*t1;
  s[ 1] = ( c2_3 - t2 ) + _0_5_p_k_fld*t0;
  s[ 2] = c1_6 * t1;

}

//----------------------------------------------------------------------------------------
// Quartic
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s4( const p_k_fld x, p_k_fld *const s ) {

  const p_k_fld c1_6   = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c1_24  = _1_0_p_k_fld / p_k_fld(24.0);
  const p_k_fld c11_24 = p_k_fld(11.0) / p_k_fld(24.0);

  const p_k_fld t0 = _0_5_p_k_fld - x;
  const p_k_fld t1 = _0_5_p_k_fld + x;
  const p_k_fld t2 = t0 * t1;

  s[-2] = c1_24 * ((t0*t0)*(t0*t0));
  s[-1] = c1_6 * ( p_k_fld(0.25) + t0 * ( _1_0_p_k_fld + t0 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 0] = c11_24 + t2 * ( _0_5_p_k_fld + p_k_fld(0.25) * t2 );
  s[ 1] = c1_6 * ( p_k_fld(0.25) + t1 * ( _1_0_p_k_fld + t1 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 2] = c1_24 * ((t1*t1)*(t1*t1));

}

__forceinline__ __device__ void splineh_s4( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  const p_k_fld c1_6   = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c1_24  = _1_0_p_k_fld / p_k_fld(24.0);
  const p_k_fld c11_24 = p_k_fld(11.0) / p_k_fld(24.0);

  const p_k_fld t0 = ( 1 - h ) - x;
  const p_k_fld t1 = (     h ) + x;
  const p_k_fld t2 = t0 * t1;

  s[-2] = c1_24 * ((t0*t0)*(t0*t0));
  s[-1] = c1_6 * ( p_k_fld(0.25) + t0 * ( _1_0_p_k_fld + t0 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 0] = c11_24 + t2 * ( _0_5_p_k_fld + p_k_fld(0.25) * t2 );
  s[ 1] = c1_6 * ( p_k_fld(0.25) + t1 * ( _1_0_p_k_fld + t1 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 2] = c1_24 * ((t1*t1)*(t1*t1));

}

//----------------------------------------------------------------------------------------
// Signbit function
//----------------------------------------------------------------------------------------
__forceinline__ __device__ int signbit_custom( const p_k_fld x ) {

  if ( x < 0 ) {
    return 1;
  } else {
    return 0;
  }

}

//----------------------------------------------------------------------------------------
// Longitudinal current weights
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void wl_s1( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  wl[0] = qnx * ( x1 - x0 );

}

__forceinline__ __device__ void wl_s2( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  const p_k_fld d    = x1 - x0;
  const p_k_fld s1_2 = _0_5_p_k_fld - x0;
  const p_k_fld p1_2 = _0_5_p_k_fld + x0;

  wl[-1] = s1_2 - _0_5_p_k_fld * d;
  wl[ 0] = p1_2 + _0_5_p_k_fld * d;

  const p_k_fld n = qnx * d;
  wl[-1] *= n;
  wl[ 0]  *= n;

}

__forceinline__ __device__ void wl_s3( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);
  const p_k_fld d    = x1 - x0;
  const p_k_fld s1_2 = _0_5_p_k_fld - x0;
  const p_k_fld p1_2 = _0_5_p_k_fld + x0;
  const p_k_fld d_3  = c1_3 * d;

  wl[-1] =  _0_5_p_k_fld * ( s1_2*s1_2 - d * ( s1_2 - d_3 ) );
  wl[ 0] =  ( p_k_fld(0.75) - x0*x0 ) - d * ( x0 + d_3 );
  wl[ 1] =  _0_5_p_k_fld * ( p1_2*p1_2 + d * ( p1_2 + d_3 ) );

  const p_k_fld n = qnx * d;
  wl[-1] *= n;
  wl[ 0] *= n;
  wl[ 1] *= n;

}

__forceinline__ __device__ void wl_s4( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  const p_k_fld c1_6 = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);
  const p_k_fld c2_3 = _2_0_p_k_fld / p_k_fld(3.0);

  const p_k_fld d = x1 - x0;
  const p_k_fld s = _0_5_p_k_fld - x0;
  const p_k_fld p = _0_5_p_k_fld + x0;
  const p_k_fld t = x0 - c1_6;
  const p_k_fld u = x0 + c1_6;

  const p_k_fld s2 = s * s;
  const p_k_fld s3 = s2 * s;
  const p_k_fld p2 = p * p;
  const p_k_fld p3 = p2 * p;

  wl[-2] =  c1_6 * ( d * ( d * ( s - p_k_fld(0.25)*d ) - p_k_fld(1.5) * s2 ) + s3 );
  wl[-1] =    d * ( ( _0_5_p_k_fld * d ) * (t + p_k_fld(0.25)*d) + ( p_k_fld(0.75)*(t*t) - c1_3 ) ) +
            ( _0_5_p_k_fld * p3 + ( c2_3 - p2 ) );
  wl[ 0] =  - d * ( ( _0_5_p_k_fld * d ) * (u + p_k_fld(0.25)*d) + ( p_k_fld(0.75)*(u*u) - c1_3 ) ) +
            ( _0_5_p_k_fld * s3 + ( c2_3 - s2 ) );
  wl[ 1] =  c1_6 * ( d * ( d * ( p + p_k_fld(0.25)*d ) + p_k_fld(1.5) * p2 ) + p3 );

  const p_k_fld n = qnx * d;
  wl[-2] *= n;
  wl[-1] *= n;
  wl[ 0] *= n;
  wl[ 1] *= n;

}


//----------------------------------------------------------------------------------------
// Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5]
// range. This is the fastest implementation (twice as fast as a sequence of ifs) because
// the two if structures compile as conditional moves and can be processed independently.
// This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
// routine that would break for x = nearest( 0.5, -1.0 )
//----------------------------------------------------------------------------------------
__forceinline__ __device__ int ntrim( const p_k_part x ) {
  int a, b;

  if ( x < -.5 ) {
    a = -1;
  } else {
    a = 0;
  }

  if ( x >= .5 ) {
    b = +1;
  } else {
    b = 0;
  }

  return a + b;

}


//----------------------------------------------------------------------------------------
// Zero out the current array
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void zero_jay_1d( const int *const prefix_size, const int avg_block_tile, const t_f1<p_k_fld> *const jay_f1){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  __shared__ int til, ind_start, size;
  __shared__ t_f1<p_k_fld> jay_f1_til;

  if (tx==0) {
    til = bx / avg_block_tile; // floor, best guess, could be wrong

    // Find the correct tile for this block number
    if ( bx < prefix_size[til] ) {
      while( bx < prefix_size[til] ) {
        til--;
      }
    } else {
      while( bx >= prefix_size[til+1] ) {
        til++;
      }
    }

    ind_start = (bx-prefix_size[til]) * BLOCK_SIZE;

    jay_f1_til = jay_f1[til];

    size = jay_f1_til.size();

  }
  __syncthreads();

  if (tx+ind_start < size) jay_f1_til.f[tx+ind_start] = 0;

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void zero_jay_2d( const int *const prefix_size, const int avg_block_tile, const t_f2<p_k_fld> *const jay_f2){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  __shared__ int til, ind_start, size;
  __shared__ t_f2<p_k_fld> jay_f2_til;

  if (tx==0) {
    til = bx / avg_block_tile; // floor, best guess, could be wrong

    // Find the correct tile for this block number
    if ( bx < prefix_size[til] ) {
      while( bx < prefix_size[til] ) {
        til--;
      }
    } else {
      while( bx >= prefix_size[til+1] ) {
        til++;
      }
    }

    ind_start = (bx-prefix_size[til]) * BLOCK_SIZE;

    jay_f2_til = jay_f2[til];

    size = jay_f2_til.size();

  }
  __syncthreads();

  if (tx+ind_start < size) jay_f2_til.f[tx+ind_start] = 0;

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void zero_jay_3d( const int *const prefix_size, const int avg_block_tile, const t_f3<p_k_fld> *const jay_f3){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  __shared__ int til, ind_start, size;
  __shared__ t_f3<p_k_fld> jay_f3_til;

  if (tx==0) {
    til = bx / avg_block_tile; // floor, best guess, could be wrong

    // Find the correct tile for this block number
    if ( bx < prefix_size[til] ) {
      while( bx < prefix_size[til] ) {
        til--;
      }
    } else {
      while( bx >= prefix_size[til+1] ) {
        til++;
      }
    }

    ind_start = (bx-prefix_size[til]) * BLOCK_SIZE;

    jay_f3_til = jay_f3[til];

    size = jay_f3_til.size();

  }
  __syncthreads();

  if (tx+ind_start < size) jay_f3_til.f[tx+ind_start] = 0;

}


//----------------------------------------------------------------------------------------
// Function stub to call kernel
//----------------------------------------------------------------------------------------
void zero_jay_wrapper( const int *const prefix_size, const int avg_block_tile, 
                       const t_current_c *const jay, const int nblocks,
                       const cudaStream_t stream ){

  switch (p_x_dim) {
    case 1:
      zero_jay_1d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( prefix_size, avg_block_tile,
                                                              jay->d_jd1 );
      break;
    case 2:
      zero_jay_2d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( prefix_size, avg_block_tile,
                                                              jay->d_jd2 );
      break;
    case 3:
      zero_jay_3d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( prefix_size, avg_block_tile,
                                                              jay->d_jd3 );
      break;
  }

    chkLastErr("zero_jay");
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif

}


//----------------------------------------------------------------------------------------
// Shift particle positions for moving window
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void shift_positions( t_chunk ***chunks, int *block_start_tiles,
                                 int *block_start_chunks, int *block_num_par,
                                 int nmove0, int nmove1, int nmove2, int chunk_size ){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;

  if (tx==0) {
    til = block_start_tiles[bx];
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int nmove[3] = { nmove0, nmove1, nmove2 };
  int num_par_proc_chunk, part_idx;

  t_f1<int> ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t ix_offset = ( p_x_dim + p_p_dim + 1 ) * sizeof( p_k_part ) * chunk_size;

  while( num_par_proc_tot > 0 ) {

    ix.g = (int *) ( *chunk + ix_offset);

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      for (int i=0; i<p_x_dim; i++) {
        ix.fd(part_idx,i) -= nmove[i];
      }

      // Note: this only works if the chunks size is divisible by the block size (which we enforce upon initialization)
      // See sort for how to do it more generally if we decide we want to do that here too
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

}


//----------------------------------------------------------------------------------------
// Function stub to call kernel
//----------------------------------------------------------------------------------------
void shift_positions_wrapper( t_chunk ***chunks, int *block_start_tiles,
                              int *block_start_chunks, int *block_num_par,
                              int nmove0, int nmove1, int nmove2,
                              int nblocks, int chunk_size ){

  shift_positions<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>>( chunks, block_start_tiles,
                                                   block_start_chunks, block_num_par,
                                                   nmove0, nmove1, nmove2, chunk_size );

  chkLastErr("shift_positions");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}


// The following two procedures are kept (commented out) for performance debugging.
// These routines should give the peak possible performance. They manually unroll
// all loops, and remove all 0 additions in the loops. The downside to writing the deposition
// in this way is that the interpolation is handled explicitly, thus one would need to write 
// seperate kernels to handle each interpolation order. In practice, we found that loop unrolling 
// did not lead to any performance increase on NVIDIA A100 GPus. A performance increase was seen,
// but this seemed to be due to the removal of 0 additions, rather than due to the loop unrolling.
// Therefore, the kernels are written in such a way as to still abstract the interpolation order
// (see  DEPOSIT_CURRENT_ROLLED_2D) and still remove the 0 additions, thereby still getting the same 
// performance. 
////----------------------------------------------------------------------------------------
////----------------------------------------------------------------------------------------
//__forceinline__ __device__ void deposit_current_unrolled_2d_s1(
//        const int ix, const int iy,
//        const p_k_part x0, const p_k_part y0,
//        const p_k_part x1, const p_k_part y1, 
//        const p_k_part qnx, const p_k_part qny, const p_k_part qvz,
//        p_k_fld *j,
//        const int offset, const int n0, const int n1 )
//{
//    const p_k_fld S0x0 = _0_5_p_k_fld - x0;
//    const p_k_fld S0x1 = _0_5_p_k_fld + x0;

//    const p_k_fld S1x0 = _0_5_p_k_fld - x1;
//    const p_k_fld S1x1 = _0_5_p_k_fld + x1;

//    const p_k_fld S0y0 = _0_5_p_k_fld - y0;
//    const p_k_fld S0y1 = _0_5_p_k_fld + y0;

//    const p_k_fld S1y0 = _0_5_p_k_fld - y1;
//    const p_k_fld S1y1 = _0_5_p_k_fld + y1;

//    const p_k_fld wl1 = qnx * (x1 - x0);
//    const p_k_fld wl2 = qnx * (y1 - y0);
    
//    const p_k_fld wp10 = (S0y0 + S1y0);
//    const p_k_fld wp11 = (S0y1 + S1y1);
    
//    const p_k_fld wp20 = (S0x0 + S1x0);
//    const p_k_fld wp21 = (S0x1 + S1x1);

//    // old shared memory indexing
//    //     atomicAdd( &s[ offset + ( n1*(jj[i]+k2) + (ii[i]+k1) ) * n0 + 0 ], wl1[k1] * wp1[k2] );
//    //     atomicAdd( &s[ offset + ( n1*(jj[i]+k2) + (ii[i]+k1) ) * n0 + 1 ], wp2[k1] * wl2[k2] );
//    //     atomicAdd( &s[ offset + ( n1*(jj[i]+k2) + (ii[i]+k1) ) * n0 + 2 ], qvz * ( tmp1 + _0_5_p_k_fld * tmp2 ) );

//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix  ) ) * n0 + 0 ], wl1 * wp10 );
//    atomicAdd( &j[ offset + ( n1*(iy+1) + (ix  ) ) * n0 + 0 ], wl1 * wp11 );    

//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix  ) ) * n0 + 1 ], wl2 * wp20 );    
//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix+1) ) * n0 + 1 ], wl2 * wp21 );    

//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix  ) ) * n0 + 2 ],
//               qvz * ( S0x0 * S0y0 + S1x0 * S1y0 + _0_5_p_k_fld * (S0x0 * S1y0 + S1x0 * S0y0) ));
//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix+1) ) * n0 + 2 ],    
//               qvz * ( S0x1 * S0y0 + S1x1 * S1y0 + _0_5_p_k_fld * (S0x1 * S1y0 + S1x1 * S0y0) ));
//    atomicAdd( &j[ offset + ( n1*(iy+1) + (ix  ) ) * n0 + 2 ],
//               qvz * ( S0x0 * S0y1 + S1x0 * S1y1 + _0_5_p_k_fld * (S0x0 * S1y1 + S1x0 * S0y1) ));
//    atomicAdd( &j[ offset + ( n1*(iy+1) + (ix+1) ) * n0 + 2 ],               
//               qvz * ( S0x1 * S0y1 + S1x1 * S1y1 + _0_5_p_k_fld * (S0x1 * S1y1 + S1x1 * S0y1) ));

//}


////----------------------------------------------------------------------------------------
////----------------------------------------------------------------------------------------
//__forceinline__ __device__ void deposit_current_unrolled_2d_s2(
//        const int ix, const int iy,
//        const p_k_part x0, const p_k_part y0,
//        const p_k_part x1, const p_k_part y1, 
//        const p_k_part qnx, const p_k_part qny, const p_k_part qvz,
//        p_k_fld *j,
//        const int offset, const int n0, const int n1 )
//{
//    const p_k_fld S0x0 = _0_5_p_k_fld - x0;
//    const p_k_fld S0x1 = _0_5_p_k_fld + x0;

//    const p_k_fld S1x0 = _0_5_p_k_fld - x1;
//    const p_k_fld S1x1 = _0_5_p_k_fld + x1;

//    const p_k_fld S0y0 = _0_5_p_k_fld - y0;
//    const p_k_fld S0y1 = _0_5_p_k_fld + y0;

//    const p_k_fld S1y0 = _0_5_p_k_fld - y1;
//    const p_k_fld S1y1 = _0_5_p_k_fld + y1;

//    const p_k_fld wl1 = qnx * (x1 - x0);
//    const p_k_fld wl2 = qnx * (y1 - y0);
    
//    const p_k_fld wp10 = (S0y0 + S1y0);
//    const p_k_fld wp11 = (S0y1 + S1y1);
    
//    const p_k_fld wp20 = (S0x0 + S1x0);
//    const p_k_fld wp21 = (S0x1 + S1x1);

//    // old shared memory indexing
//    //     atomicAdd( &s[ offset + ( n1*(jj[i]+k2) + (ii[i]+k1) ) * n0 + 0 ], wl1[k1] * wp1[k2] );
//    //     atomicAdd( &s[ offset + ( n1*(jj[i]+k2) + (ii[i]+k1) ) * n0 + 1 ], wp2[k1] * wl2[k2] );
//    //     atomicAdd( &s[ offset + ( n1*(jj[i]+k2) + (ii[i]+k1) ) * n0 + 2 ], qvz * ( tmp1 + _0_5_p_k_fld * tmp2 ) );

//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix  ) ) * n0 + 0 ], wl1 * wp10 );
//    atomicAdd( &j[ offset + ( n1*(iy+1) + (ix  ) ) * n0 + 0 ], wl1 * wp11 );    

//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix  ) ) * n0 + 1 ], wl2 * wp20 );    
//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix+1) ) * n0 + 1 ], wl2 * wp21 );    

//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix  ) ) * n0 + 2 ],
//               qvz * ( S0x0 * S0y0 + S1x0 * S1y0 + _0_5_p_k_fld * (S0x0 * S1y0 + S1x0 * S0y0) ));
//    atomicAdd( &j[ offset + ( n1*(iy  ) + (ix+1) ) * n0 + 2 ],    
//               qvz * ( S0x1 * S0y0 + S1x1 * S1y0 + _0_5_p_k_fld * (S0x1 * S1y0 + S1x1 * S0y0) ));
//    atomicAdd( &j[ offset + ( n1*(iy+1) + (ix  ) ) * n0 + 2 ],
//               qvz * ( S0x0 * S0y1 + S1x0 * S1y1 + _0_5_p_k_fld * (S0x0 * S1y1 + S1x0 * S0y1) ));
//    atomicAdd( &j[ offset + ( n1*(iy+1) + (ix+1) ) * n0 + 2 ],               
//               qvz * ( S0x1 * S0y1 + S1x1 * S1y1 + _0_5_p_k_fld * (S0x1 * S1y1 + S1x1 * S0y1) ));

//}


#define __TEMPLATE__

//********************************* Linear interpolation *********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s1
#define DUDT_BORIS_2D   dudt_boris_2d_s1
#define DUDT_BORIS_3D   dudt_boris_3d_s1

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs1
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs1
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs1

#define SPLINE spline_s1
#define SPLINEH splineh_s1

// Lower/upper point limits for dudt
#define LP 0
#define UP 1

#define WL wl_s1

// Perpendicular current cell  limits
#define P0 0
#define P1 1

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s1
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s1
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s1

// Rolled deposit_current functions which are compiled for all dimensions
#define DEPOSIT_CURRENT_ROLLED_1D   deposit_current_rolled_1d_s1
#define DEPOSIT_CURRENT_ROLLED_2D   deposit_current_rolled_2d_s1
#define DEPOSIT_CURRENT_ROLLED_3D   deposit_current_rolled_3d_s1

// This is the current deposition that is actually called by advance_deposit_*d.
// One can optionally change this (e.g. the hard-coded routine deposit_current_unrolled_*d) 
// to use a custom current deposition function instead of the default rolled current deposition.
#define DEPOSIT_CURRENT_1D   deposit_current_rolled_1d_s1
#define DEPOSIT_CURRENT_2D   deposit_current_rolled_2d_s1
#define DEPOSIT_CURRENT_3D   deposit_current_rolled_3d_s1

#include __FILE__

//******************************* Quadratic interpolation ********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s2
#define DUDT_BORIS_2D   dudt_boris_2d_s2
#define DUDT_BORIS_3D   dudt_boris_3d_s2

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs2
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs2
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs2

#define SPLINE spline_s2
#define SPLINEH splineh_s2

// Lower/upper point limits for dudt
#define LP -1
#define UP 1

#define WL wl_s2

// Perpendicular current cell  limits
#define P0 -1
#define P1 1

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s2
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s2
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s2

// Rolled deposit_current functions which are compiled for all dimensions
#define DEPOSIT_CURRENT_ROLLED_1D   deposit_current_rolled_1d_s2
#define DEPOSIT_CURRENT_ROLLED_2D   deposit_current_rolled_2d_s2
#define DEPOSIT_CURRENT_ROLLED_3D   deposit_current_rolled_3d_s2

// This is the current deposition that is actually called by advance_deposit_*d.
// One can optionally change this (e.g. the hard-coded routine deposit_current_unrolled_*d) 
// to use a custom current deposition function instead of the default rolled current deposition.
#define DEPOSIT_CURRENT_1D   deposit_current_rolled_1d_s2
#define DEPOSIT_CURRENT_2D   deposit_current_rolled_2d_s2
#define DEPOSIT_CURRENT_3D   deposit_current_rolled_3d_s2

#include __FILE__

//********************************* Cubic interpolation **********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s3
#define DUDT_BORIS_2D   dudt_boris_2d_s3
#define DUDT_BORIS_3D   dudt_boris_3d_s3

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs3
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs3
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs3

#define SPLINE spline_s3
#define SPLINEH splineh_s3

// Lower/upper point limits for dudt
#define LP -1
#define UP 2

#define WL wl_s3

// Perpendicular current cell  limits
#define P0 -1
#define P1 2

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s3
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s3
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s3

// Rolled deposit_current functions which are compiled for all dimensions
#define DEPOSIT_CURRENT_ROLLED_1D   deposit_current_rolled_1d_s3
#define DEPOSIT_CURRENT_ROLLED_2D   deposit_current_rolled_2d_s3
#define DEPOSIT_CURRENT_ROLLED_3D   deposit_current_rolled_3d_s3

// This is the current deposition that is actually called by advance_deposit_*d.
// One can optionally change this (e.g. the hard-coded routine deposit_current_unrolled_*d) 
// to use a custom current deposition function instead of the default rolled current deposition.
#define DEPOSIT_CURRENT_1D   deposit_current_rolled_1d_s3
#define DEPOSIT_CURRENT_2D   deposit_current_rolled_2d_s3
#define DEPOSIT_CURRENT_3D   deposit_current_rolled_3d_s3

#include __FILE__

//******************************** Quartic interpolation *********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s4
#define DUDT_BORIS_2D   dudt_boris_2d_s4
#define DUDT_BORIS_3D   dudt_boris_3d_s4

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs4
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs4
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs4

#define SPLINE spline_s4
#define SPLINEH splineh_s4

// Lower/upper point limits for dudt
#define LP -2
#define UP 2

#define WL wl_s4

// Perpendicular current cell  limits
#define P0 -2
#define P1 2

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s4
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s4
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s4

// Rolled deposit_current functions which are compiled for all dimensions
#define DEPOSIT_CURRENT_ROLLED_1D   deposit_current_rolled_1d_s4
#define DEPOSIT_CURRENT_ROLLED_2D   deposit_current_rolled_2d_s4
#define DEPOSIT_CURRENT_ROLLED_3D   deposit_current_rolled_3d_s4

// This is the current deposition that is actually called by advance_deposit_*d.
// One can optionally change this (e.g. the hard-coded routine deposit_current_unrolled_*d) 
// to use a custom current deposition function instead of the default rolled current deposition.
#define DEPOSIT_CURRENT_1D   deposit_current_rolled_1d_s4
#define DEPOSIT_CURRENT_2D   deposit_current_rolled_2d_s4
#define DEPOSIT_CURRENT_3D   deposit_current_rolled_3d_s4

#include __FILE__

#undef __TEMPLATE__

#else


// // ---------------------------------------------------------------------------------------
// // NON-SHARED MEMORY VERSION OF DUDT
// // ---------------------------------------------------------------------------------------
// template <int BLOCK_SIZE>
// __global__ void DUDT_BORIS_1D( t_chunk ***chunks, int *block_start_tiles,
//                                int *block_start_chunks, int *block_num_par,
//                                p_k_part rqm, t_f1<p_k_fld> *b_f1,
//                                t_f1<p_k_fld> *e_f1, p_double *energy,
//                                bool report_energy, p_double dt, int chunk_size ){

//   const int bx = blockIdx.x;
//   const int tx = threadIdx.x;

//   // Shared variables that are the same between all threads
//   __shared__ int til;
//   __shared__ p_k_part tem;
//   __shared__ t_f1<p_k_fld> b_f1_til, e_f1_til;

//   // Shared variables that will be prefix scanned
//   __shared__ p_double accum_energy[BLOCK_SIZE];

//   if (tx==0) {
//     til = block_start_tiles[bx];

//     b_f1_til = b_f1[til];
//     e_f1_til = e_f1[til];

//     // get factor that includes timestep and charge-to-mass ratio
//     // note that the charge-to-mass ratio is used since the
//     // momentum p is not the total momentum of a particles but
//     // the momentum per unit restmass (which is the electron mass)
//     tem = p_k_part( _0_5_p_double * dt / rqm );
//   }
//   __syncthreads();

//   int num_par_proc_tot = block_num_par[bx];
//   t_chunk **chunk = chunks[til] + block_start_chunks[bx];

//   p_k_fld w1_buff[UP-LP+1], w1h_buff[UP-LP+1];

//   p_k_fld *const w1  = w1_buff  - LP;
//   p_k_fld *const w1h = w1h_buff - LP;

//   t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
//   t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
//   t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
//   t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

//   const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
//   const size_t x_offset = 0;
//   const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
//   const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
//   const size_t ix_offset = q_offset + chunk_size_t;

//   p_double accum_energy_last, loc_ene = 0;

//   // loop over all chunks assigned to this block 
//   while( num_par_proc_tot > 0 ) {

//     x.g  = (p_k_part *) ( *chunk + x_offset );
//     p.g  = (p_k_part *) ( *chunk + p_offset );
//     q.g  = (p_k_part *) ( *chunk + q_offset );
//     ix.g = (int *)      ( *chunk + ix_offset);

//     // modify bp & ep to include timestep and charge-to-mass ratio
//     // and perform half the electric field acceleration.

//     const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

//     // nloops is the number of times that we have repeated the inner while loop.
//     // this allows us to declare part_idx as a const. 
//     int nloops = 0; 

//     // loop over particles in a chunk, BLOCK_SIZE at a time
//     while( 1 ) {

//        const int part_idx = tx + nloops * BLOCK_SIZE;

//       __syncwarp();

//       if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

//       if( part_idx >= num_par_proc_chunk ) {
//         nloops++;
//         continue;
//       }

//       const int i1 = ix.fd(part_idx,0) - 1;

//       const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );

//       p_k_fld utemp2_0 = p.fd(part_idx,0);
//       p_k_fld utemp2_1 = p.fd(part_idx,1);
//       p_k_fld utemp2_2 = p.fd(part_idx,2);

//       const int h1 = signbit_custom(dx1);

//       const int i1h = i1 - h1;

//       // get spline weights for x
//       SPLINE( dx1, w1 );
//       SPLINEH( dx1, h1, w1h );

//       p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
      
//       // Interpolate E - Field
//       {
//         p_k_fld f1 = 0;
//         p_k_fld f2 = 0;
//         p_k_fld f3 = 0;

//         for (int k1=LP; k1<=UP; ++k1){
//           f1 += e_f1_til.fd(0,i1h+k1) * w1h[k1];
//           f2 += e_f1_til.fd(1,i1 +k1) * w1[ k1];
//           f3 += e_f1_til.fd(2,i1 +k1) * w1[ k1];
//         }
        
//         ep_0 = f1;
//         ep_1 = f2;
//         ep_2 = f3;
//       }
      
//       // Interpolate B - Field
//       {
//         p_k_fld f1 = 0;
//         p_k_fld f2 = 0;
//         p_k_fld f3 = 0;
        
//         for (int k1=LP; k1<=UP; ++k1){
//           f1 += b_f1_til.fd(0,i1 +k1) * w1[ k1];
//           f2 += b_f1_til.fd(1,i1h+k1) * w1h[k1];
//           f3 += b_f1_til.fd(2,i1h+k1) * w1h[k1];
//         }
        
//         bp_0 = f1;
//         bp_1 = f2;
//         bp_2 = f3;
//       }
      
//       ep_0 *= tem;
//       ep_1 *= tem;
//       ep_2 *= tem;
      
//       p_k_fld utemp_0 = utemp2_0 + ep_0;
//       p_k_fld utemp_1 = utemp2_1 + ep_1;
//       p_k_fld utemp_2 = utemp2_2 + ep_2;

//       // Get time centered gamma
//       const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

//       const p_k_part gamma = sqrt(u2+1);

//       // accumulate time centered energy
//       // this is done in double precision always
//       if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

//       const p_k_part gam_tem = tem / gamma;

//       bp_0 *= gam_tem;
//       bp_1 *= gam_tem;
//       bp_2 *= gam_tem;

//       utemp2_0 = utemp_0 + utemp_1 * bp_2;
//       utemp2_1 = utemp_1 + utemp_2 * bp_0;
//       utemp2_2 = utemp_2 + utemp_0 * bp_1;

//       utemp2_0 -= utemp_2 * bp_1;
//       utemp2_1 -= utemp_0 * bp_2;
//       utemp2_2 -= utemp_1 * bp_0;

//       const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

//       bp_0 *= otsq;
//       bp_1 *= otsq;
//       bp_2 *= otsq;

//       utemp_0 += utemp2_1 * bp_2;
//       utemp_1 += utemp2_2 * bp_0;
//       utemp_2 += utemp2_0 * bp_1;

//       utemp_0 -= utemp2_2 * bp_1;
//       utemp_1 -= utemp2_0 * bp_2;
//       utemp_2 -= utemp2_1 * bp_0;

//       // Perform second half of electric field acceleration.
//       p.fd(part_idx,0) = utemp_0 + ep_0;
//       p.fd(part_idx,1) = utemp_1 + ep_1;
//       p.fd(part_idx,2) = utemp_2 + ep_2;

//       nloops++;

//     } // end loop over chunk of particles ``while ( 1 )``

//     num_par_proc_tot -= num_par_proc_chunk;
//     chunk++;

//   } // end loop over chunks

//   if (report_energy) {
//     accum_energy[tx] = loc_ene;
//     if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
//     prescan<BLOCK_SIZE>( accum_energy );
//     __syncthreads();
//     if (tx==0) {
//       atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
//     }
//   }

// }


// ---------------------------------------------------------------------------------------
// Advance velocities and calculate time-centered total energy
// with fields on a staggered Yee mesh
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_1D( t_chunk ***chunks, int *block_start_tiles,
                               int *block_start_chunks, int *block_num_par,
                               p_k_part rqm, t_f1<p_k_fld> *b_f1,
                               t_f1<p_k_fld> *e_f1, p_double *energy,
                               bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til, n0, n1, offset;
  __shared__ p_k_part tem;
  __shared__ t_f1<p_k_fld> b_f1_til, e_f1_til;

  extern __shared__ p_k_fld s[];
  __shared__ p_k_fld *e_s;
  __shared__ p_k_fld *b_s;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f1_til = b_f1[til];
    e_f1_til = e_f1[til];

    n0 = b_f1_til.n0;
    n1 = b_f1_til.n1;
    offset = b_f1_til.g - b_f1_til.f;

    // divvy up shared mem array
    b_s = (p_k_fld *)s;
    e_s = (p_k_fld *)&b_s[n0*n1];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // populate shared mem arrays
  int i = tx;
  while ( i < n0 * n1 ) {
    b_s[i] = b_f1_til.f[i];
    e_s[i] = e_f1_til.f[i];
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_fld w1_buff[UP-LP+1], w1h_buff[UP-LP+1];

  p_k_fld *const w1  = w1_buff  - LP;
  p_k_fld *const w1h = w1h_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  // loop over all chunks assigned to this block 
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

    // nloops is the number of times that we have repeated the inner while loop.
    // this allows us to declare part_idx as a const. 
    int nloops = 0; 

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while( 1 ) {

       const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }

      const int i1 = ix.fd(part_idx,0) - 1;

      const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );

      p_k_fld utemp2_0 = p.fd(part_idx,0);
      p_k_fld utemp2_1 = p.fd(part_idx,1);
      p_k_fld utemp2_2 = p.fd(part_idx,2);

      const int h1 = signbit_custom(dx1);

      const int i1h = i1 - h1;

      // get spline weights for x
      SPLINE( dx1, w1 );
      SPLINEH( dx1, h1, w1h );

      p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
      
      // Interpolate E - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k1=LP; k1<=UP; ++k1){
          // keep these lines here for reference for the complicated indexing of shmem buff
          // f1 += e_f1_til.fd(0,i1h+k1) * w1h[k1];
          // f2 += e_f1_til.fd(1,i1 +k1) * w1[ k1];
          // f3 += e_f1_til.fd(2,i1 +k1) * w1[ k1];
          f1 += e_s[ offset + (i1h+k1)*n0 + 0 ] * w1h[k1];
          f2 += e_s[ offset + (i1 +k1)*n0 + 1 ] * w1[ k1];
          f3 += e_s[ offset + (i1 +k1)*n0 + 2 ] * w1[ k1];
        }
        
        ep_0 = f1;
        ep_1 = f2;
        ep_2 = f3;
      }
      
      // Interpolate B - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;
        
        for (int k1=LP; k1<=UP; ++k1){
          // keep these lines here for reference for the complicated indexing of shmem buff
          // f1 += b_f1_til.fd(0,i1 +k1) * w1[ k1];
          // f2 += b_f1_til.fd(1,i1h+k1) * w1h[k1];
          // f3 += b_f1_til.fd(2,i1h+k1) * w1h[k1];
          f1 += b_s[ offset + (i1 +k1)*n0 + 0 ] * w1[ k1];
          f2 += b_s[ offset + (i1h+k1)*n0 + 1 ] * w1h[k1];
          f3 += b_s[ offset + (i1h+k1)*n0 + 2 ] * w1h[k1];
        }
        
        bp_0 = f1;
        bp_1 = f2;
        bp_2 = f3;
      }
      
      ep_0 *= tem;
      ep_1 *= tem;
      ep_2 *= tem;
      
      p_k_fld utemp_0 = utemp2_0 + ep_0;
      p_k_fld utemp_1 = utemp2_1 + ep_1;
      p_k_fld utemp_2 = utemp2_2 + ep_2;

      // Get time centered gamma
      const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

      const p_k_part gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precision always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      const p_k_part gam_tem = tem / gamma;

      bp_0 *= gam_tem;
      bp_1 *= gam_tem;
      bp_2 *= gam_tem;

      utemp2_0 = utemp_0 + utemp_1 * bp_2;
      utemp2_1 = utemp_1 + utemp_2 * bp_0;
      utemp2_2 = utemp_2 + utemp_0 * bp_1;

      utemp2_0 -= utemp_2 * bp_1;
      utemp2_1 -= utemp_0 * bp_2;
      utemp2_2 -= utemp_1 * bp_0;

      const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

      bp_0 *= otsq;
      bp_1 *= otsq;
      bp_2 *= otsq;

      utemp_0 += utemp2_1 * bp_2;
      utemp_1 += utemp2_2 * bp_0;
      utemp_2 += utemp2_0 * bp_1;

      utemp_0 -= utemp2_2 * bp_1;
      utemp_1 -= utemp2_0 * bp_2;
      utemp_2 -= utemp2_1 * bp_0;

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp_0 + ep_0;
      p.fd(part_idx,1) = utemp_1 + ep_1;
      p.fd(part_idx,2) = utemp_2 + ep_2;

      nloops++;

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// // ---------------------------------------------------------------------------------------
// // NON-SHARED MEMORY VERSION OF DUDT
// // ---------------------------------------------------------------------------------------
// template <int BLOCK_SIZE>
// __global__ void DUDT_BORIS_2D( t_chunk ***chunks, int *block_start_tiles,
//                                int *block_start_chunks, int *block_num_par,
//                                p_k_part rqm, t_f2<p_k_fld> *b_f2,
//                                t_f2<p_k_fld> *e_f2, p_double *energy,
//                                bool report_energy, p_double dt, int chunk_size ){

//   const int bx = blockIdx.x;
//   const int tx = threadIdx.x;

//   __shared__ int til;
//   __shared__ p_k_part tem;
//   __shared__ t_f2<p_k_fld> b_f2_til, e_f2_til;

//   __shared__ p_double accum_energy[BLOCK_SIZE];

//   if (tx==0) {
//     til = block_start_tiles[bx];

//     b_f2_til = b_f2[til];
//     e_f2_til = e_f2[til];

//     tem = p_k_part( _0_5_p_double * dt / rqm );
//   }
//   __syncthreads();

//   int num_par_proc_tot = block_num_par[bx];
//   t_chunk **chunk = chunks[til] + block_start_chunks[bx];

//   p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1], w1h_buff[UP-LP+1], w2h_buff[UP-LP+1];

//   p_k_fld *const w1  = w1_buff  - LP;
//   p_k_fld *const w2  = w2_buff  - LP;
//   p_k_fld *const w1h = w1h_buff - LP;
//   p_k_fld *const w2h = w2h_buff - LP;

//   t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
//   t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
//   t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
//   t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

//   const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
//   const size_t x_offset = 0;
//   const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
//   const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
//   const size_t ix_offset = q_offset + chunk_size_t;

//   p_double accum_energy_last, loc_ene = 0;

//   // loop over all chunks assigned to this block
//   while( num_par_proc_tot > 0 ) {

//     x.g  = (p_k_part *) ( *chunk + x_offset );
//     p.g  = (p_k_part *) ( *chunk + p_offset );
//     q.g  = (p_k_part *) ( *chunk + q_offset );
//     ix.g = (int *)      ( *chunk + ix_offset);

//     const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

//     int nloops = 0;
    
//     // loop over particles in a chunk, BLOCK_SIZE at a time
//     while ( 1 ) {
      
//       const int part_idx = tx + nloops * BLOCK_SIZE;

//       __syncwarp();

//       if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

//       if( part_idx >= num_par_proc_chunk ) {
//         nloops++;
//         continue;
//       }
    
//       const int i1 = ix.fd(part_idx,0) - 1;
//       const int i2 = ix.fd(part_idx,1) - 1;

//       const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );
//       const p_k_part dx2 = p_k_fld( x.fd(part_idx,1) );

//       p_k_fld utemp2_0 = p.fd(part_idx,0);
//       p_k_fld utemp2_1 = p.fd(part_idx,1);
//       p_k_fld utemp2_2 = p.fd(part_idx,2);
    
//       const int h1 = signbit_custom(dx1);
//       const int h2 = signbit_custom(dx2);

//       const int i1h = i1 - h1;
//       const int i2h = i2 - h2;
    
//       SPLINE( dx1, w1 );
//       SPLINEH( dx1, h1, w1h );

//       SPLINE( dx2, w2 );
//       SPLINEH( dx2, h2, w2h );

//       p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
      
//       // Interpolate E - Field
//       {
//         p_k_fld f1 = 0;
//         p_k_fld f2 = 0;
//         p_k_fld f3 = 0;
        
//         for (int k2=LP; k2<=UP; ++k2){
//           p_k_fld f1line = 0;
//           p_k_fld f2line = 0;
//           p_k_fld f3line = 0;
          
//           for (int k1=LP; k1<=UP; ++k1){
//             f1line += e_f2_til.fd(0,i1h+k1,i2 +k2) * w1h[k1];
//             f2line += e_f2_til.fd(1,i1 +k1,i2h+k2) * w1[ k1];
//             f3line += e_f2_til.fd(2,i1 +k1,i2 +k2) * w1[ k1];
//           }
          
//           f1 += f1line * w2[ k2];
//           f2 += f2line * w2h[k2];
//           f3 += f3line * w2[ k2];
//         }
        
//         ep_0 = f1;
//         ep_1 = f2;
//         ep_2 = f3;
//       }
      
//       // Interpolate B - Field
//       {
//         p_k_fld f1 = 0;
//         p_k_fld f2 = 0;
//         p_k_fld f3 = 0;

//         for (int k2=LP; k2<=UP; ++k2){
//           p_k_fld f1line = 0;
//           p_k_fld f2line = 0;
//           p_k_fld f3line = 0;
          
//           for (int k1=LP; k1<=UP; ++k1){
//             f1line += b_f2_til.fd(0,i1 +k1,i2h+k2) * w1[ k1];
//             f2line += b_f2_til.fd(1,i1h+k1,i2 +k2) * w1h[k1];
//             f3line += b_f2_til.fd(2,i1h+k1,i2h+k2) * w1h[k1];
//           }
          
//           f1 += f1line * w2h[k2];
//           f2 += f2line * w2[ k2];
//           f3 += f3line * w2h[k2];
//         }

//         bp_0 = f1;
//         bp_1 = f2;
//         bp_2 = f3;
//       }
        
//       ep_0 *= tem;
//       ep_1 *= tem;
//       ep_2 *= tem;

//       p_k_fld utemp_0 = utemp2_0 + ep_0;
//       p_k_fld utemp_1 = utemp2_1 + ep_1;
//       p_k_fld utemp_2 = utemp2_2 + ep_2;

//       const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

//       const p_k_part gamma = sqrt(u2+1);

//       if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

//       const p_k_part gam_tem = tem / gamma;

//       bp_0 *= gam_tem;
//       bp_1 *= gam_tem;
//       bp_2 *= gam_tem;

//       utemp2_0 = utemp_0 + utemp_1 * bp_2;
//       utemp2_1 = utemp_1 + utemp_2 * bp_0;
//       utemp2_2 = utemp_2 + utemp_0 * bp_1;

//       utemp2_0 -= utemp_2 * bp_1;
//       utemp2_1 -= utemp_0 * bp_2;
//       utemp2_2 -= utemp_1 * bp_0;

//       const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

//       bp_0 *= otsq;
//       bp_1 *= otsq;
//       bp_2 *= otsq;

//       utemp_0 += utemp2_1 * bp_2;
//       utemp_1 += utemp2_2 * bp_0;
//       utemp_2 += utemp2_0 * bp_1;

//       utemp_0 -= utemp2_2 * bp_1;
//       utemp_1 -= utemp2_0 * bp_2;
//       utemp_2 -= utemp2_1 * bp_0;

//       p.fd(part_idx,0) = utemp_0 + ep_0;
//       p.fd(part_idx,1) = utemp_1 + ep_1;
//       p.fd(part_idx,2) = utemp_2 + ep_2;

//       nloops++;

//     } // end loop over chunk of particles ``while ( 1 )``

//     num_par_proc_tot -= num_par_proc_chunk;
//     chunk++;

//   } // end loop over chunks

//   if (report_energy) {
//     accum_energy[tx] = loc_ene;
//     if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
//     prescan<BLOCK_SIZE>( accum_energy );
//     __syncthreads();
//     if (tx==0) {
//       atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
//     }
//   }

// }


// ---------------------------------------------------------------------------------------
// See DUDT_BORIS_1D for implementation with comments 
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_2D( t_chunk ***chunks, int *block_start_tiles,
                               int *block_start_chunks, int *block_num_par,
                               p_k_part rqm, t_f2<p_k_fld> *b_f2,
                               t_f2<p_k_fld> *e_f2, p_double *energy,
                               bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  __shared__ int til, n0, n1, n2, offset;
  __shared__ p_k_part tem;
  __shared__ t_f2<p_k_fld> b_f2_til, e_f2_til;

  extern __shared__ p_k_fld s[];
  __shared__ p_k_fld *e_s;
  __shared__ p_k_fld *b_s;

  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f2_til = b_f2[til];
    e_f2_til = e_f2[til];

    n0 = b_f2_til.n0;
    n1 = b_f2_til.n1;
    n2 = b_f2_til.n2;
    offset = b_f2_til.g - b_f2_til.f;

    // divvy up shared mem array
    b_s = (p_k_fld *)s;
    e_s = (p_k_fld *)&b_s[n0*n1*n2];

    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // populate shared mem arrays
  int i = tx;
  while ( i < n0 * n1 * n2 ) {
    b_s[i] = b_f2_til.f[i];
    e_s[i] = e_f2_til.f[i];
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1], w1h_buff[UP-LP+1], w2h_buff[UP-LP+1];

  p_k_fld *const w1  = w1_buff  - LP;
  p_k_fld *const w2  = w2_buff  - LP;
  p_k_fld *const w1h = w1h_buff - LP;
  p_k_fld *const w2h = w2h_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  // loop over all chunks assigned to this block
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

    int nloops = 0;
    
    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }
    
      const int i1 = ix.fd(part_idx,0) - 1;
      const int i2 = ix.fd(part_idx,1) - 1;

      const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );
      const p_k_part dx2 = p_k_fld( x.fd(part_idx,1) );

      p_k_fld utemp2_0 = p.fd(part_idx,0);
      p_k_fld utemp2_1 = p.fd(part_idx,1);
      p_k_fld utemp2_2 = p.fd(part_idx,2);
    
      const int h1 = signbit_custom(dx1);
      const int h2 = signbit_custom(dx2);

      const int i1h = i1 - h1;
      const int i2h = i2 - h2;
    
      SPLINE( dx1, w1 );
      SPLINEH( dx1, h1, w1h );

      SPLINE( dx2, w2 );
      SPLINEH( dx2, h2, w2h );

      p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
      
      // Interpolate E - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;
        
        for (int k2=LP; k2<=UP; ++k2){
          p_k_fld f1line = 0;
          p_k_fld f2line = 0;
          p_k_fld f3line = 0;
          
          for (int k1=LP; k1<=UP; ++k1){
            // keep these lines here for reference for the complicated indexing of shmem buff
            // f1line += e_f2_til.fd(0,i1h+k1,i2 +k2) * w1h[k1];
            // f2line += e_f2_til.fd(1,i1 +k1,i2h+k2) * w1[ k1];
            // f3line += e_f2_til.fd(2,i1 +k1,i2 +k2) * w1[ k1];
            f1line += e_s[ offset + ( n1*(i2 +k2) + (i1h+k1) ) * n0 + 0 ] * w1h[k1];
            f2line += e_s[ offset + ( n1*(i2h+k2) + (i1 +k1) ) * n0 + 1 ] * w1[ k1];
            f3line += e_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 2 ] * w1[ k1];
          }
          
          f1 += f1line * w2[ k2];
          f2 += f2line * w2h[k2];
          f3 += f3line * w2[ k2];
        }
        
        ep_0 = f1;
        ep_1 = f2;
        ep_2 = f3;
      }
      
      // Interpolate B - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k2=LP; k2<=UP; ++k2){
          p_k_fld f1line = 0;
          p_k_fld f2line = 0;
          p_k_fld f3line = 0;
          
          for (int k1=LP; k1<=UP; ++k1){
            // keep these lines here for reference for the complicated indexing of shmem buff
            // f1line += b_f2_til.fd(0,i1 +k1,i2h+k2) * w1[ k1];
            // f2line += b_f2_til.fd(1,i1h+k1,i2 +k2) * w1h[k1];
            // f3line += b_f2_til.fd(2,i1h+k1,i2h+k2) * w1h[k1];
            f1line += b_s[ offset + ( n1*(i2h+k2) + (i1 +k1) ) * n0 + 0 ] * w1[ k1];
            f2line += b_s[ offset + ( n1*(i2 +k2) + (i1h+k1) ) * n0 + 1 ] * w1h[k1];
            f3line += b_s[ offset + ( n1*(i2h+k2) + (i1h+k1) ) * n0 + 2 ] * w1h[k1];
          }
          
          f1 += f1line * w2h[k2];
          f2 += f2line * w2[ k2];
          f3 += f3line * w2h[k2];
        }

        bp_0 = f1;
        bp_1 = f2;
        bp_2 = f3;
      }
        
      ep_0 *= tem;
      ep_1 *= tem;
      ep_2 *= tem;

      p_k_fld utemp_0 = utemp2_0 + ep_0;
      p_k_fld utemp_1 = utemp2_1 + ep_1;
      p_k_fld utemp_2 = utemp2_2 + ep_2;

      const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

      const p_k_part gamma = sqrt(u2+1);

      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      const p_k_part gam_tem = tem / gamma;

      bp_0 *= gam_tem;
      bp_1 *= gam_tem;
      bp_2 *= gam_tem;

      utemp2_0 = utemp_0 + utemp_1 * bp_2;
      utemp2_1 = utemp_1 + utemp_2 * bp_0;
      utemp2_2 = utemp_2 + utemp_0 * bp_1;

      utemp2_0 -= utemp_2 * bp_1;
      utemp2_1 -= utemp_0 * bp_2;
      utemp2_2 -= utemp_1 * bp_0;

      const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

      bp_0 *= otsq;
      bp_1 *= otsq;
      bp_2 *= otsq;

      utemp_0 += utemp2_1 * bp_2;
      utemp_1 += utemp2_2 * bp_0;
      utemp_2 += utemp2_0 * bp_1;

      utemp_0 -= utemp2_2 * bp_1;
      utemp_1 -= utemp2_0 * bp_2;
      utemp_2 -= utemp2_1 * bp_0;

      p.fd(part_idx,0) = utemp_0 + ep_0;
      p.fd(part_idx,1) = utemp_1 + ep_1;
      p.fd(part_idx,2) = utemp_2 + ep_2;

      nloops++;

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// // ---------------------------------------------------------------------------------------
// // NON-SHARED MEMORY VERSION OF DUDT
// // ---------------------------------------------------------------------------------------
// template <int BLOCK_SIZE>
// __global__ void DUDT_BORIS_3D( t_chunk ***chunks, int *block_start_tiles,
//                                int *block_start_chunks, int *block_num_par,
//                                p_k_part rqm, t_f3<p_k_fld> *b_f3,
//                                t_f3<p_k_fld> *e_f3, p_double *energy,
//                                bool report_energy, p_double dt, int chunk_size ){

//   const int bx = blockIdx.x;
//   const int tx = threadIdx.x;

//   __shared__ int til;
//   __shared__ p_k_part tem;
//   __shared__ t_f3<p_k_fld> b_f3_til, e_f3_til;

//   __shared__ p_double accum_energy[BLOCK_SIZE];

//   if (tx==0) {
//     til = block_start_tiles[bx];

//     b_f3_til = b_f3[til];
//     e_f3_til = e_f3[til];

//     tem = p_k_part( _0_5_p_double * dt / rqm );
//   }
//   __syncthreads();

//   int num_par_proc_tot = block_num_par[bx];
//   t_chunk **chunk = chunks[til] + block_start_chunks[bx];

//   p_k_fld  w1_buff[UP-LP+1],  w2_buff[UP-LP+1],  w3_buff[UP-LP+1];
//   p_k_fld w1h_buff[UP-LP+1], w2h_buff[UP-LP+1], w3h_buff[UP-LP+1];

//   p_k_fld *const w1  = w1_buff  - LP;
//   p_k_fld *const w2  = w2_buff  - LP;
//   p_k_fld *const w3  = w3_buff  - LP;
//   p_k_fld *const w1h = w1h_buff - LP;
//   p_k_fld *const w2h = w2h_buff - LP;
//   p_k_fld *const w3h = w3h_buff - LP;

//   t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
//   t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
//   t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
//   t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

//   const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
//   const size_t x_offset = 0;
//   const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
//   const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
//   const size_t ix_offset = q_offset + chunk_size_t;

//   p_double accum_energy_last, loc_ene = 0;

//   // loop over all chunks assigned to this block 
//   while( num_par_proc_tot > 0 ) {

//     x.g  = (p_k_part *) ( *chunk + x_offset );
//     p.g  = (p_k_part *) ( *chunk + p_offset );
//     q.g  = (p_k_part *) ( *chunk + q_offset );
//     ix.g = (int *)      ( *chunk + ix_offset);

//     const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

//     int nloops = 0; 

//     // loop over particles in a chunk, BLOCK_SIZE at a time
//     while ( 1 ) {
      
//       const int part_idx = tx + nloops * BLOCK_SIZE;

//       __syncwarp();

//       if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

//       if( part_idx >= num_par_proc_chunk ) {
//         nloops++;
//         continue;
//       }
      
//       const int i1 = ix.fd(part_idx,0) - 1;
//       const int i2 = ix.fd(part_idx,1) - 1;
//       const int i3 = ix.fd(part_idx,2) - 1;

//       const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );
//       const p_k_part dx2 = p_k_fld( x.fd(part_idx,1) );
//       const p_k_part dx3 = p_k_fld( x.fd(part_idx,2) );

//       p_k_fld utemp2_0 = p.fd(part_idx,0);
//       p_k_fld utemp2_1 = p.fd(part_idx,1);
//       p_k_fld utemp2_2 = p.fd(part_idx,2);

//       const int h1 = signbit_custom(dx1);
//       const int h2 = signbit_custom(dx2);
//       const int h3 = signbit_custom(dx3);

//       const int i1h = i1 - h1;
//       const int i2h = i2 - h2;
//       const int i3h = i3 - h3;

//       SPLINE( dx1, w1 );
//       SPLINEH( dx1, h1, w1h );

//       SPLINE( dx2, w2 );
//       SPLINEH( dx2, h2, w2h );

//       SPLINE( dx3, w3 );
//       SPLINEH( dx3, h3, w3h );

//       p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
            
//       // Interpolate E - Field
//       {
//         p_k_fld f1 = 0;
//         p_k_fld f2 = 0;
//         p_k_fld f3 = 0;

//         for (int k3=LP; k3<=UP; ++k3){
//           p_k_fld f1plane = 0;
//           p_k_fld f2plane = 0;
//           p_k_fld f3plane = 0;
          
//           for (int k2=LP; k2<=UP; ++k2){
//             p_k_fld f1line = 0;
//             p_k_fld f2line = 0;
//             p_k_fld f3line = 0;
            
//             for (int k1=LP; k1<=UP; ++k1){
//               f1line += e_f3_til.fd(0,i1h+k1,i2 +k2,(i3 +k3)) * w1h[k1];
//               f2line += e_f3_til.fd(1,i1 +k1,i2h+k2,(i3 +k3)) * w1[ k1];
//               f3line += e_f3_til.fd(2,i1 +k1,i2 +k2,(i3h+k3)) * w1[ k1];
//             }

//             f1plane += f1line * w2[ k2];
//             f2plane += f2line * w2h[k2];
//             f3plane += f3line * w2[ k2];
//           }
          
//           f1 += f1plane * w3[ k3];
//           f2 += f2plane * w3[ k3];
//           f3 += f3plane * w3h[k3];
          
//         }
        
//         ep_0 = f1;
//         ep_1 = f2;
//         ep_2 = f3;
//       }
    
//       // Interpolate B - Field
//       {
//         p_k_part f1 = 0;
//         p_k_part f2 = 0;
//         p_k_part f3 = 0;

//         for (int k3=LP; k3<=UP; ++k3){
//           p_k_part f1plane = 0;
//           p_k_part f2plane = 0;
//           p_k_part f3plane = 0;

//           for (int k2=LP; k2<=UP; ++k2){
//             p_k_part f1line = 0;
//             p_k_part f2line = 0;
//             p_k_part f3line = 0;

//             for (int k1=LP; k1<=UP; ++k1){
//               f1line += b_f3_til.fd(0,i1 +k1,i2h+k2,i3h+k3) * w1[ k1];
//               f2line += b_f3_til.fd(1,i1h+k1,i2 +k2,i3h+k3) * w1h[k1];
//               f3line += b_f3_til.fd(2,i1h+k1,i2h+k2,i3 +k3) * w1h[k1];
//             }

//             f1plane += f1line * w2h[k2];
//             f2plane += f2line * w2[ k2];
//             f3plane += f3line * w2h[k2];
//           }

//           f1 += f1plane * w3h[k3];
//           f2 += f2plane * w3h[k3];
//           f3 += f3plane * w3[ k3];

//         }

//         bp_0 = f1;
//         bp_1 = f2;
//         bp_2 = f3;
//       }

//       ep_0 *= tem;
//       ep_1 *= tem;
//       ep_2 *= tem;

//       p_k_fld utemp_0 = utemp2_0 + ep_0;
//       p_k_fld utemp_1 = utemp2_1 + ep_1;
//       p_k_fld utemp_2 = utemp2_2 + ep_2;

//       const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

//       const p_k_part gamma = sqrt(u2+1);

//       if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

//       const p_k_part gam_tem = tem / gamma;

//       bp_0 *= gam_tem;
//       bp_1 *= gam_tem;
//       bp_2 *= gam_tem;

//       utemp2_0 = utemp_0 + utemp_1 * bp_2;
//       utemp2_1 = utemp_1 + utemp_2 * bp_0;
//       utemp2_2 = utemp_2 + utemp_0 * bp_1;

//       utemp2_0 -= utemp_2 * bp_1;
//       utemp2_1 -= utemp_0 * bp_2;
//       utemp2_2 -= utemp_1 * bp_0;

//       const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

//       bp_0 *= otsq;
//       bp_1 *= otsq;
//       bp_2 *= otsq;

//       utemp_0 += utemp2_1 * bp_2;
//       utemp_1 += utemp2_2 * bp_0;
//       utemp_2 += utemp2_0 * bp_1;

//       utemp_0 -= utemp2_2 * bp_1;
//       utemp_1 -= utemp2_0 * bp_2;
//       utemp_2 -= utemp2_1 * bp_0;

//       p.fd(part_idx,0) = utemp_0 + ep_0;
//       p.fd(part_idx,1) = utemp_1 + ep_1;
//       p.fd(part_idx,2) = utemp_2 + ep_2;

//       nloops++;

//     } // end loop over chunk of particles ``while ( 1 )``

//     num_par_proc_tot -= num_par_proc_chunk;
//     chunk++;

//   } // end loop over chunks

//   if (report_energy) {
//     accum_energy[tx] = loc_ene;
//     if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
//     prescan<BLOCK_SIZE>( accum_energy );
//     __syncthreads();
//     if (tx==0) {
//       atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
//     }
//   }

// }


// ---------------------------------------------------------------------------------------
// See DUDT_BORIS_1D for implementation with comments 
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_3D( t_chunk ***chunks, int *block_start_tiles,
                               int *block_start_chunks, int *block_num_par,
                               p_k_part rqm, t_f3<p_k_fld> *b_f3,
                               t_f3<p_k_fld> *e_f3, p_double *energy,
                               bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  __shared__ int til, n0, n1, n2, n3, offset;
  __shared__ p_k_part tem;
  __shared__ t_f3<p_k_fld> b_f3_til, e_f3_til;

  extern __shared__ p_k_fld s[];
  __shared__ p_k_fld *e_s;
  __shared__ p_k_fld *b_s;

  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f3_til = b_f3[til];
    e_f3_til = e_f3[til];

    n0 = b_f3_til.n0;
    n1 = b_f3_til.n1;
    n2 = b_f3_til.n2;
    n3 = b_f3_til.n3;
    offset = b_f3_til.g - b_f3_til.f;

    // divvy up shared mem array
    b_s = (p_k_fld *)s;
    e_s = (p_k_fld *)&b_s[n0*n1*n2*n3];

    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // populate shared mem arrays
  int i = tx;
  while ( i < n0 * n1 * n2 * n3 ) {
    b_s[i] = b_f3_til.f[i];
    e_s[i] = e_f3_til.f[i];
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_fld  w1_buff[UP-LP+1],  w2_buff[UP-LP+1],  w3_buff[UP-LP+1];
  p_k_fld w1h_buff[UP-LP+1], w2h_buff[UP-LP+1], w3h_buff[UP-LP+1];

  p_k_fld *const w1  = w1_buff  - LP;
  p_k_fld *const w2  = w2_buff  - LP;
  p_k_fld *const w3  = w3_buff  - LP;
  p_k_fld *const w1h = w1h_buff - LP;
  p_k_fld *const w2h = w2h_buff - LP;
  p_k_fld *const w3h = w3h_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  // loop over all chunks assigned to this block 
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

    int nloops = 0; 

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }
      
      const int i1 = ix.fd(part_idx,0) - 1;
      const int i2 = ix.fd(part_idx,1) - 1;
      const int i3 = ix.fd(part_idx,2) - 1;

      const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );
      const p_k_part dx2 = p_k_fld( x.fd(part_idx,1) );
      const p_k_part dx3 = p_k_fld( x.fd(part_idx,2) );

      p_k_fld utemp2_0 = p.fd(part_idx,0);
      p_k_fld utemp2_1 = p.fd(part_idx,1);
      p_k_fld utemp2_2 = p.fd(part_idx,2);

      const int h1 = signbit_custom(dx1);
      const int h2 = signbit_custom(dx2);
      const int h3 = signbit_custom(dx3);

      const int i1h = i1 - h1;
      const int i2h = i2 - h2;
      const int i3h = i3 - h3;

      SPLINE( dx1, w1 );
      SPLINEH( dx1, h1, w1h );

      SPLINE( dx2, w2 );
      SPLINEH( dx2, h2, w2h );

      SPLINE( dx3, w3 );
      SPLINEH( dx3, h3, w3h );

      p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
            
      // Interpolate E - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k3=LP; k3<=UP; ++k3){
          p_k_fld f1plane = 0;
          p_k_fld f2plane = 0;
          p_k_fld f3plane = 0;
          
          for (int k2=LP; k2<=UP; ++k2){
            p_k_fld f1line = 0;
            p_k_fld f2line = 0;
            p_k_fld f3line = 0;
            
            for (int k1=LP; k1<=UP; ++k1){
              // f1line += e_f3_til.fd(0,i1h+k1,i2 +k2,(i3 +k3)) * w1h[k1];
              // f2line += e_f3_til.fd(1,i1 +k1,i2h+k2,(i3 +k3)) * w1[ k1];
              // f3line += e_f3_til.fd(2,i1 +k1,i2 +k2,(i3h+k3)) * w1[ k1];
              f1line += e_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1h+k1 ) * n0 + 0 ] * w1h[k1];
              f2line += e_s[ offset + ( ( n2*(i3 +k3) + i2h+k2 ) * n1 + i1 +k1 ) * n0 + 1 ] * w1[ k1];
              f3line += e_s[ offset + ( ( n2*(i3h+k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 2 ] * w1[ k1];
            }

            f1plane += f1line * w2[ k2];
            f2plane += f2line * w2h[k2];
            f3plane += f3line * w2[ k2];
          }
          
          f1 += f1plane * w3[ k3];
          f2 += f2plane * w3[ k3];
          f3 += f3plane * w3h[k3];
          
        }
        
        ep_0 = f1;
        ep_1 = f2;
        ep_2 = f3;
      }
    
      // Interpolate B - Field
      {
        p_k_part f1 = 0;
        p_k_part f2 = 0;
        p_k_part f3 = 0;

        for (int k3=LP; k3<=UP; ++k3){
          p_k_part f1plane = 0;
          p_k_part f2plane = 0;
          p_k_part f3plane = 0;

          for (int k2=LP; k2<=UP; ++k2){
            p_k_part f1line = 0;
            p_k_part f2line = 0;
            p_k_part f3line = 0;

            for (int k1=LP; k1<=UP; ++k1){
              // f1line += b_f3_til.fd(0,i1 +k1,i2h+k2,i3h+k3) * w1[ k1];
              // f2line += b_f3_til.fd(1,i1h+k1,i2 +k2,i3h+k3) * w1h[k1];
              // f3line += b_f3_til.fd(2,i1h+k1,i2h+k2,i3 +k3) * w1h[k1];
              f1line += b_s[ offset + ( ( n2*(i3h+k3) + i2h+k2 ) * n1 + i1 +k1 ) * n0 + 0 ] * w1[ k1];
              f2line += b_s[ offset + ( ( n2*(i3h+k3) + i2 +k2 ) * n1 + i1h+k1 ) * n0 + 1 ] * w1h[k1];
              f3line += b_s[ offset + ( ( n2*(i3 +k3) + i2h+k2 ) * n1 + i1h+k1 ) * n0 + 2 ] * w1h[k1];
            }

            f1plane += f1line * w2h[k2];
            f2plane += f2line * w2[ k2];
            f3plane += f3line * w2h[k2];
          }

          f1 += f1plane * w3h[k3];
          f2 += f2plane * w3h[k3];
          f3 += f3plane * w3[ k3];

        }

        bp_0 = f1;
        bp_1 = f2;
        bp_2 = f3;
      }

      ep_0 *= tem;
      ep_1 *= tem;
      ep_2 *= tem;

      p_k_fld utemp_0 = utemp2_0 + ep_0;
      p_k_fld utemp_1 = utemp2_1 + ep_1;
      p_k_fld utemp_2 = utemp2_2 + ep_2;

      const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

      const p_k_part gamma = sqrt(u2+1);

      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      const p_k_part gam_tem = tem / gamma;

      bp_0 *= gam_tem;
      bp_1 *= gam_tem;
      bp_2 *= gam_tem;

      utemp2_0 = utemp_0 + utemp_1 * bp_2;
      utemp2_1 = utemp_1 + utemp_2 * bp_0;
      utemp2_2 = utemp_2 + utemp_0 * bp_1;

      utemp2_0 -= utemp_2 * bp_1;
      utemp2_1 -= utemp_0 * bp_2;
      utemp2_2 -= utemp_1 * bp_0;

      const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

      bp_0 *= otsq;
      bp_1 *= otsq;
      bp_2 *= otsq;

      utemp_0 += utemp2_1 * bp_2;
      utemp_1 += utemp2_2 * bp_0;
      utemp_2 += utemp2_0 * bp_1;

      utemp_0 -= utemp2_2 * bp_1;
      utemp_1 -= utemp2_0 * bp_2;
      utemp_2 -= utemp2_1 * bp_0;

      p.fd(part_idx,0) = utemp_0 + ep_0;
      p.fd(part_idx,1) = utemp_1 + ep_1;
      p.fd(part_idx,2) = utemp_2 + ep_2;

      nloops++;

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// Advance velocities and calculate time-centered total energy
// with fields centered at cell corner
// See DUDT_BORIS_1D for implementation with comments for the parts that are similar 
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_1D_CS( t_chunk ***chunks, int *block_start_tiles,
                                  int *block_start_chunks, int *block_num_par,
                                  p_k_part rqm, t_f1<p_k_fld> *b_f1,
                                  t_f1<p_k_fld> *e_f1, p_double *energy,
                                  bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  __shared__ int til, n0, n1, offset;
  __shared__ p_k_part tem;
  __shared__ t_f1<p_k_fld> b_f1_til, e_f1_til;
  __shared__ p_double accum_energy[BLOCK_SIZE];

  extern __shared__ p_k_fld s[];
  __shared__ p_k_fld *e_s;
  __shared__ p_k_fld *b_s;

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f1_til = b_f1[til];
    e_f1_til = e_f1[til];

    n0 = b_f1_til.n0;
    n1 = b_f1_til.n1;
    offset = b_f1_til.g - b_f1_til.f;

    // divvy up shared mem array
    b_s = (p_k_fld *)s;
    e_s = (p_k_fld *)&b_s[n0*n1];

    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // populate shared mem arrays
  int i = tx;
  while ( i < n0 * n1 ) {
    b_s[i] = b_f1_til.f[i];
    e_s[i] = e_f1_til.f[i];
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_fld w1_buff[UP-LP+1];

  p_k_fld *const w1 = w1_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  
  // loop over all chunks assigned to this block 
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // loop over chunks of particles until all particles in tile have been processed
    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    int nloops = 0;

    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }
      
      const int i1 = ix.fd(part_idx,0) - 1;

      const p_k_fld dx1 = p_k_fld( x.fd(part_idx,0) );

      p_k_fld utemp2_0 = p.fd(part_idx,0);
      p_k_fld utemp2_1 = p.fd(part_idx,1);
      p_k_fld utemp2_2 = p.fd(part_idx,2);

      SPLINE( dx1, w1 );

      p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;

      // Interpolate E - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k1=LP; k1<=UP; ++k1){
          // f1 += e_f1_til.fd(0,i1+k1) * w1[k1];
          // f2 += e_f1_til.fd(1,i1+k1) * w1[k1];
          // f3 += e_f1_til.fd(2,i1+k1) * w1[k1];
          f1 += e_s[ offset + (i1 +k1)*n0 + 0 ] * w1[ k1];
          f2 += e_s[ offset + (i1 +k1)*n0 + 1 ] * w1[ k1];
          f3 += e_s[ offset + (i1 +k1)*n0 + 2 ] * w1[ k1];
        }

        ep_0 = f1;
        ep_1 = f2;
        ep_2 = f3;
      }
      
      // Interpolate B - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k1=LP; k1<=UP; ++k1){
          // f1 += b_f1_til.fd(0,i1+k1) * w1[k1];
          // f2 += b_f1_til.fd(1,i1+k1) * w1[k1];
          // f3 += b_f1_til.fd(2,i1+k1) * w1[k1];
          f1 += b_s[ offset + (i1 +k1)*n0 + 0 ] * w1[ k1];
          f2 += b_s[ offset + (i1 +k1)*n0 + 1 ] * w1[ k1];
          f3 += b_s[ offset + (i1 +k1)*n0 + 2 ] * w1[ k1];
        }

        bp_0 = f1;
        bp_1 = f2;
        bp_2 = f3;
      }

      ep_0 *= tem;
      ep_1 *= tem;
      ep_2 *= tem;

      p_k_fld utemp_0 = utemp2_0 + ep_0;
      p_k_fld utemp_1 = utemp2_1 + ep_1;
      p_k_fld utemp_2 = utemp2_2 + ep_2;

      const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

      const p_k_part gamma = sqrt(u2+1);

      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      const p_k_part gam_tem = tem / gamma;

      bp_0 *= gam_tem;
      bp_1 *= gam_tem;
      bp_2 *= gam_tem;

      utemp2_0 = utemp_0 + utemp_1 * bp_2;
      utemp2_1 = utemp_1 + utemp_2 * bp_0;
      utemp2_2 = utemp_2 + utemp_0 * bp_1;

      utemp2_0 -= utemp_2 * bp_1;
      utemp2_1 -= utemp_0 * bp_2;
      utemp2_2 -= utemp_1 * bp_0;

      const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

      bp_0 *= otsq;
      bp_1 *= otsq;
      bp_2 *= otsq;

      utemp_0 += utemp2_1 * bp_2;
      utemp_1 += utemp2_2 * bp_0;
      utemp_2 += utemp2_0 * bp_1;

      utemp_0 -= utemp2_2 * bp_1;
      utemp_1 -= utemp2_0 * bp_2;
      utemp_2 -= utemp2_1 * bp_0;

      p.fd(part_idx,0) = utemp_0 + ep_0;
      p.fd(part_idx,1) = utemp_1 + ep_1;
      p.fd(part_idx,2) = utemp_2 + ep_2;

      nloops++;

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// Advance velocities and calculate time-centered total energy
// with fields centered at cell corner
// See DUDT_BORIS_1D for implementation with comments for the parts that are similar 
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_2D_CS( t_chunk ***chunks, int *block_start_tiles,
                                  int *block_start_chunks, int *block_num_par,
                                  p_k_part rqm, t_f2<p_k_fld> *b_f2,
                                  t_f2<p_k_fld> *e_f2, p_double *energy,
                                  bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  __shared__ int til, n0, n1, n2, offset;
  __shared__ p_k_part tem;
  __shared__ t_f2<p_k_fld> b_f2_til, e_f2_til;

  extern __shared__ p_k_fld s[];
  __shared__ p_k_fld *e_s;
  __shared__ p_k_fld *b_s;

  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f2_til = b_f2[til];
    e_f2_til = e_f2[til];

    n0 = b_f2_til.n0;
    n1 = b_f2_til.n1;
    n2 = b_f2_til.n2;
    offset = b_f2_til.g - b_f2_til.f;

    // divvy up shared mem array
    b_s = (p_k_fld *)s;
    e_s = (p_k_fld *)&b_s[n0*n1*n2];

    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // populate shared mem arrays
  int i = tx;
  while ( i < n0 * n1 * n2 ) {
    b_s[i] = b_f2_til.f[i];
    e_s[i] = e_f2_til.f[i];
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1];

  p_k_fld *const w1 = w1_buff - LP;
  p_k_fld *const w2 = w2_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  
  // loop over all chunks assigned to this block
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);

    int nloops = 0; 

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }
    
      const int i1 = ix.fd(part_idx,0) - 1;
      const int i2 = ix.fd(part_idx,1) - 1;

      const p_k_part dx1 = p_k_fld( x.fd(part_idx,0) );
      const p_k_part dx2 = p_k_fld( x.fd(part_idx,1) );

      p_k_fld utemp2_0 = p.fd(part_idx,0);
      p_k_fld utemp2_1 = p.fd(part_idx,1);
      p_k_fld utemp2_2 = p.fd(part_idx,2);
      
      SPLINE( dx1, w1 );
      SPLINE( dx2, w2 );

      p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
      
      // Interpolate E - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;
        
        for (int k2=LP; k2<=UP; ++k2){
          p_k_fld f1line = 0;
          p_k_fld f2line = 0;
          p_k_fld f3line = 0;
          
          for (int k1=LP; k1<=UP; ++k1){
            // f1line += e_f2_til.fd(0,i1+k1,i2+k2) * w1[k1];
            // f2line += e_f2_til.fd(1,i1+k1,i2+k2) * w1[k1];
            // f3line += e_f2_til.fd(2,i1+k1,i2+k2) * w1[k1];
            f1line += e_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 0 ] * w1[ k1];
            f2line += e_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 1 ] * w1[ k1];
            f3line += e_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 2 ] * w1[ k1];
          }
          
          f1 += f1line * w2[k2];
          f2 += f2line * w2[k2];
          f3 += f3line * w2[k2];
        }
        
        ep_0 = f1;
        ep_1 = f2;
        ep_2 = f3;
      }
      
      // Interpolate B - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;
      
        for (int k2=LP; k2<=UP; ++k2){
          p_k_fld f1line = 0;
          p_k_fld f2line = 0;
          p_k_fld f3line = 0;

          for (int k1=LP; k1<=UP; ++k1){
            // f1line += b_f2_til.fd(0,i1+k1,i2+k2) * w1[k1];
            // f2line += b_f2_til.fd(1,i1+k1,i2+k2) * w1[k1];
            // f3line += b_f2_til.fd(2,i1+k1,i2+k2) * w1[k1];
            f1line += b_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 0 ] * w1[ k1];
            f2line += b_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 1 ] * w1[ k1];
            f3line += b_s[ offset + ( n1*(i2 +k2) + (i1 +k1) ) * n0 + 2 ] * w1[ k1];
          }
          
        f1 += f1line * w2[k2];
        f2 += f2line * w2[k2];
        f3 += f3line * w2[k2];
        }
        
        bp_0 = f1;
        bp_1 = f2;
        bp_2 = f3;
      }

      ep_0 *= tem;
      ep_1 *= tem;
      ep_2 *= tem;

      p_k_fld utemp_0 = utemp2_0 + ep_0;
      p_k_fld utemp_1 = utemp2_1 + ep_1;
      p_k_fld utemp_2 = utemp2_2 + ep_2;

      const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

      const p_k_part gamma = sqrt(u2+1);

      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      const p_k_part gam_tem = tem / gamma;

      bp_0 *= gam_tem;
      bp_1 *= gam_tem;
      bp_2 *= gam_tem;

      utemp2_0 = utemp_0 + utemp_1 * bp_2;
      utemp2_1 = utemp_1 + utemp_2 * bp_0;
      utemp2_2 = utemp_2 + utemp_0 * bp_1;

      utemp2_0 -= utemp_2 * bp_1;
      utemp2_1 -= utemp_0 * bp_2;
      utemp2_2 -= utemp_1 * bp_0;

      const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

      bp_0 *= otsq;
      bp_1 *= otsq;
      bp_2 *= otsq;

      utemp_0 += utemp2_1 * bp_2;
      utemp_1 += utemp2_2 * bp_0;
      utemp_2 += utemp2_0 * bp_1;

      utemp_0 -= utemp2_2 * bp_1;
      utemp_1 -= utemp2_0 * bp_2;
      utemp_2 -= utemp2_1 * bp_0;

      p.fd(part_idx,0) = utemp_0 + ep_0;
      p.fd(part_idx,1) = utemp_1 + ep_1;
      p.fd(part_idx,2) = utemp_2 + ep_2;

      nloops++;
      
    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// Advance velocities and calculate time-centered total energy
// with fields centered at cell corner
// See DUDT_BORIS_1D for implementation with comments for the parts that are similar 
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_3D_CS( t_chunk ***chunks, int *block_start_tiles,
                                  int *block_start_chunks, int *block_num_par,
                                  p_k_part rqm, t_f3<p_k_fld> *b_f3,
                                  t_f3<p_k_fld> *e_f3, p_double *energy,
                                  bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  __shared__ int til, n0, n1, n2, n3, offset;
  __shared__ p_k_part tem;
  __shared__ t_f3<p_k_fld> b_f3_til, e_f3_til;

  extern __shared__ p_k_fld s[];
  __shared__ p_k_fld *e_s;
  __shared__ p_k_fld *b_s;

  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f3_til = b_f3[til];
    e_f3_til = e_f3[til];

    n0 = b_f3_til.n0;
    n1 = b_f3_til.n1;
    n2 = b_f3_til.n2;
    n3 = b_f3_til.n3;
    offset = b_f3_til.g - b_f3_til.f;

    // divvy up shared mem array
    b_s = (p_k_fld *)s;
    e_s = (p_k_fld *)&b_s[n0*n1*n2*n3];

    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // populate shared mem arrays
  int i = tx;
  while ( i < n0 * n1 * n2 * n3 ) {
    b_s[i] = b_f3_til.f[i];
    e_s[i] = e_f3_til.f[i];
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int num_par_proc_chunk; //, part_idx;

  p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1], w3_buff[UP-LP+1];

  p_k_fld *const w1 = w1_buff - LP;
  p_k_fld *const w2 = w2_buff - LP;
  p_k_fld *const w3 = w3_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  
  // loop over all chunks assigned to this block
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
  
    int nloops = 0; 

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }

      const int i1 = ix.fd(part_idx,0) - 1;
      const int i2 = ix.fd(part_idx,1) - 1;
      const int i3 = ix.fd(part_idx,2) - 1;

      const p_k_fld dx1 = p_k_fld( x.fd(part_idx,0) );
      const p_k_fld dx2 = p_k_fld( x.fd(part_idx,1) );
      const p_k_fld dx3 = p_k_fld( x.fd(part_idx,2) );

      p_k_fld utemp2_0 = p.fd(part_idx,0);
      p_k_fld utemp2_1 = p.fd(part_idx,1);
      p_k_fld utemp2_2 = p.fd(part_idx,2);

      SPLINE( dx1, w1 );
      SPLINE( dx2, w2 );
      SPLINE( dx3, w3 );

      p_k_fld ep_0, ep_1, ep_2, bp_0, bp_1, bp_2;
            
      // Interpolate E - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k3=LP; k3<=UP; ++k3){
          p_k_fld f1plane = 0;
          p_k_fld f2plane = 0;
          p_k_fld f3plane = 0;

          for (int k2=LP; k2<=UP; ++k2){
            p_k_fld f1line = 0;
            p_k_fld f2line = 0;
            p_k_fld f3line = 0;

            for (int k1=LP; k1<=UP; ++k1){
              // f1line += e_f3_til.fd(0,i1+k1,i2+k2,i3+k3) * w1[k1];
              // f2line += e_f3_til.fd(1,i1+k1,i2+k2,i3+k3) * w1[k1];
              // f3line += e_f3_til.fd(2,i1+k1,i2+k2,i3+k3) * w1[k1];
              f1line += e_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 0 ] * w1[ k1];
              f2line += e_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 1 ] * w1[ k1];
              f3line += e_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 2 ] * w1[ k1];
            }

            f1plane += f1line * w2[k2];
            f2plane += f2line * w2[k2];
            f3plane += f3line * w2[k2];
          }

          f1 += f1plane * w3[k3];
          f2 += f2plane * w3[k3];
          f3 += f3plane * w3[k3];

        }

        ep_0 = f1;
        ep_1 = f2;
        ep_2 = f3;
      }
      
      // Interpolate B - Field
      {
        p_k_fld f1 = 0;
        p_k_fld f2 = 0;
        p_k_fld f3 = 0;

        for (int k3=LP; k3<=UP; ++k3){
          p_k_fld f1plane = 0;
          p_k_fld f2plane = 0;
          p_k_fld f3plane = 0;

          for (int k2=LP; k2<=UP; ++k2){
            p_k_fld f1line = 0;
            p_k_fld f2line = 0;
            p_k_fld f3line = 0;

            for (int k1=LP; k1<=UP; ++k1){
              // f1line += b_f3_til.fd(0,i1+k1,i2+k2,i3+k3) * w1[k1];
              // f2line += b_f3_til.fd(1,i1+k1,i2+k2,i3+k3) * w1[k1];
              // f3line += b_f3_til.fd(2,i1+k1,i2+k2,i3+k3) * w1[k1];
              f1line += b_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 0 ] * w1[ k1];
              f2line += b_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 1 ] * w1[ k1];
              f3line += b_s[ offset + ( ( n2*(i3 +k3) + i2 +k2 ) * n1 + i1 +k1 ) * n0 + 2 ] * w1[ k1];
            }

            f1plane += f1line * w2[k2];
            f2plane += f2line * w2[k2];
            f3plane += f3line * w2[k2];
          }

          f1 += f1plane * w3[k3];
          f2 += f2plane * w3[k3];
          f3 += f3plane * w3[k3];

        }

        bp_0 = f1;
        bp_1 = f2;
        bp_2 = f3;
      }

      ep_0 *= tem;
      ep_1 *= tem;
      ep_2 *= tem;

      p_k_fld utemp_0 = utemp2_0 + ep_0;
      p_k_fld utemp_1 = utemp2_1 + ep_1;
      p_k_fld utemp_2 = utemp2_2 + ep_2;

      const p_k_part u2 = utemp_0*utemp_0 + utemp_1*utemp_1 + utemp_2*utemp_2;

      const p_k_part gamma = sqrt(u2+1);

      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      const p_k_part gam_tem = tem / gamma;

      bp_0 *= gam_tem;
      bp_1 *= gam_tem;
      bp_2 *= gam_tem;

      utemp2_0 = utemp_0 + utemp_1 * bp_2;
      utemp2_1 = utemp_1 + utemp_2 * bp_0;
      utemp2_2 = utemp_2 + utemp_0 * bp_1;

      utemp2_0 -= utemp_2 * bp_1;
      utemp2_1 -= utemp_0 * bp_2;
      utemp2_2 -= utemp_1 * bp_0;

      const p_k_part otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp_0*bp_0) + bp_1*bp_1 ) + bp_2*bp_2 );

      bp_0 *= otsq;
      bp_1 *= otsq;
      bp_2 *= otsq;

      utemp_0 += utemp2_1 * bp_2;
      utemp_1 += utemp2_2 * bp_0;
      utemp_2 += utemp2_0 * bp_1;

      utemp_0 -= utemp2_2 * bp_1;
      utemp_1 -= utemp2_0 * bp_2;
      utemp_2 -= utemp2_1 * bp_0;

      p.fd(part_idx,0) = utemp_0 + ep_0;
      p.fd(part_idx,1) = utemp_1 + ep_1;
      p.fd(part_idx,2) = utemp_2 + ep_2;

      nloops++;

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


//----------------------------------------------------------------------------------------
// use loop. this current deposition works for any interpolation order. 
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void DEPOSIT_CURRENT_ROLLED_1D(
        const int ix, 
        const p_k_part x0, 
        const p_k_part x1, 
        const p_k_part qnx, const p_k_part qvy, const p_k_part qvz,
        p_k_fld *j,
        const int offset, const int n0 )
{
  p_k_fld S0x_buff[P1-P0+1], S1x_buff[P1-P0+1];
  p_k_fld wl1_buff[P1-P0+1];

  // Obtain offset pointers to the above buffers to enable non-zero indexing
  p_k_fld *const S0x = S0x_buff - P0;
  p_k_fld *const S1x = S1x_buff - P0;
  p_k_fld *const wl1 = wl1_buff - P0;

  // get spline weights for x and y
  SPLINE( x0, S0x );
  SPLINE( x1, S1x );

  // get longitudinal motion weights
  // the last vaiue is set to 0 so we can accumulate all current components in a single
  // pass
  WL( qnx, x0, x1, wl1 );
  wl1[P1] = 0;

  // j1 
  for (int k1=P0; k1<P1; ++k1){

    // atomicAdd( jay_f1_til.fdp(0,ii+k1), wl1[k1] );
    atomicAdd( &j[ offset + n0*(ix+k1) + 0 ], wl1[k1] );

  }

  // j2, j3 
  for (int k1=P0; k1<=P1; ++k1){

    p_k_fld pw1 = _0_5_p_k_fld * (S0x[k1]+S1x[k1]);

    // atomicAdd( jay_f1_til.fdp(1,ii+k1), qvy * pw1 );
    // atomicAdd( jay_f1_til.fdp(2,ii+k1), qvz * pw1 );
    atomicAdd( &j[ offset + n0*(ix+k1) + 1 ], qvy * pw1 );
    atomicAdd( &j[ offset + n0*(ix+k1) + 2 ], qvz * pw1 );

  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void ADVANCE_DEPOSIT_1D( t_chunk ***chunks, int *block_start_tiles,
                                    int *block_start_chunks, int *block_num_par,
                                    p_double dx1,
                                    t_f1<p_k_fld> *jay_f1, p_double dt,
                                    int chunk_size, int *n_hole, int *n_comm_ext,
                                    int *n_recv_loc ){

  // Block, thread index
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til, n0, n1, offset;
  __shared__ t_f1<p_k_fld> jay_f1_til;
  __shared__ p_k_part dt_dx1;

  extern __shared__ p_k_fld s[];

  if (tx==0) {
    til = block_start_tiles[bx];

    jay_f1_til = jay_f1[til];

    n0 = jay_f1_til.n0;
    n1 = jay_f1_til.n1;
    offset = jay_f1_til.g - jay_f1_til.f;

    dt_dx1 = p_k_part( dt/dx1 );
  }
  __syncthreads();

  // zero out shared memory array
  int i = tx;
  while ( i < n0 * n1 ) {
    s[i] = 0;
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  const p_k_fld dx1_dt = p_k_fld( dx1/dt );

  // const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);

  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    
    int nloops = 0;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }

      // read in all data 
      const p_k_part q0 = q.fd(part_idx);
      const p_k_part x0 = x.fd(part_idx,0); 
      const int ix0 = ix.fd(part_idx,0) - 1;
      const p_k_part px = p.fd(part_idx,0);
      const p_k_part py = p.fd(part_idx,1);
      const p_k_part pz = p.fd(part_idx,2);

      const p_k_part rgamma = _1_0_p_k_part / sqrt( _1_0_p_k_part + px*px + py*py + pz*pz );
      
      const p_k_part delta_x = px * rgamma * dt_dx1;

      const p_k_part x1 = x0 + delta_x;

      const int delta_ix = ntrim(x1);

      // virtual particle data
      int v_ix[2];
      p_k_part v_x0[2];
      p_k_part v_x1[2];
      p_k_part v_qvy[2];
      p_k_part v_qvz[2];
      
      p_k_part xint, eps;

      // number of virtual particles
      int nvp = 1;

      const p_k_part qnx = q0 * dx1_dt;
      const p_k_part qvy = q0 * py * rgamma;
      const p_k_part qvz = q0 * pz * rgamma;

      // initial position for virtual particle 0 is the same in all cases 
      v_ix[0] = ix0;
      v_x0[0] = x0;
      
      // Splits particle trajectories so that all virtual particles have a motion starting and ending
      // in the same cell. This routines also calculates vz for each virtual particle.

      if ( delta_ix == 0 ) {

        nvp = 1; 
        
        v_x1[0] = x1;
        v_qvy[0] = qvy;
        v_qvz[0] = qvz;

      } else {

        nvp = 2; 

        // interp 
        xint  = _0_5_p_k_fld * delta_ix;
        eps = ( xint - x0 ) / delta_x;

        // v0 
        v_x1[0] = xint;
        v_qvy[0] = qvy * eps;
        v_qvz[0] = qvz * eps;

        // v1 
        v_ix[1] = ix0 + delta_ix ;
        v_x0[1] = -xint;
        v_x1[1] = x1 - delta_ix;
        v_qvy[1] = qvy * (1-eps);
        v_qvz[1] = qvz * (1-eps);
      }

      for (int i=0; i<nvp; ++i){

        DEPOSIT_CURRENT_1D( v_ix[i], v_x0[i], v_x1[i],
                            qnx, v_qvy[i], v_qvz[i],
                            s, offset, n0 );
        
      }

      // update particle positions
      x.fd(part_idx,0) = x1 - delta_ix;
      ix.fd(part_idx,0) += delta_ix;

      nloops++; 

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks  ``while ( num_par_remainin > 0 )``

  __syncthreads();
      
  // Add shared memory current array back to global
  i = tx;
  while ( i < n0 * n1 ) {
    atomicAdd( &jay_f1_til.f[i], s[i] );
    i += BLOCK_SIZE;
  }

  // one thread on each tile zero out counters to be used in sort
  if ( tx==0 and block_start_chunks[bx]==0 ) {
    n_hole[til] = 0;
    n_comm_ext[til] = 0;
    n_recv_loc[til] = 0;
  }

}



//----------------------------------------------------------------------------------------
// use loop. this current deposition works for any interpolation order. 
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void DEPOSIT_CURRENT_ROLLED_2D(
        const int ix, const int iy,
        const p_k_part x0, const p_k_part y0,
        const p_k_part x1, const p_k_part y1, 
        const p_k_part qnx, const p_k_part qny, const p_k_part qvz,
        p_k_fld *j,
        const int offset, const int n0, const int n1 )
{

  p_k_fld S0x_buff[P1-P0+1], S1x_buff[P1-P0+1], S0y_buff[P1-P0+1], S1y_buff[P1-P0+1];
  p_k_fld wp1_buff[P1-P0+1], wp2_buff[P1-P0+1];
  p_k_fld wl1_buff[P1-P0+1], wl2_buff[P1-P0+1];

  p_k_fld tmp1, tmp2;
  
  // Obtain offset pointers to the above buffers to enable non-zero indexing
  p_k_fld *const S0x = S0x_buff - P0;
  p_k_fld *const S1x = S1x_buff - P0;
  p_k_fld *const S0y = S0y_buff - P0;
  p_k_fld *const S1y = S1y_buff - P0;
  p_k_fld *const wp1 = wp1_buff - P0;
  p_k_fld *const wp2 = wp2_buff - P0;
  p_k_fld *const wl1 = wl1_buff - P0;
  p_k_fld *const wl2 = wl2_buff - P0;

  // get spline weights for x and y
  SPLINE( x0, S0x );
  SPLINE( x1, S1x );

  SPLINE( y0, S0y );
  SPLINE( y1, S1y );

  // get longitudinal motion weights
  // the last vaiue is set to 0 so we can accumulate all current components in a single
  // pass
  WL( qnx, x0, x1, wl1 );
  wl1[P1] = 0;

  WL( qny, y0, y1, wl2 );
  wl2[P1] = 0;

  // get perpendicular motion weights
  // (a division by 2 was moved to the jnorm vars. above)
  for (int k1=P0; k1<=P1; ++k1){
    wp1[k1] = S0y[k1] + S1y[k1];
    wp2[k1] = S0x[k1] + S1x[k1];
  }

  // j1 
  for (int k2=P0; k2<=P1; ++k2){
    for (int k1=P0; k1<P1; ++k1){
      atomicAdd( &j[ offset + ( n1*(iy+k2) + (ix+k1) ) * n0 + 0 ], wl1[k1] * wp1[k2] );
    }
  }

  // j2 
  for (int k2=P0; k2<P1; ++k2){
    for (int k1=P0; k1<=P1; ++k1){
      atomicAdd( &j[ offset + ( n1*(iy+k2) + (ix+k1) ) * n0 + 1 ], wp2[k1] * wl2[k2] );
    }
  }

  // j3 
  for (int k2=P0; k2<=P1; ++k2){
    for (int k1=P0; k1<=P1; ++k1){
      tmp1 = S0x[k1]*S0y[k2] + S1x[k1]*S1y[k2];
      tmp2 = S0x[k1]*S1y[k2] + S1x[k1]*S0y[k2];
      atomicAdd( &j[ offset + ( n1*(iy+k2) + (ix+k1) ) * n0 + 2 ], qvz * ( tmp1 + _0_5_p_k_fld * tmp2 ) );
    }
  }
}



// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void ADVANCE_DEPOSIT_2D( t_chunk ***chunks, int *block_start_tiles,
                                    int *block_start_chunks, int *block_num_par,
                                    p_double dx1, p_double dx2,
                                    t_f2<p_k_fld> *jay_f2, p_double dt,
                                    int chunk_size, int *n_hole, int *n_comm_ext,
                                    int *n_recv_loc ){

  // Block, thread index
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til, n0, n1, n2, offset;
  __shared__ t_f2<p_k_fld> jay_f2_til;
  __shared__ p_k_part dt_dx1, dt_dx2;

  extern __shared__ p_k_fld s[];

  if (tx==0) {
    til = block_start_tiles[bx];

    jay_f2_til = jay_f2[til];

    n0 = jay_f2_til.n0;
    n1 = jay_f2_til.n1;
    n2 = jay_f2_til.n2;
    offset = jay_f2_til.g - jay_f2_til.f;

    dt_dx1 = p_k_part( dt/dx1 );
    dt_dx2 = p_k_part( dt/dx2 );
  }
  __syncthreads();

  // zero out shared memory array
  int i = tx;
  while ( i < n0 * n1 * n2 ) {
    s[i] = 0;
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);
  const p_k_fld jnorm1 = p_k_fld( dx1/dt/2 );
  const p_k_fld jnorm2 = p_k_fld( dx2/dt/2 );

  
  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    __syncwarp();
    
    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    
    int nloops = 0;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }

      // read in all data 
      const p_k_part q0 = q.fd(part_idx);
      const p_k_part x0 = x.fd(part_idx,0); 
      const p_k_part y0 = x.fd(part_idx,1); 
      const int ix0 = ix.fd(part_idx,0) - 1;
      const int iy0 = ix.fd(part_idx,1) - 1;
      const p_k_part px = p.fd(part_idx,0);
      const p_k_part py = p.fd(part_idx,1);
      const p_k_part pz = p.fd(part_idx,2);
        
      const p_k_part rgamma = _1_0_p_k_part / sqrt( _1_0_p_k_part + px*px + py*py + pz*pz );
      
      const p_k_part delta_x = px * rgamma * dt_dx1;
      const p_k_part delta_y = py * rgamma * dt_dx2;
      
      const p_k_part x1 = x0 + delta_x;
      const p_k_part y1 = y0 + delta_y; 

      const int delta_ix = ntrim(x1);
      const int delta_iy = ntrim(y1);

      int v_ix[3];
      int v_iy[3];
      p_k_part v_x0[3];
      p_k_part v_y0[3];
      p_k_part v_x1[3];
      p_k_part v_y1[3];
      p_k_part v_qvz[3];
      
      p_k_part xint, yint, xint2, yint2, eps, eps2;

      // number of virtual particles
      int nvp = 1;

      const p_k_part qnx = q0 * jnorm1;
      const p_k_part qny = q0 * jnorm2;
      const p_k_part qvz = q0 * pz * rgamma * c1_3;

      // initial position for virtual particle 0 is the same in all cases 
      v_ix[0] = ix0;
      v_iy[0] = iy0;
      v_x0[0] = x0;
      v_y0[0] = y0;

      // Split particle trajectories so that all virtual particles have a motion starting and ending
      // in the same cell. This routines also calculates vz for each virtual particle.

      switch ( abs(delta_ix) + 2 * abs(delta_iy) ) {
        // no cross
        case 0:

          // v0 
          v_x1[0] = x1;
          v_y1[0] = y1;
          v_qvz[0] = qvz; 

          // printf("--CUDA DEBUG: til %d, part %d k %d x0 %f, y0 %f, x1 %f, y1 %f, qq %f, vz %f, ii %d, jj %d \n",
          //     bx+1, part_idx+1, k, x0[k], y0[k], x1[k], y1[k], qq[k], vz[k], ii[k], jj[k]);

          break;

        // x cross only
        case 1:

          nvp = 2;

          // interp 
          xint  = _0_5_p_k_fld * delta_ix;
          eps = ( xint - x0 ) / delta_x;
          yint  =  y0 + delta_y * eps;

          // v0 
          v_x1[0] = xint;
          v_y1[0] = yint;
          v_qvz[0] = qvz * eps;

          // v1 
          v_ix[1] = ix0 + delta_ix;
          v_iy[1] = iy0;
          v_x0[1] = -xint;
          v_y0[1] = yint;
          v_x1[1] = x1 - delta_ix;
          v_y1[1] = y1;
          v_qvz[1] = qvz * (1-eps);

          break;

        // y cross only
        case 2:

          nvp = 2; 

          // interp 
          yint = _0_5_p_k_fld * delta_iy;
          eps = ( yint - y0 ) / delta_y;
          xint =  x0 + delta_x * eps;

          // v0
          v_x1[0] = xint;
          v_y1[0] = yint;
          v_qvz[0] = qvz * eps;

          // v1
          v_ix[1] = ix0;
          v_iy[1] = iy0 + delta_iy;
          v_x0[1] = xint;
          v_y0[1] = -yint;
          v_x1[1] = x1;
          v_y1[1] = y1 - delta_iy;
          v_qvz[1] = qvz * (1-eps);

          break;

        // x, y cross
        case 3:

          nvp = 3;
          
          // split in x direction first
          xint = _0_5_p_k_fld * delta_ix;
          eps = ( xint - x0 ) / delta_x;
          yint =  y0 + delta_y * eps;

          // check if y intersection occured for 1st or 2nd split
          if ( yint >= -_0_5_p_k_fld && yint < _0_5_p_k_fld ){

            // v0
            v_x1[0] = xint;
            v_y1[0] = yint;
            v_qvz[0] = qvz * eps;

            // y split 2nd vp
            yint2 = _0_5_p_k_fld * delta_iy;
            eps2 = ( yint2 - yint ) / ( y1 - yint );
            xint2 = -xint + ( x1 - xint ) * eps2;

            // v1 
            v_ix[1] = ix0 + delta_ix;
            v_iy[1] = iy0;
            v_x0[1] = -xint;
            v_y0[1] =  yint;
            v_x1[1] = xint2;
            v_y1[1] = yint2;
            v_qvz[1] = qvz * (1-eps) * eps2;

            // v2 
            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0 + delta_iy;
            v_x0[2] = xint2;
            v_y0[2] = -yint2;
            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1 - delta_iy;
            v_qvz[2] = qvz * (1-eps) * (1-eps2);


          } else {

            // y split 1st vp
            yint2 = _0_5_p_k_fld * delta_iy;
            eps2 = ( yint2 - y0 ) / ( yint - y0 );
            xint2 = x0 + (xint - x0) * eps2;

            // v0 
            v_x1[0] = xint2;
            v_y1[0] = yint2;
            v_qvz[0] = qvz * eps * eps2; 

            // v1
            v_ix[1] = ix0;
            v_iy[1] = iy0 + delta_iy;
            v_x0[1] = xint2;
            v_y0[1] = -yint2;
            v_x1[1] = xint;
            v_y1[1] = yint - delta_iy;
            v_qvz[1] = qvz * eps * (1-eps2);

            // v2 (no y cross on second vp)
            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0 + delta_iy;
            v_x0[2] = - xint;
            v_y0[2] = yint - delta_iy;
            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1 - delta_iy;
            v_qvz[2] = qvz * (1 - eps); 
          }

          break;

      }


      // // One might think that separating the v_ix array  into its separate components
      // // as in the commented section below would give a speedup from not having to do 
      // // array accesses, and similarly for the other arrays. However, we found that this
      // // gave no speedup on NVIDIA A100 GPUs. And it makes the code ridiculously way more
      // // complex particularly in the split section above in the 3D case.

      // DEPOSIT_CURRENT_2D( v0_ix, v0_iy, v0_x0, v0_y0, v0_x1, v0_y1,
      //                     qnx, qny, v0_qvz, s, offset, n0, n1 );

      // if( nvp > 1 )
      //   DEPOSIT_CURRENT_2D( v1_ix, v1_iy, v1_x0, v1_y0, v1_x1, v1_y1,
      //                       qnx, qny, v1_qvz, s, offset, n0, n1 );
      
      // if( nvp > 2 ) 
      //   DEPOSIT_CURRENT_2D( v2_ix, v2_iy, v2_x0, v2_y0, v2_x1, v2_y1,
      //                       qnx, qny, v2_qvz, s, offset, n0, n1 );
      

      for( int i=0; i<nvp; i++ ) {

        DEPOSIT_CURRENT_2D( v_ix[i], v_iy[i], v_x0[i], v_y0[i], v_x1[i], v_y1[i],
                            qnx, qny, v_qvz[i], s, offset, n0, n1 );

      }
      
      // update particle positions 
      x.fd(part_idx,0) = x1 - delta_ix;
      x.fd(part_idx,1) = y1 - delta_iy;
      ix.fd(part_idx,0) += delta_ix;
      ix.fd(part_idx,1) += delta_iy;


      nloops++; 


    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks  ``while ( num_par_remainin > 0 )``

  __syncthreads();
      
  // Add shared memory current array back to global
  i = tx;
  while ( i < n0 * n1 * n2 ) {
    atomicAdd( &jay_f2_til.f[i], s[i] );
    i += BLOCK_SIZE;
  }

  // one thread on each tile zero out counters to be used in sort
  if ( tx==0 and block_start_chunks[bx]==0 ) {
    n_hole[til] = 0;
    n_comm_ext[til] = 0;
    n_recv_loc[til] = 0;
  }

}



//----------------------------------------------------------------------------------------
// use loop. this current deposition works for any interpolation order. 
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void DEPOSIT_CURRENT_ROLLED_3D(
        const int ix, const int iy, const int iz,
        const p_k_part x0, const p_k_part y0, const p_k_part z0,
        const p_k_part x1, const p_k_part y1, const p_k_part z1, 
        const p_k_part qnx, const p_k_part qny, const p_k_part qnz,
        p_k_fld *j,
        const int offset, const int n0, const int n1, const int n2 )
{
#define NB (P1-P0+1)

  p_k_fld S0x_buff[NB], S1x_buff[NB], S0y_buff[NB], S1y_buff[NB], S0z_buff[NB], S1z_buff[NB];
  p_k_fld wp1_buff[NB*NB], wp2_buff[NB*NB], wp3_buff[NB*NB];
  p_k_fld wl1_buff[NB], wl2_buff[NB], wl3_buff[NB];

  // Obtain offset pointers to the above buffers to enable non-zero indexing
  p_k_fld *const S0x = S0x_buff - P0;
  p_k_fld *const S1x = S1x_buff - P0;
  p_k_fld *const S0y = S0y_buff - P0;
  p_k_fld *const S1y = S1y_buff - P0;
  p_k_fld *const S0z = S0z_buff - P0;
  p_k_fld *const S1z = S1z_buff - P0;
  p_k_fld *const wp1 = wp1_buff - P0 - P0*NB;
  p_k_fld *const wp2 = wp2_buff - P0 - P0*NB;
  p_k_fld *const wp3 = wp3_buff - P0 - P0*NB;
  p_k_fld *const wl1 = wl1_buff - P0;
  p_k_fld *const wl2 = wl2_buff - P0;
  p_k_fld *const wl3 = wl3_buff - P0;
  
  // get spline weights for x, y and z
  SPLINE( x0, S0x );
  SPLINE( x1, S1x );

  SPLINE( y0, S0y );
  SPLINE( y1, S1y );

  SPLINE( z0, S0z );
  SPLINE( z1, S1z );

  // get longitudinal motion weights
  // the last value is set to 0 so we can accumulate
  // all current components in a single pass
  WL( qnx, x0, x1, wl1 );
  wl1[P1] = 0;

  WL( qny, y0, y1, wl2 );
  wl2[P1] = 0;

  WL( qnz, z0, z1, wl3 );
  wl3[P1] = 0;

  // get perpendicular motion weights
  // (a division by 2 was moved to the jnorm vars. above)
  for (int k2=P0; k2<=P1; ++k2){
    for (int k1=P0; k1<=P1; ++k1){
      wp1[k1+k2*NB] = (S0y[k1]*S0z[k2] + S1y[k1]*S1z[k2]) +
        _0_5_p_k_fld * (S0y[k1]*S1z[k2] + S1y[k1]*S0z[k2]);
      wp2[k1+k2*NB] = (S0x[k1]*S0z[k2] + S1x[k1]*S1z[k2]) +
        _0_5_p_k_fld * (S0x[k1]*S1z[k2] + S1x[k1]*S0z[k2]);
      wp3[k1+k2*NB] = (S0x[k1]*S0y[k2] + S1x[k1]*S1y[k2]) +
        _0_5_p_k_fld * (S0x[k1]*S1y[k2] + S1x[k1]*S0y[k2]);
    }
  }

  // j1 
  for (int k3=P0; k3<=P1; ++k3){
    for (int k2=P0; k2<=P1; ++k2){
      for (int k1=P0; k1<P1; ++k1){

        // atomicAdd( jay_f3_til.fdp(0,v_ix[i]+k1,v_iy[i]+k2,v_iz[i]+k3), wl1[k1] * wp1[k2+k3*NB] );
        atomicAdd( &j[ offset + ( ( n2 * (iz+k3) + iy+k2 ) * n1 + ix+k1 ) * n0     ], wl1[k1] * wp1[k2+k3*NB] );

      }
    }
  }

  // j2 
  for (int k3=P0; k3<=P1; ++k3){
    for (int k2=P0; k2<P1; ++k2){
      for (int k1=P0; k1<=P1; ++k1){

        // atomicAdd( jay_f3_til.fdp(1,v_ix[i]+k1,v_iy[i]+k2,v_iz[i]+k3), wl2[k2] * wp2[k1+k3*NB] );
        atomicAdd( &j[ offset + ( ( n2 * (iz+k3) + iy+k2 ) * n1 + ix+k1 ) * n0 + 1 ], wl2[k2] * wp2[k1+k3*NB] );

      }
    }
  }

  // j3 
  for (int k3=P0; k3<P1; ++k3){
    for (int k2=P0; k2<=P1; ++k2){
      for (int k1=P0; k1<=P1; ++k1){

        // atomicAdd( jay_f3_til.fdp(2,v_ix[i]+k1,v_iy[i]+k2,v_iz[i]+k3), wl3[k3] * wp3[k1+k2*NB] );
        atomicAdd( &j[ offset + ( ( n2 * (iz+k3) + iy+k2 ) * n1 + ix+k1 ) * n0 + 2 ], wl3[k3] * wp3[k1+k2*NB] );

      }
    }
  }

#undef NB
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void ADVANCE_DEPOSIT_3D( t_chunk ***chunks, int *block_start_tiles,
                                    int *block_start_chunks, int *block_num_par,
                                    p_double dx1, p_double dx2, p_double dx3,
                                    t_f3<p_k_fld> *jay_f3, p_double dt,
                                    int chunk_size, int *n_hole, int *n_comm_ext,
                                    int *n_recv_loc ){

  // Block, thread index
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til, n0, n1, n2, n3, offset;
  __shared__ t_f3<p_k_fld> jay_f3_til;
  __shared__ p_k_part dt_dx1, dt_dx2, dt_dx3;

  extern __shared__ p_k_fld s[];

  if (tx==0) {
    til = block_start_tiles[bx];

    jay_f3_til = jay_f3[til];

    n0 = jay_f3_til.n0;
    n1 = jay_f3_til.n1;
    n2 = jay_f3_til.n2;
    n3 = jay_f3_til.n3;
    offset = jay_f3_til.g - jay_f3_til.f;

    dt_dx1 = p_k_part( dt/dx1 );
    dt_dx2 = p_k_part( dt/dx2 );
    dt_dx3 = p_k_part( dt/dx3 );
  }
  __syncthreads();

  // zero out shared memory array
  int i = tx;
  while ( i < n0 * n1 * n2 * n3 ) {
    s[i] = 0;
    i += BLOCK_SIZE;
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  const p_k_fld jnorm1 = p_k_fld( dx1/dt/3 );
  const p_k_fld jnorm2 = p_k_fld( dx2/dt/3 );
  const p_k_fld jnorm3 = p_k_fld( dx3/dt/3 );

  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    __syncwarp();
    
    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    const int num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    
    int nloops = 0;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( 1 ) {
      
      const int part_idx = tx + nloops * BLOCK_SIZE;

      __syncwarp();

      if( nloops * BLOCK_SIZE >= num_par_proc_chunk ) break; 

      if( part_idx >= num_par_proc_chunk ) {
        nloops++;
        continue;
      }

      // read in all data 
      const p_k_part q0 = q.fd(part_idx);
      const p_k_part x0 = x.fd(part_idx,0); 
      const p_k_part y0 = x.fd(part_idx,1); 
      const p_k_part z0 = x.fd(part_idx,2); 
      const int ix0 = ix.fd(part_idx,0) - 1;
      const int iy0 = ix.fd(part_idx,1) - 1;
      const int iz0 = ix.fd(part_idx,2) - 1;
      const p_k_part px = p.fd(part_idx,0);
      const p_k_part py = p.fd(part_idx,1);
      const p_k_part pz = p.fd(part_idx,2);
        
      const p_k_part rgamma = _1_0_p_k_part / sqrt( _1_0_p_k_part + px*px + py*py + pz*pz );

      const p_k_part delta_x = px * rgamma * dt_dx1;
      const p_k_part delta_y = py * rgamma * dt_dx2;
      const p_k_part delta_z = pz * rgamma * dt_dx3;
      
      const p_k_part x1 = x0 + delta_x;
      const p_k_part y1 = y0 + delta_y; 
      const p_k_part z1 = z0 + delta_z; 

      const int delta_ix = ntrim(x1);
      const int delta_iy = ntrim(y1);
      const int delta_iz = ntrim(z1);

      // virtual particle buffers 
      int v_ix[4], v_iy[4], v_iz[4];
      p_k_part v_x0[4], v_y0[4], v_z0[4], v_x1[4], v_y1[4], v_z1[4];

      p_k_part xint, yint, zint, xint2, yint2, zint2, eps;

      const p_k_part qnx = q0 * jnorm1;
      const p_k_part qny = q0 * jnorm2;
      const p_k_part qnz = q0 * jnorm3;

      // initial position for virtual particle 0 is the same in all cases 
      v_ix[0] = ix0;
      v_iy[0] = iy0;
      v_iz[0] = iz0;
      v_x0[0] = x0;
      v_y0[0] = y0;
      v_z0[0] = z0;
      
      // number of virtual particles
      int nvp;

      // Split particle trajectories so that all virtual particles have a motion starting and ending
      // in the same cell. 
      
      switch ( abs(delta_ix) + 2 * abs(delta_iy) + 4 * abs(delta_iz) ) {
        // no cross
        case 0:

          nvp = 1;
          
          v_x1[0] = x1;
          v_y1[0] = y1;
          v_z1[0] = z1;

          break;

        // x cross only
        case 1:

          nvp = 2;
          
          xint  = _0_5_p_k_fld * delta_ix;
          eps = ( xint - x0 ) / ( x1 - x0 );
          yint  =  y0 + (y1 - y0) * eps;
          zint  =  z0 + (z1 - z0) * eps;

          v_x1[0] = xint;
          v_y1[0] = yint;
          v_z1[0] = zint;

          // v1 
          v_x0[1] = -xint;
          v_y0[1] = yint;
          v_z0[1] = zint;
          v_x1[1] = x1 - delta_ix;
          v_y1[1] = y1;
          v_z1[1] = z1;

          v_ix[1] = ix0 + delta_ix;
          v_iy[1] = iy0;
          v_iz[1] = iz0;

          break;

        // y cross only
        case 2:

          nvp = 2;
          
          yint = _0_5_p_k_fld * delta_iy;
          eps = ( yint - y0 ) / ( y1 - y0);
          xint =  x0 + (x1 - x0) * eps;
          zint =  z0 + (z1 - z0) * eps;

          // v0 
          v_x1[0] = xint;
          v_y1[0] = yint;
          v_z1[0] = zint;

          // v1 
          v_x0[1] = xint;
          v_y0[1] = -yint;
          v_z0[1] = zint;
          v_x1[1] = x1;
          v_y1[1] = y1 - delta_iy;
          v_z1[1] = z1;

          v_ix[1] = ix0;
          v_iy[1] = iy0 + delta_iy;
          v_iz[1] = iz0;

          break;

        // x, y cross
        case 3:

          nvp = 3;
          
          xint = _0_5_p_k_fld * delta_ix;
          eps = ( xint - x0 ) / ( x1 - x0);
          yint =  y0 + ( y1 - y0) * eps;
          zint =  z0 + ( z1 - z0) * eps;

          if ( yint >= -_0_5_p_k_fld && yint < _0_5_p_k_fld ){

            // yspi = 1

            // no y cross on 1st vp

            // v0 
            v_x1[0] = xint;
            v_y1[0] = yint;
            v_z1[0] = zint;

            // y split 2nd vp
            yint2 = _0_5_p_k_fld * delta_iy;

            eps = ( yint2 - yint ) / (y1 - yint);
            xint2 =  -xint + (x1 - xint) * eps;
            zint2 =   zint + (z1 - zint) * eps;

            // v1 
            v_x0[1] = -xint;
            v_y0[1] =  yint;
            v_z0[1] =  zint;

            v_x1[1] = xint2;
            v_y1[1] = yint2;
            v_z1[1] = zint2;

            v_ix[1] = ix0 + delta_ix;
            v_iy[1] = iy0;
            v_iz[1] = iz0;

            // v2 
            v_x0[2] = xint2;
            v_y0[2] = -yint2;
            v_z0[2] = zint2;

            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1 - delta_iy;
            v_z1[2] = z1;

            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0 + delta_iy;
            v_iz[2] = iz0;

          } else {

            // yspi = 2

            // y split 1st vp
            yint2 = _0_5_p_k_fld * delta_iy;
            eps = ( yint2 - y0 ) / ( yint - y0);
            xint2 = x0 + (xint - x0) * eps;
            zint2 = z0 + (zint - z0) * eps;

            // v0
            v_x1[0] = xint2;
            v_y1[0] = yint2;
            v_z1[0] = zint2;

            // v1
            v_x0[1] = xint2;
            v_y0[1] = -yint2;
            v_z0[1] = zint2;

            v_x1[1] = xint;
            v_y1[1] = yint - delta_iy;
            v_z1[1] = zint;

            v_ix[1] = ix0;
            v_iy[1] = iy0 + delta_iy;
            v_iz[1] = iz0;

            // v2
            v_x0[2] = -xint;
            v_y0[2] = yint - delta_iy;
            v_z0[2] = zint;

            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1 - delta_iy;
            v_z1[2] = z1;

            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0 + delta_iy;
            v_iz[2] = iz0;

          }

          break;

        // z cross only
        case 4:

          nvp = 2;
          
          zint = _0_5_p_k_fld * delta_iz;
          eps = ( zint - z0 ) / ( z1 - z0);
          xint =  x0 + (x1 - x0) * eps;
          yint =  y0 + (y1 - y0) * eps;

          // v0 
          v_x1[0] = xint;
          v_y1[0] = yint;
          v_z1[0] = zint;

          // v1 
          v_x0[1] = xint;
          v_y0[1] = yint;
          v_z0[1] = -zint;

          v_x1[1] = x1;
          v_y1[1] = y1;
          v_z1[1] = z1 - delta_iz;

          v_ix[1] = ix0;
          v_iy[1] = iy0;
          v_iz[1] = iz0 + delta_iz;

          break;

        // x, z split - 3 vp
        case 5:

          nvp = 3; 
          
          xint = _0_5_p_k_fld * delta_ix;
          eps = ( xint - x0 ) / ( x1 - x0);
          yint =  y0 + ( y1 - y0) * eps;
          zint =  z0 + ( z1 - z0) * eps;

          zint2 = _0_5_p_k_fld * delta_iz;

          if ( zint >= -_0_5_p_k_fld && zint < _0_5_p_k_fld ){

            // zspi = 1

            // no z cross on 1st vp

            // v0 
            v_x1[0] = xint;
            v_y1[0] = yint;
            v_z1[0] = zint;

            // z split 2nd vp
            eps = ( zint2 - zint ) / (z1 - zint);
            xint2 =  -xint + (x1 - xint) * eps;
            yint2 =   yint + (y1 - yint) * eps;

            // v1 
            v_x0[1] = -xint;
            v_y0[1] =  yint;
            v_z0[1] =  zint;

            v_x1[1] = xint2;
            v_y1[1] = yint2;
            v_z1[1] = zint2;

            v_ix[1] = ix0 + delta_ix;
            v_iy[1] = iy0;
            v_iz[1] = iz0;

            // v2 
            v_x0[2] = xint2;
            v_y0[2] = yint2;
            v_z0[2] = -zint2;

            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1;
            v_z1[2] = z1 - delta_iz;

            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0;
            v_iz[2] = iz0 + delta_iz;

          } else {

            // zspi = 2

            // z split 1st vp
            eps = ( zint2 - z0 ) / ( zint - z0);
            xint2 = x0 + (xint - x0) * eps;
            yint2 = y0 + (yint - y0) * eps;

            // v0 
            v_x1[0] = xint2;
            v_y1[0] = yint2;
            v_z1[0] = zint2;

            // v1 
            v_x0[1] = xint2;
            v_y0[1] = yint2;
            v_z0[1] = -zint2;

            v_x1[1] = xint;
            v_y1[1] = yint;
            v_z1[1] = zint - delta_iz;

            v_ix[1] = ix0;
            v_iy[1] = iy0;
            v_iz[1] = iz0 + delta_iz;

            // v2 
            v_x0[2] = -xint;
            v_y0[2] = yint;
            v_z0[2] = zint - delta_iz;

            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1;
            v_z1[2] = z1 - delta_iz;

            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0;
            v_iz[2] = iz0 + delta_iz;

          }

          break;

        // y, z split - 3 vp
        case 6:

          nvp = 3;
          
          yint = _0_5_p_k_fld * delta_iy;
          eps = ( yint - y0 ) / ( y1 - y0);
          xint =  x0 + ( x1 - x0) * eps; // no x cross
          zint =  z0 + ( z1 - z0) * eps;

          if ( zint >= -_0_5_p_k_fld && zint < _0_5_p_k_fld ){

            // zspi = 1

            // no z cross on 1st vp

            // v0 
            v_x1[0] = xint;
            v_y1[0] = yint;
            v_z1[0] = zint;

            // z split 2nd vp
            zint2 = _0_5_p_k_fld * delta_iz;
            eps = ( zint2 - zint ) / (z1 - zint);
            xint2 =  xint + (x1 - xint) * eps;
            yint2 =  -yint + (y1 - yint) * eps;

            // v1 
            v_x0[1] = xint;
            v_y0[1] = -yint;
            v_z0[1] = zint;

            v_x1[1] = xint2;
            v_y1[1] = yint2;
            v_z1[1] = zint2;

            v_ix[1] = ix0;
            v_iy[1] = iy0 + delta_iy;
            v_iz[1] = iz0;

            // v2 
            v_x0[2] = xint2;
            v_y0[2] = yint2;
            v_z0[2] = -zint2;

            v_x1[2] = x1;
            v_y1[2] = y1 - delta_iy;
            v_z1[2] = z1 - delta_iz;

            v_ix[2] = ix0;
            v_iy[2] = iy0 + delta_iy;
            v_iz[2] = iz0 + delta_iz;

          } else {

            // zspi = 2

            // z split 1st vp
            zint2 = _0_5_p_k_fld * delta_iz;
            eps = ( zint2 - z0 ) / ( zint - z0);
            xint2 = x0 + (xint - x0) * eps;
            yint2 = y0 + (yint - y0) * eps;

            // v0 
            v_x1[0] = xint2;
            v_y1[0] = yint2;
            v_z1[0] = zint2;

            // v1 
            v_x0[1] = xint2;
            v_y0[1] = yint2;
            v_z0[1] = -zint2;

            v_x1[1] = xint;
            v_y1[1] = yint;
            v_z1[1] = zint - delta_iz;

            v_ix[1] = ix0;
            v_iy[1] = iy0;
            v_iz[1] = iz0 + delta_iz;

            // v2 
            v_x0[2] = xint;
            v_y0[2] = -yint;
            v_z0[2] = zint - delta_iz;

            v_x1[2] = x1;
            v_y1[2] = y1 - delta_iy;
            v_z1[2] = z1 - delta_iz;

            v_ix[2] = ix0;
            v_iy[2] = iy0 + delta_iy;
            v_iz[2] = iz0 + delta_iz;

          }

          break;

        // x, y, z cross
        case 7:

          nvp = 4;
          
          xint = _0_5_p_k_fld * delta_ix;
          eps = ( xint - x0 ) / ( x1 - x0);
          yint =  y0 + ( y1 - y0) * eps;
          zint =  z0 + ( z1 - z0) * eps;

          if ( yint >= -_0_5_p_k_fld && yint < _0_5_p_k_fld ){

            // no y cross on 1st vp

            // v0 
            v_x1[0] = xint;
            v_y1[0] = yint;
            v_z1[0] = zint;

            // y split 2nd vp
            yint2 = _0_5_p_k_fld * delta_iy;

            eps = ( yint2 - yint ) / (y1 - yint);
            xint2 =  -xint + (x1 - xint) * eps;
            zint2 =   zint + (z1 - zint) * eps;

            // v1 
            v_x0[1] = -xint;
            v_y0[1] =  yint;
            v_z0[1] =  zint;

            v_x1[1] = xint2;
            v_y1[1] = yint2;
            v_z1[1] = zint2;

            v_ix[1] = ix0 + delta_ix;
            v_iy[1] = iy0;
            v_iz[1] = iz0;

            // v2 
            v_x0[2] = xint2;
            v_y0[2] = -yint2;
            v_z0[2] = zint2;

            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1 - delta_iy;
            v_z1[2] = z1;

            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0 + delta_iy;
            v_iz[2] = iz0;

          } else {

            // y split 1st vp
            yint2 = _0_5_p_k_fld * delta_iy;
            eps = ( yint2 - y0 ) / ( yint - y0);
            xint2 = x0 + (xint - x0) * eps;
            zint2 = z0 + (zint - z0) * eps;

            // v0 
            v_x1[0] = xint2;
            v_y1[0] = yint2;
            v_z1[0] = zint2;

            // v1
            v_x0[1] = xint2;
            v_y0[1] = -yint2;
            v_z0[1] = zint2;

            v_x1[1] = xint;
            v_y1[1] = yint - delta_iy;
            v_z1[1] = zint;

            v_ix[1] = ix0;
            v_iy[1] = iy0 + delta_iy;
            v_iz[1] = iz0;

            // v2 
            v_x0[2] = -xint;
            v_y0[2] = yint - delta_iy;
            v_z0[2] = zint;

            v_x1[2] = x1 - delta_ix;
            v_y1[2] = y1 - delta_iy;
            v_z1[2] = z1;

            v_ix[2] = ix0 + delta_ix;
            v_iy[2] = iy0 + delta_iy;
            v_iz[2] = iz0;

          }

          // one of the 3 vp requires an additional z split
          zint = _0_5_p_k_fld * delta_iz;

          for ( int vpidx=0; vpidx<3; vpidx++ ) {
            if ( v_z1[vpidx] < -_0_5_p_k_fld or v_z1[vpidx] >= _0_5_p_k_fld ){
              eps = ( zint - v_z0[vpidx] ) / ( v_z1[vpidx] - v_z0[vpidx] );
              xint =  v_x0[vpidx] + (v_x1[vpidx] - v_x0[vpidx]) * eps;
              yint =  v_y0[vpidx] + (v_y1[vpidx] - v_y0[vpidx]) * eps;

              v_x0[3] = xint;
              v_y0[3] = yint;
              v_z0[3] = -zint;

              v_x1[3] = v_x1[vpidx];
              v_y1[3] = v_y1[vpidx];
              v_z1[3] = v_z1[vpidx] - delta_iz;

              v_ix[3] = v_ix[vpidx];
              v_iy[3] = v_iy[vpidx];
              v_iz[3] = v_iz[vpidx] + delta_iz;

              // correct old vp
              v_x1[vpidx] = xint;
              v_y1[vpidx] = yint;
              v_z1[vpidx] = zint;

              // correct remaining vp
              for ( int v_iyj=vpidx+1; v_iyj<3; v_iyj++ ) {
                v_z0[v_iyj] -= delta_iz;
                v_z1[v_iyj] -= delta_iz;
                v_iz[v_iyj] += delta_iz;
              }

              break;
            }
          }

          break;

      }
      
      // --------------

      // // One might think that separating the v_ix array  into its separate components
      // // as in the commented section below would give a speedup from not having to do 
      // // array accesses, and similarly for the other arrays. However, we found that this
      // // gave no speedup on NVIDIA A100 GPUs. And it makes the code ridiculously way more
      // // complex particularly in the split section above in the 3D case.

      // now accumulate jay iooping through all virtual particles
      // and shifting grid indexes

      // DEPOSIT_CURRENT_3D( v0_ix, v0_iy, v0_x0, v0_y0, v0_x1, v0_y1,
      //                     qnx, qny, v0_qvz, s, offset, n0, n1 );

      // if( nvp > 1 )
      //   DEPOSIT_CURRENT_3D( v1_ix, v1_iy, v1_x0, v1_y0, v1_x1, v1_y1,
      //                       qnx, qny, v1_qvz, s, offset, n0, n1 );
      
      // if( nvp > 2 ) 
      //   DEPOSIT_CURRENT_3D( v2_ix, v2_iy, v2_x0, v2_y0, v2_x1, v2_y1,
      //                       qnx, qny, v2_qvz, s, offset, n0, n1 );

      // if( nvp > 3 ) 
      //   DEPOSIT_CURRENT_3D( v2_ix, v2_iy, v2_x0, v2_y0, v2_x1, v2_y1,
      //                       qnx, qny, v2_qvz, s, offset, n0, n1 );      

      
      for (int i=0; i<nvp; i++){

        DEPOSIT_CURRENT_3D( v_ix[i], v_iy[i], v_iz[i],
                            v_x0[i], v_y0[i], v_z0[i],
                            v_x1[i], v_y1[i], v_z1[i],
                            qnx, qny, qnz,
                            s, offset, n0, n1, n2 ) ;
        
      }

      x.fd(part_idx,0) = x1 - delta_ix;
      x.fd(part_idx,1) = y1 - delta_iy;
      x.fd(part_idx,2) = z1 - delta_iz;
      ix.fd(part_idx,0) += delta_ix;
      ix.fd(part_idx,1) += delta_iy;
      ix.fd(part_idx,2) += delta_iz;

      nloops++;

    } // end loop over chunk of particles ``while ( 1 )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks  ``while ( num_par_remainin > 0 )``

  __syncthreads();
      
  // Write shared memory current array back to global
  i = tx;
  while ( i < n0 * n1 * n2 * n3 ) {
    atomicAdd( &jay_f3_til.f[i], s[i] );
    i += BLOCK_SIZE;
  }

  // one thread on each tile zero out counters to be used in sort
  if ( tx==0 and block_start_chunks[bx]==0 ) {
    n_hole[til] = 0;
    n_comm_ext[til] = 0;
    n_recv_loc[til] = 0;
  }

}




#undef DUDT_BORIS_1D
#undef DUDT_BORIS_2D
#undef DUDT_BORIS_3D
#undef DUDT_BORIS_1D_CS
#undef DUDT_BORIS_2D_CS
#undef DUDT_BORIS_3D_CS
#undef SPLINE
#undef SPLINEH
#undef LP
#undef UP
#undef WL
#undef P0
#undef P1
#undef ADVANCE_DEPOSIT_1D
#undef ADVANCE_DEPOSIT_2D
#undef ADVANCE_DEPOSIT_3D
#undef DEPOSIT_CURRENT_ROLLED_1D
#undef DEPOSIT_CURRENT_ROLLED_2D
#undef DEPOSIT_CURRENT_ROLLED_3D
#undef DEPOSIT_CURRENT_1D
#undef DEPOSIT_CURRENT_2D
#undef DEPOSIT_CURRENT_3D

#endif

#ifndef __TEMPLATE__

// ---------------------------------------------------------------------------------------
// Function stub to call kernel
// ---------------------------------------------------------------------------------------
void dudt_boris_wrapper( t_chunk ***chunks, int *block_start_tiles,
                         int *block_start_chunks, int *block_num_par,
                         p_k_part rqm, t_emf_c *emf, int interpolation, bool grid_center,
                         p_double *energy, bool report_energy, p_double dt,
                         int nblocks, int chunk_size, cudaStream_t stream,
                         size_t shmem_size ){

  if (grid_center) {
    // Interpolate fields centered at the cell corner
    switch( interpolation ) {
      case p_linear:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_cs1<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_cs1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_cs1<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_cs1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_cs1<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_cs1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quadratic:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_cs2<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_cs2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_cs2<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_cs2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_cs2<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_cs2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_cubic:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_cs3<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_cs3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_cs3<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_cs3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_cs3<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_cs3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quartic:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_cs4<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_cs4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_cs4<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_cs4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_cs4<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_cs4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
    }
  } else {
    // Interpolate fields on a staggered Yee mesh
    switch( interpolation ) {
      case p_linear:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_s1<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_s1<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_s1<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quadratic:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_s2<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_s2<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute(dudt_boris_3d_s2<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, shmem_size);

            dudt_boris_3d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_cubic:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_s3<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_s3<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_s3<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quartic:
        switch( p_x_dim ) {
          case 1:
            cudaFuncSetAttribute( dudt_boris_1d_s4<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_1d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            cudaFuncSetAttribute( dudt_boris_2d_s4<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_2d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            cudaFuncSetAttribute( dudt_boris_3d_s4<__BLOCK_SIZE__>, 
                cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

            dudt_boris_3d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
    }
  }

  chkLastErr("dudt_boris");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}


// ---------------------------------------------------------------------------------------
// Function stub to call kernel
// ---------------------------------------------------------------------------------------
void advance_deposit_wrapper( t_chunk ***chunks, int *block_start_tiles,
                              int *block_start_chunks, int *block_num_par,
                              p_double *dx, t_current_c *jay, int interpolation,
                              p_double dt, int nblocks,
                              int chunk_size, int *n_hole, int *n_comm_ext,
                              int *n_recv_loc, cudaStream_t stream, size_t shmem_size ){

  switch( interpolation ) {
    case p_linear:
      switch( p_x_dim ) {
        case 1:
          cudaFuncSetAttribute( advance_deposit_1d_s1<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_1d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          cudaFuncSetAttribute( advance_deposit_2d_s1<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_2d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          cudaFuncSetAttribute( advance_deposit_3d_s1<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_3d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
    case p_quadratic:
      switch( p_x_dim ) {
        case 1:
          cudaFuncSetAttribute( advance_deposit_1d_s2<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_1d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          cudaFuncSetAttribute( advance_deposit_2d_s2<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_2d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          cudaFuncSetAttribute( advance_deposit_3d_s2<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_3d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
    case p_cubic:
      switch( p_x_dim ) {
        case 1:
          cudaFuncSetAttribute( advance_deposit_1d_s3<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_1d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          cudaFuncSetAttribute( advance_deposit_2d_s3<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_2d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          cudaFuncSetAttribute( advance_deposit_3d_s3<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_3d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
    case p_quartic:
      switch( p_x_dim ) {
        case 1:
          cudaFuncSetAttribute( advance_deposit_1d_s4<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_1d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          cudaFuncSetAttribute( advance_deposit_2d_s4<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_2d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          cudaFuncSetAttribute( advance_deposit_3d_s4<__BLOCK_SIZE__>,
              cudaFuncAttributeMaxDynamicSharedMemorySize, int(shmem_size) );

          advance_deposit_3d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, shmem_size, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
  }

  chkLastErr("advance_deposit");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}

#endif
