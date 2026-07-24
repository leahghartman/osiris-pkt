#ifndef OS_VDF_CUDA

#include "os-mem-transfer-manager-cuda.h"
#include <cuda_runtime.h>

// if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

// Workhorse classes for multi-dimensional arrays

#define __TEMPLATE__

//**************************************** 1D ********************************************

#define T_FN t_f0
#define DIMS n0
#define DIM_ARGS int nx0
#define COMP c
#define COMP_ARGS int c
#define INP_COPY n0 = nx0;

#include __FILE__

//**************************************** 2D ********************************************

#define T_FN t_f1
#define DIMS n0, n1
#define DIM_ARGS int nx0, int nx1
#define COMP n0*i + c
#define COMP_ARGS int c, int i
#define INP_COPY n0 = nx0; n1 = nx1;

#include __FILE__

//**************************************** 3D ********************************************

#define T_FN t_f2
#define DIMS n0, n1, n2
#define DIM_ARGS int nx0, int nx1, int nx2
#define COMP ( n1*j + i ) * n0 + c
#define COMP_ARGS int c, int i, int j
#define INP_COPY n0 = nx0; n1 = nx1; n2 = nx2;

#include __FILE__

//**************************************** 4D ********************************************

#define T_FN t_f3
#define DIMS n0, n1, n2, n3
#define DIM_ARGS int nx0, int nx1, int nx2, int nx3
#define COMP ( ( n2*k + j ) * n1 + i ) * n0 + c
#define COMP_ARGS int c, int i, int j, int k
#define INP_COPY n0 = nx0; n1 = nx1; n2 = nx2; n3 = nx3;

#include __FILE__

#define OS_VDF_CUDA
#undef __TEMPLATE__

#else

//----------------------------------------------------------------------------------------

template<class T>
struct T_FN {
public:
  int DIMS;

  T *f; // Pointer to the beginning of the data array
  T *g; // Pointer to certain location of data array (e.g., first non-guard-cell quantity)

  // True if the data array was allocated upon construction of this T_FN object. Used to
  // determine whether or not that buffer should also be cleaned up when this object is cleaned up
  bool unique_buff;

  __host__ __device__ T_FN( void ) {}
  __host__ __device__ T_FN( const T_FN &obj );
  __host__ __device__ T_FN& operator=( const T_FN &obj );
  __host__ T_FN( DIM_ARGS, T *f_, bool on_dev );
  __device__ __forceinline__ T_FN( DIM_ARGS ) { INP_COPY f = NULL; g = NULL; unique_buff = false; };
  __host__ T_FN( DIM_ARGS, bool on_dev );
  __host__ __device__ void init( DIM_ARGS, T *f_, bool on_dev );
  __host__ void copy_to_device( T_FN *d_obj ) const;
  __host__ void copy_to_device( T *d_f ) const;
  __host__ T* malloc_copy_to_device( T_FN *d_obj ) const;
  __host__ void copy_to_device_async( T_FN *d_obj, cudaStream_t stream ) const;
  __host__ void copy_to_device_async( T *d_f, cudaStream_t stream ) const;
  // __host__ void buffered_copy_to_device( T_FN *d_obj, t_mem_transfer_manager *mem_transfer_manager ) const;
  __host__ void buffered_copy_to_device( T *d_f, t_mem_transfer_manager *mem_transfer_manager ) const;
  __host__ T* malloc_copy_to_device_async( T_FN *d_obj, cudaStream_t stream ) const;
  __host__ T* buffered_malloc_copy_to_device( T_FN *d_obj, t_mem_transfer_manager *mem_transfer_manager, bool copy_f ) const;
  __host__ void copy_from_device( T_FN *d_obj ) const;
  __host__ void copy_free_from_device( T_FN *d_obj ) const;
  __host__ void copy_from_device_async( T_FN *d_obj, cudaStream_t stream ) const;
  // __host__ void buffered_copy_from_device( T_FN *d_obj, t_mem_transfer_manager *mem_transfer_manager ) const;
  __host__ void copy_free_from_device_async( T_FN *d_obj, cudaStream_t stream ) const;
  __host__ void copy_from_device( T *d_f ) const;
  __host__ void copy_free_from_device( T *d_f ) const;
  __host__ void copy_from_device_async( T *d_f, cudaStream_t stream ) const;
  __host__ void buffered_copy_from_device( T *d_f, t_mem_transfer_manager *mem_transfer_manager ) const;
  __host__ void copy_free_from_device_async( T *d_f, cudaStream_t stream ) const;
  __host__ __device__ int size( void ) const;
  // __host__ __device__ T& operator()( COMP_ARGS );
  __host__ T& fh( COMP_ARGS ) const { return g[ COMP ]; }
  __device__ __forceinline__ T& fd( COMP_ARGS ) const { return g[ COMP ]; }
  __host__ T* fhp( COMP_ARGS ) const { return g + COMP; }
  __device__ __forceinline__ T* fdp( COMP_ARGS ) const { return g + COMP; }
  __host__ __device__ void offset( COMP_ARGS );
  __host__ void cleanup( bool on_dev = false );
  // __host__ void cleanup_from_device( void );
  __host__ __device__ ~T_FN( void ) {}
};
//----------------------------------------------------------------------------------------

#undef T_FN
#undef DIMS
#undef DIM_ARGS
#undef COMP
#undef COMP_ARGS
#undef INP_COPY

#endif

#endif
