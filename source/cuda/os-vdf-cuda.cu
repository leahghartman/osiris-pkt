// if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

#include <stdio.h>
#include <cuda_runtime.h>
#include <algorithm>
using namespace std;

#include "os-param-cuda.h"
#include "os-mem-transfer-manager-cuda.h"
#include "os-vdf-cuda.h"
#include "os-sys-cuda.h"

#define __TEMPLATE__

//**************************************** 1D ********************************************

#define T_FN t_f0
#define DIMS n0
#define DIM_ARGS int nx0
#define COMP c
#define COMP_ARGS int c
#define OBJ_ASSGN n0(obj.n0)
#define OBJ_COPY n0 = obj.n0;
#define INP_ASSGN n0(nx0)
#define INP_COPY n0 = nx0;
#define SIZE n0

#include __FILE__

//**************************************** 2D ********************************************

#define T_FN t_f1
#define DIMS n0, n1
#define DIM_ARGS int nx0, int nx1
#define COMP n0*i + c
#define COMP_ARGS int c, int i
#define OBJ_ASSGN n0(obj.n0), n1(obj.n1)
#define OBJ_COPY n0 = obj.n0; n1 = obj.n1;
#define INP_ASSGN n0(nx0), n1(nx1)
#define INP_COPY n0 = nx0; n1 = nx1;
#define SIZE n0 * n1

#include __FILE__

//**************************************** 3D ********************************************

#define T_FN t_f2
#define DIMS n0, n1, n2
#define DIM_ARGS int nx0, int nx1, int nx2
#define COMP ( n1*j + i ) * n0 + c
#define COMP_ARGS int c, int i, int j
#define OBJ_ASSGN n0(obj.n0), n1(obj.n1), n2(obj.n2)
#define OBJ_COPY n0 = obj.n0; n1 = obj.n1; n2 = obj.n2;
#define INP_ASSGN n0(nx0), n1(nx1), n2(nx2)
#define INP_COPY n0 = nx0; n1 = nx1; n2 = nx2;
#define SIZE n0 * n1 * n2

#include __FILE__

//**************************************** 4D ********************************************

#define T_FN t_f3
#define DIMS n0, n1, n2, n3
#define DIM_ARGS int nx0, int nx1, int nx2, int nx3
#define COMP ( ( n2*k + j ) * n1 + i ) * n0 + c
#define COMP_ARGS int c, int i, int j, int k
#define OBJ_ASSGN n0(obj.n0), n1(obj.n1), n2(obj.n2), n3(obj.n3)
#define OBJ_COPY n0 = obj.n0; n1 = obj.n1; n2 = obj.n2; n3 = obj.n3;
#define INP_ASSGN n0(nx0), n1(nx1), n2(nx2), n3(nx3)
#define INP_COPY n0 = nx0; n1 = nx1; n2 = nx2; n3 = nx3;
#define SIZE n0 * n1 * n2 * n3

#include __FILE__

#else

template<class T>
__host__ __device__ T_FN<T>::T_FN( const T_FN &obj ) : f(obj.f), g(obj.g),
  unique_buff(obj.unique_buff), OBJ_ASSGN {}

template<class T>
__host__ __device__ T_FN<T>& T_FN<T>::operator=( const T_FN &obj ) {
  OBJ_COPY
  f = obj.f;
  g = obj.g;
  unique_buff = obj.unique_buff;
  return *this;
}

template<class T>
__host__ T_FN<T>::T_FN( DIM_ARGS, T *f_, bool on_dev ) :
  unique_buff(false), INP_ASSGN {
  if (on_dev) {
    // do nothing yet
  }
  else {
    f = f_;
    g = f;
  }
}

template<class T>
__host__ T_FN<T>::T_FN( DIM_ARGS, bool on_dev ) :
  unique_buff(true), INP_ASSGN {
  if (on_dev) {
    // do nothing yet
  }
  else {
    f = new T[ SIZE ];
    g = f;
  }
}

template<class T>
__host__ __device__ void T_FN<T>::init( DIM_ARGS, T *f_, bool on_dev ) {

  INP_COPY
  f = f_;
  g = f;

}

template<class T>
__host__ void T_FN<T>::copy_to_device( T_FN *d_obj ) const {
  T_FN h_obj;
  chkErr(cudaMemcpy(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost));
  chkErr(cudaMemcpy(h_obj.f, f, sizeof(T)*size(), cudaMemcpyHostToDevice));
}

template<class T>
__host__ void T_FN<T>::copy_to_device( T *d_f ) const {
  chkErr(cudaMemcpy(d_f, f, sizeof(T)*size(), cudaMemcpyHostToDevice));
}

template<class T>
__host__ T* T_FN<T>::malloc_copy_to_device( T_FN *d_obj ) const {
  T_FN h_obj = *this;
  size_t mem_size = sizeof(T) * size();
  chkErr(cudaMalloc(reinterpret_cast<void **>(&(h_obj.f)), mem_size));
  chkErr(cudaMemcpy(h_obj.f, f, mem_size, cudaMemcpyHostToDevice));
  h_obj.g = h_obj.f + (g - f);
  chkErr(cudaMemcpy(d_obj, &h_obj, sizeof(T_FN), cudaMemcpyHostToDevice));
  return h_obj.f;
}

template<class T>
__host__ T* T_FN<T>::buffered_malloc_copy_to_device( T_FN *d_obj, t_mem_transfer_manager *mem_transfer_manager, bool copy_f ) const {
  T_FN h_obj = *this;
  size_t mem_size = sizeof(T) * size();
  chkErr(cudaMalloc(reinterpret_cast<void **>(&(h_obj.f)), mem_size));
  if( copy_f ) mem_transfer_manager->bufferedCudaMemcpy(h_obj.f, f, mem_size, cudaMemcpyHostToDevice);
  h_obj.g = h_obj.f + (g - f);
  mem_transfer_manager->bufferedCudaMemcpy(d_obj, &h_obj, sizeof(T_FN), cudaMemcpyHostToDevice);
  return h_obj.f;
}

template<class T>
__host__ void T_FN<T>::copy_to_device_async( T_FN *d_obj, cudaStream_t stream ) const {
  T_FN h_obj;
  chkErr(cudaMemcpyAsync(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost, stream));
  chkErr(cudaMemcpyAsync(h_obj.f, f, sizeof(T)*size(), cudaMemcpyHostToDevice, stream));
}

template<class T>
__host__ void T_FN<T>::copy_to_device_async( T *d_f, cudaStream_t stream ) const {
  chkErr(cudaMemcpyAsync(d_f, f, sizeof(T)*size(), cudaMemcpyHostToDevice, stream));
}

// template<class T>
// __host__ void T_FN<T>::buffered_copy_to_device( T_FN *d_obj, t_mem_transfer_manager *mem_transfer_manager ) const {
//   T_FN h_obj;
//   mem_transfer_manager->bufferedCudaMemcpy(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost );
//   mem_transfer_manager->bufferedCudaMemcpy(h_obj.f, f, sizeof(T)*size(), cudaMemcpyHostToDevice );
// }

template<class T>
__host__ void T_FN<T>::buffered_copy_to_device( T *d_f, t_mem_transfer_manager *mem_transfer_manager ) const {
  mem_transfer_manager->bufferedCudaMemcpy(d_f, f, sizeof(T)*size(), cudaMemcpyHostToDevice );
}

template<class T>
__host__ T* T_FN<T>::malloc_copy_to_device_async( T_FN *d_obj, cudaStream_t stream ) const {
  T_FN h_obj = *this;
  size_t mem_size = sizeof(T) * size();
  chkErr(cudaMalloc(reinterpret_cast<void **>(&(h_obj.f)), mem_size));
  chkErr(cudaMemcpyAsync(h_obj.f, f, mem_size, cudaMemcpyHostToDevice, stream));
  h_obj.g = h_obj.f + (g - f);
  chkErr(cudaMemcpyAsync(d_obj, &h_obj, sizeof(T_FN), cudaMemcpyHostToDevice, stream));
  return h_obj.f;
}

template<class T>
__host__ void T_FN<T>::copy_from_device( T_FN *d_obj ) const {
  T_FN h_obj;
  chkErr(cudaMemcpy(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost));
  chkErr(cudaMemcpy(f, h_obj.f, sizeof(T)*size(), cudaMemcpyDeviceToHost));
}

template<class T>
__host__ void T_FN<T>::copy_free_from_device( T_FN *d_obj ) const {
  T_FN h_obj;
  chkErr(cudaMemcpy(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost));
  chkErr(cudaMemcpy(f, h_obj.f, sizeof(T)*size(), cudaMemcpyDeviceToHost));
  chkErr(cudaFree(h_obj.f));
}

template<class T>
__host__ void T_FN<T>::copy_from_device_async( T_FN *d_obj, cudaStream_t stream ) const {
  T_FN h_obj;
  chkErr(cudaMemcpyAsync(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost, stream));
  chkErr(cudaMemcpyAsync(f, h_obj.f, sizeof(T)*size(), cudaMemcpyDeviceToHost, stream));
}

// template<class T>
// __host__ void T_FN<T>::buffered_copy_from_device( T_FN *d_obj, t_mem_transfer_manager *mem_transfer_manager ) const {
//   T_FN h_obj;
//   mem_transfer_manager->bufferedCudaMemcpy(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost );
//   mem_transfer_manager->bufferedCudaMemcpy(f, h_obj.f, sizeof(T)*size(), cudaMemcpyDeviceToHost );
// }


template<class T>
__host__ void T_FN<T>::copy_free_from_device_async( T_FN *d_obj, cudaStream_t stream ) const {
  T_FN h_obj;
  chkErr(cudaMemcpyAsync(&h_obj, d_obj, sizeof(T_FN), cudaMemcpyDeviceToHost, stream));
  chkErr(cudaMemcpyAsync(f, h_obj.f, sizeof(T)*size(), cudaMemcpyDeviceToHost, stream));
  chkErr(cudaFree(h_obj.f));
}

template<class T>
__host__ void T_FN<T>::copy_from_device( T *d_f ) const {
  chkErr(cudaMemcpy(f, d_f, sizeof(T)*size(), cudaMemcpyDeviceToHost));
}

template<class T>
__host__ void T_FN<T>::copy_free_from_device( T *d_f ) const {
  chkErr(cudaMemcpy(f, d_f, sizeof(T)*size(), cudaMemcpyDeviceToHost));
  chkErr(cudaFree(d_f));
}

template<class T>
__host__ void T_FN<T>::copy_from_device_async( T *d_f, cudaStream_t stream ) const {
  chkErr(cudaMemcpyAsync(f, d_f, sizeof(T)*size(), cudaMemcpyDeviceToHost, stream));
}

template<class T>
__host__ void T_FN<T>::buffered_copy_from_device( T *d_f, t_mem_transfer_manager *mem_transfer_manager ) const {
  mem_transfer_manager->bufferedCudaMemcpy( f, d_f, sizeof(T)*size(), cudaMemcpyDeviceToHost );
}

template<class T>
__host__ void T_FN<T>::copy_free_from_device_async( T *d_f, cudaStream_t stream ) const {
  chkErr(cudaMemcpyAsync(f, d_f, sizeof(T)*size(), cudaMemcpyDeviceToHost, stream));
  chkErr(cudaFree(d_f));
}

template<class T>
__host__ __device__ int T_FN<T>::size( void ) const {
  return SIZE;
}

// template<class T>
// __host__ __device__ T& T_FN<T>::operator()( COMP_ARGS ) {
//   return g[ COMP ];
// }

template<class T>
__host__ __device__ void T_FN<T>::offset( COMP_ARGS ) {
  g += COMP;
}

template<class T>
__host__ void T_FN<T>::cleanup( bool on_dev ) {
  if( ! on_dev ) {
    if (unique_buff) {
      delete[] f;
      f = NULL;
      g = NULL;
    }
  } else {
    T_FN h_obj;
    chkErr(cudaMemcpy(&h_obj, this, sizeof(T_FN), cudaMemcpyDeviceToHost));
    chkErr(cudaFree(h_obj.f));
  }
}


// Declare allowed class templates
template class T_FN<float>;
template class T_FN<double>;
template class T_FN<int>;

#undef T_FN
#undef DIMS
#undef DIM_ARGS
#undef COMP
#undef COMP_ARGS
#undef OBJ_ASSGN
#undef OBJ_COPY
#undef INP_ASSGN
#undef INP_COPY
#undef SIZE

#endif
