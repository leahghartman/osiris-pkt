// system related funcs for cuda

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "os-param-cuda.h"
#include "os-sys-cuda.h"

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void __abort_program( const char *file, const int line, bool abort) {
  fprintf(stderr, "ABORTING: returned from file %s:line %d\n", file, line);
  if (abort) {
    exit(-1);
  }
}


// --------------------------------------------------------------------------------------
// Assumes single device when calling cudaDeviceReset(); and exit(code);
// In some cases a more lengthy program clean up / termination may be needed
// --------------------------------------------------------------------------------------
void checkError(cudaError_t code, char const * func, const char *file, const int line,
                bool abort) {
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    fprintf(stderr, "CUDA error returned from \"%s\" at %s:%d, Error code: %d (%s)\n",
            func, file, line, code, errorMessage);
    if (abort) {
      cudaDeviceReset();
      exit(code);
    }
  }
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void checkMallocError( void *x, char const * func, const char *file, const int line,
                       bool abort) {
  if ( !x ) {
    fprintf(stderr, "ERROR calling malloc from \"%s\" at %s:%d\n", func, file, line );
    if (abort) {
      abort_program();
    }
  }
}


// --------------------------------------------------------------------------------------
// Assumes single device when calling cudaDeviceReset(); and exit(code);
// In some cases a more lengthy program clean up / termination may be needed
// --------------------------------------------------------------------------------------
void checkErrorCustom(cudaError_t code, char const * func, const char *file, const int line,
                      bool abort, const char *cstmMsg) {
  if (code != cudaSuccess) {
    const char *errorMessage = cudaGetErrorString(code);
    fprintf(stderr, "%s Error code: %d (%s)\n", cstmMsg, code, errorMessage);
    if (abort) {
      cudaDeviceReset();
      exit(code);
    }
  }
}


// --------------------------------------------------------------------------------------
// Assumes single device when calling cudaDeviceReset(); and exit(code);
// In some cases a more lengthy program clean up / termination may be needed
// --------------------------------------------------------------------------------------
void checkErrorSetDevice(cudaError_t code, char const * func, const char *file, const int line,
                         bool abort, int device_count) {
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    fprintf(stderr, "CUDA error returned from \"%s\" at %s:%d, Error code: %d (%s)\n",
            func, file, line, code, errorMessage);

    if (code==cudaErrorInvalidDevice){
      fprintf(stderr, "Make sure the number of ranks per node you requested doesn't "
                      "exceed the number of GPUs per node, which is %d\n", device_count);
    }

    if (abort) {
      cudaDeviceReset();
      exit(code);
    }
  }
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void checkLastError(char const * func, const char *file, const int line, bool abort) {
  cudaError_t code = cudaGetLastError();
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    fprintf(stderr, "CUDA error returned from \"%s\" at %s:%d, Error code: %d (%s)\n",
            func, file, line, code, errorMessage);
    if (abort) {
      cudaDeviceReset();
      exit(code);
    }
  }
}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
void strmem( char *dest, size_t mem_size ) {

  if ( mem_size > 1073741824 ) {
    sprintf( dest, "%zu GB", mem_size / 1073741824 );
  }
  else if ( mem_size > 1048576 ) {
    sprintf( dest, "%zu MB", mem_size / 1048576 );
  }
  else if ( mem_size > 1024 ) {
    sprintf( dest, "%zu kB", mem_size / 1024 );
  }
  else {
    if( mem_size == 1 ) {
      sprintf( dest, "1 byte" );
    }
    else {
      sprintf( dest, "%zu bytes", mem_size );
    }
  }

}


// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// void checkStatus( char const * func, const char *file, const int line, int status ) {

//     if ( ! code ) {
//     fprintf(stderr, "Error returned from \"%s\" at %s:%d \n",
//             func, file, line, code );
//     if (abort) {
//       cudaDeviceReset();
//       exit(code);
//     }
//   }
// }
