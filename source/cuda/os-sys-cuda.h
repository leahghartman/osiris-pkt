#ifndef OS_SYS_CUDA_H
#define OS_SYS_CUDA_H

// ----------------------
// Mirror os-sys-multi.f03
// ----------------------

#define p_single float
#define p_double double
#define p_max_dim 3

// -----------------
// os-sys-cuda.cu
// -----------------

void __abort_program( const char *file, const int line, bool abort);
#define abort_program() __abort_program(__FILE__, __LINE__, true)

// To be used around calls that return an error code, ex. cudaDeviceSynchronize or
// cudaMallocManaged
void checkError(cudaError_t code, char const * func, const char *file, const int line,
                bool abort);
#define chkErr(val) checkError((val), #val, __FILE__, __LINE__, true)

void checkErrorCustom(cudaError_t code, char const * func, const char *file, const int line,
                      bool abort, const char *cstmMsg);
#define chkErrCstm(val, cstmMsg) checkErrorCustom((val), #val, __FILE__, __LINE__, true, cstmMsg)

// same thing for malloc
void checkMallocError(void *x, char const * func, const char *file, const int line,
                      bool abort);
#define chkMallocErr(x) checkMallocError((x), #x, __FILE__, __LINE__, true)

void strmem( char *dest, size_t mem_size );

// To be used specifically around one particular cudaSetDevice call, so I can return
// an error message specific to that scenaaario
void checkErrorSetDevice(cudaError_t code, char const * func, const char *file,
                         const int line, bool abort, int device_count);
#define chkErrStDv(val, device_count) checkErrorSetDevice((val), #val, __FILE__, \
                   __LINE__, true, device_count)

// To be used after calls that do not return an error code, ex. kernels to check kernel
// launch errors
void checkLastError(char const * func, const char *file, const int line, bool abort);
#define chkLastErr(func) checkLastError(func, __FILE__, __LINE__, true)
#define chkLastErr_noAbort(func) checkLastError(func, __FILE__, __LINE__, false)

// Gives ceil(a/b)
#define CEIL_DIV(a,b) ((a) + (b) - 1) / (b)


// --------------
// Debug
// --------------

#define PRINT( x ) printf("%s\n", x )
#define PRINT_INT( x ) printf("%s = %d\n", #x, x )
#define PRINT_HEX( x ) printf("%s = %x\n", #x, x )
#define PRINT_FLOAT( x ) printf("%s = %f\n", #x, x )
#define PRINT_EFLOAT( x ) printf("%s = %.5e\n", #x, x )
#define PRINT_BOOL( x ) printf("%s = %d\n", #x, (x) )
#define PRINT_PTR( x ) printf("%s = %p\n", #x, x )
#define PRINT_SIZET( x ) printf("%s = %zu\n", #x, x )

#define DEBUG_BANNER() printf( "@@@@@@@@@@@@@@@@\n@@@@@@@@@@@@@@@@\n" )
#define REACHED() printf("REACHED\n"); exit(0);

#endif
