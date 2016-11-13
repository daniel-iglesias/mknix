#ifndef CUDA_HELPER_H
#define CUDA_HELPER_H
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cuda.h>
//#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#define CUDA_ERROR_CHECK

//namespace CudaHelper{

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line ){
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );exit( -1 );}
#endif
  return;}
///////////////////////////////////////////////////////////////////////////////
inline void __cudaCheckError( const char *file, const int line ){
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err )
    {  fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );  exit( -1 ); }
    // More careful checking. However, this will affect performance.// Comment away if needed.
    /*err = cudaDeviceSynchronize();
    if( cudaSuccess != err ){fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) ); exit( -1 );}*/
#endif
  return;}

  /*__device__ double atomicAdd(double* address, double val)
  {
   unsigned long long int* address_as_ull = (unsigned long long int*)address;
   unsigned long long int old = *address_as_ull, assumed;
   do {
     assumed = old;
     old = atomicCAS(address_as_ull, assumed,
     __double_as_longlong(val +__longlong_as_double(assumed)));
   } while (assumed != old);
   return __longlong_as_double(old);
 }*/




/////////////////////////////////////////////////////////////////////////////
void checkGPUMemory();

struct cudaClock
{
  cudaEvent_t start, stop;
  float gpu_time, last_measure;
  double elapsedMicroseconds; //to have same as chTimer
};

void cudaTick(cudaClock *ck);

void cudaTock(cudaClock *ck, std::string function_name);
//namespace CudaHelper{

   bool allocate_gpu_array(double *array,int size);
    bool allocate_gpu_array(float *array,int size);
     bool allocate_gpu_array(int *array,int size);

   template <typename T>
   bool copy_to_gpu(T *gpu_array, T* cpu_array, int size);

   template <typename T>
   bool copy_to_gpu(T *gpu_array, const T* cpu_array, int size);

   template <typename T>
   bool copy_from_gpu(T *cpu_array, T*gpu_array, int size);

    template <typename T>
    bool divide_array_by_value(T*  array,
                               T value,
                               int size,
                               int threads_per_block,
                               cudaStream_t stream = 0);

    bool exclusive_scan(int *input_first,
                        int *input_last,
                        int *output,
                        int init);

    size_t free_memory_gpu();

    template <typename T>
    bool init_array_to_value(T *array,
                             T value,
                             int size,
                             int threads_per_block,
                             cudaStream_t stream = 0);

    bool init_index_array(int *index_array,
                          int size,
                          int threads_per_block,
                          cudaStream_t stream = 0);

    template <typename T, typename T2>
    bool multiply_arrays_by_element(T *array,
                                    T2 *array2,
                                    int size,
                                    int threads_per_block,
                                    cudaStream_t stream = 0);

    void print_gpu_memory();


    __global__ void k_init_array_to_value(double *array,
                                          double value,
                                          int size);

//}


#endif //CUDA_HELPER_H
