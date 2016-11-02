#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/inner_product.h>
#include <thrust/reduce.h>

#include "cuda_helper.h"


/**
 * allocate GPU memory
 * @param  {[type]} T*                    array         array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 */
template <typename T>
bool allocate_gpu_array(T *array,int size)
{
    cudaMalloc((void**)&array, size * sizeof(T));
    return true;
}


/**
 * Divides each element of the array by the given value
 * @param  {[type]} T*                    array         array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 */
template <typename T>
__global__ void k_divide_array_by_value(T *array,
                                        T value,
                                        int size)
{
    int tidx = threadIdx.x + blockDim.x * blockIdx.x;
    if (tidx >= size) return;

    array[tidx] /= value;
}
/**
 * Divides each element of the array by the given value
 * @param  {[type]} T*                    array         array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 * @param  {[type]} int threads_per_block threads per cuda block
 */
template <typename T>
bool CudaHelper::divide_array_by_value(T*  array,
                                       T value,
                                       int size,
                                       int threads_per_block,
                                       cudaStream_t stream)
{
    dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
    dim3 blockDim = dim3(threads_per_block,1,1);

    k_divide_array_by_value<<<gridDim, blockDim, 0, stream>>>(array,
                                                   value,
                                                   size);

    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess)
    {
      std::cout << "GPU.CudaError k_divide_array_by_value"
                    "{\"cuda_error:\"" << cudaError <<
                    ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
                    ",\"array_size:\"" << size <<
                    ",\"stream:\"" << stream <<
                    ",\"threads_per_block:\"" << threads_per_block << "}";
      return false;
    }
    return true;
}


/**
 * Performs the exclusive scan starting from first to last, counting from init
 * @param  {[type]} int *input_first  First element of the array to start the exclusive scan
 * @param  {[type]} int *input_last   Last element to be used in the exclusive scan
 * @param  {[type]} int *output       Output array
 * @param  {[type]} int init          Initial value
 */
bool CudaHelper::exclusive_scan(int *input_first,
                                int *input_last,
                                int *output,
                                int init)
{
    thrust::plus<int> binary_op;
  	thrust::exclusive_scan(thrust::device,
                           input_first,
                           input_last,
                           output,
                           init,
                           binary_op);

  cudaError_t cudaError = cudaGetLastError();
  if (cudaError != cudaSuccess)
  {
    std::cout << "GPU.CudaError thrust::exclusive_scan"
                  "{\"cuda_error:\"" << cudaError <<
                  ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" << "}";
    return false;
  }
  return true;
}

/**
 * Gives the free memory of the gpu in bytes
 * @return {[type]} free bytes in the gpu
 */
size_t CudaHelper::free_memory_gpu()
{
  size_t free_byte;
  size_t total_byte;
  cudaMemGetInfo( &free_byte, &total_byte);
  return free_byte;
}

/**
 * Cuda kernel that initializes every element of an array to a value
 * @param  {[type]} T   *array        array of values
 * @param  {[type]} T   value         init value
 * @param  {[type]} int size          size of the array
 */
template <typename T>
__global__ void k_init_array_to_value(T *array,
                                      T value,
                                      int size)
{
    int tidx = threadIdx.x + blockDim.x * blockIdx.x;
    if (tidx >= size) return;

    array[tidx] = value;
}
/**
 * Initializes every element of an array to a value
 * @param  {[type]} T   *array        array of values
 * @param  {[type]} T   value         init value
 * @param  {[type]} int size          size of the array
 * @param  {[type]} int threads_per_block threads per cuda block
 */
 template <typename T>
 bool CudaHelper::init_array_to_value(T *array,
                                      T value,
                                      int size,
                                      int threads_per_block,
                                      cudaStream_t stream)
 {

     dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
     dim3 blockDim = dim3(threads_per_block,1,1);

     k_init_array_to_value<<<gridDim, blockDim, 0, stream>>>(array,
                                   value,
                                   size);

     cudaError_t cudaError = cudaGetLastError();
     if (cudaError != cudaSuccess)
     {
       std::cout << "GPU.CudaError k_init_array_to_value"
                     "{\"cuda_error:\"" << cudaError <<
                     ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
                     ",\"array_size:\"" << size <<
                     ",\"value:\"" << value <<
                     ",\"stream:\"" << stream <<
                     ",\"threads_per_block:\"" << threads_per_block << "}";
       return false;
     }
     return true;
 }

 /**
  * Initializes an array to the value of its index
  * @param  {[type]} int *index_array  array of values
  * @param  {[type]} int size          size of the array
  */
 __global__ void k_init_index_array(int *index_array,
                                    int size)
 {
     int tidx = threadIdx.x + blockDim.x * blockIdx.x;
     if (tidx >= size) return;

     index_array[tidx] = tidx;
 }

 /**
  * Initializes an array to the value of its index
  * @param  {[type]} int *index_array  array of values
  * @param  {[type]} int size          size of the array
  * @param  {[type]} int threads_per_block threads per cuda block
  */
 bool CudaHelper::init_index_array(int *index_array,
                                   int size,
                                   int threads_per_block,
                                   cudaStream_t stream)
 {

     dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
     dim3 blockDim = dim3(threads_per_block,1,1);

     k_init_index_array<<<gridDim, blockDim, 0, stream>>>(index_array,
                                                          size);

     cudaError_t cudaError = cudaGetLastError();
     if (cudaError != cudaSuccess)
     {
       std::cout << "GPU.CudaError k_init_index_array"
                     "{\"cuda_error:\"" << cudaError <<
                     ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
                     ",\"array_size:\"" << size <<
                     ",\"threads_per_block:\"" << threads_per_block << "}";
       return false;
     }
     return true;

 }


 template <typename T, typename T2>
 __global__ void k_multiply_arrays_by_element(T *array,
                                              T2 *array2,
                                              int size)

 {
   int tidx = threadIdx.x + blockDim.x * blockIdx.x;
   if (tidx >= size) return;

   array[tidx] *= array2[tidx];

 }

 template <typename T, typename T2>
 bool CudaHelper::multiply_arrays_by_element(T *array,
                                 T2 *array2,
                                 int size,
                                 int threads_per_block,
                                 cudaStream_t stream)

 {
   dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
   dim3 blockDim = dim3(threads_per_block,1,1);

   k_multiply_arrays_by_element<<<gridDim, blockDim,0,stream>>>(array,
                                                                array2,
                                                                size);

   cudaError_t cudaError = cudaGetLastError();
   if (cudaError != cudaSuccess)
   {
     std::cout << "GPU.CudaError k_multiply_arrays_by_element"
                   "{\"cuda_error:\"" << cudaError <<
                   ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
                   ",\"array_size:\"" << size <<
                   ",\"threads_per_block:\"" << threads_per_block << "}";
     return false;
   }
   return true;

 }

 /**
  * Prints GPU memory information
  */
 void CudaHelper::print_gpu_memory()
 {
   size_t free_byte0;
   size_t total_byte0;

   if (cudaSuccess != cudaMemGetInfo( &free_byte0, &total_byte0 ) )
   {
     	 std::cout <<"GPU.CudaError cudaMemGetInfo fails";
        return;
   }

   double free_db0 = ((double)free_byte0)/(1024.0*1024.0) ;
   double total_db0 = ((double)total_byte0)/(1024.0*1024.0) ;
   double used_db0 = total_db0 - free_db0 ;

   double memperc = (100 * used_db0)/total_db0;

   std::cout << "---- GPU memory usage : "<< memperc << "% ---" << std::endl
             << " used = " << used_db0 <<" MB" << std::endl
             << " free = " << free_db0 <<" MB" << std::endl
             << " total = " << total_db0 << " MB" << std::endl
             << "----------------------------------" << std::endl;
 }
