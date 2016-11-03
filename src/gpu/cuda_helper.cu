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
//#include "device_helper.h"



/**
 * allocate GPU memory
 * @param  {[type]} T*                    array         array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 */
//template <typename T>
bool allocate_gpu_array(double *array,int size)
{
    CudaSafeCall(cudaMalloc((void**)&array, size * sizeof(double)));
    return true;
}
bool allocate_gpu_array(float *array,int size)
{
    CudaSafeCall(cudaMalloc((void**)&array, size * sizeof(float)));
    return true;
}

bool allocate_gpu_array(int *array,int size)
{
    CudaSafeCall(cudaMalloc((void**)&array, size * sizeof(int)));
    return true;
}


/**
 * copies data from CPU(host) array to GPU(device) array
 * @param  {[type]} T*                    gpu_array        device array of values
 * @param  {[type]} T*                    cpu_array        host array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 */
template <typename T>
bool copy_to_gpu(T *gpu_array, T*cpu_array, int size)
{
    CudaSafeCall(cudaMemcpy(gpu_array, cpu_array,size * sizeof(T), cudaMemcpyDeviceToHost));
    return true;
}

/**
 * copies data from const CPU(host) array to GPU(device) array
 * @param  {[type]} T*                    gpu_array        device array of values
 * @param  {[type]} T*                    cpu_array        host const array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 */
template <typename T>
bool copy_to_gpu(T *gpu_array,const T*cpu_array, int size)
{
    CudaSafeCall(cudaMemcpy(gpu_array, cpu_array,size * sizeof(T), cudaMemcpyDeviceToHost));
    return true;
}

/**
 * copies data from GPU(device) array to CPU(host) array
 * @param  {[type]} T*                    cpu_array        host array of values
 * @param  {[type]} T*                    gpu_array        device array of values
 * @param  {[type]} T   value             dividend
 * @param  {[type]} int size              size of the array
 */
template <typename T>
bool copy_from_gpu(T *cpu_array, T*gpu_array, int size)
{
    CudaSafeCall(cudaMemcpy(cpu_array, gpu_array,size * sizeof(T), cudaMemcpyHostToDevice));
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
bool divide_array_by_value(T*  array,
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
bool exclusive_scan(int *input_first,
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
size_t free_memory_gpu()
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
__global__ void k_init_array_to_value(double *array,
                                      double value,
                                      int size)
{
    int tidx = threadIdx.x + blockDim.x * blockIdx.x;
    if (tidx >= size) return;

    array[tidx] = value;
}
__global__ void k_init_array_to_value(float *array,
                                      float value,
                                      int size)
{
    int tidx = threadIdx.x + blockDim.x * blockIdx.x;
    if (tidx >= size) return;

    array[tidx] = value;
}
__global__ void k_init_array_to_value(int *array,
                                      int value,
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
 bool init_array_to_value(T **array,
                                      T value,
                                      int size,
                                      int threads_per_block,
                                      cudaStream_t stream)
 {
     std::cout<<"Initializing array of "<< size << " elements with value " << value << std::endl;
     dim3 gridDim = dim3((size + threads_per_block-1) / threads_per_block, 1, 1);
     dim3 blockDim = dim3(threads_per_block,1,1);

     k_init_array_to_value<<<gridDim, blockDim, 0, stream>>>(*array,
                                   value,
                                   size);

     cudaError_t cudaError = cudaGetLastError();
     if (cudaError != cudaSuccess)
     {
       std::cout << "GPU.CudaError k_init_array_to_value"
                     "{\"cuda_error: \"" << cudaError <<
                     ",\"cuda_error_message \":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
                     ",\"array_size: \"" << size <<
                     ",\"value: \"" << value <<
                     ",\"stream: \"" << stream <<
                     ",\"gridDim.x: \"" << gridDim.x <<
                     ",\"threads_per_block: \"" << threads_per_block << "}";
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
 bool init_index_array(int *index_array,
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
 bool multiply_arrays_by_element(T *array,
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
 void print_gpu_memory()
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





 /////////////////////////////////////////////////////////////////////////////
 void checkGPUMemory()
 { size_t free_byte0 ;
 	size_t total_byte0 ;
 	if ( cudaSuccess != cudaMemGetInfo( &free_byte0, &total_byte0 ) ){
 				std::cout <<"Error: cudaMemGetInfo fails" << std::endl; exit(1);}
 	double free_db0 = ((double)free_byte0)/(1024.0*1024.0) ;
 	double total_db0 = ((double)total_byte0)/(1024.0*1024.0) ;
 	double used_db0 = total_db0 - free_db0 ;
 	double memperc = (100 * used_db0)/total_db0;
 	std::cout << std::fixed << std::setprecision(1);
 	std::cout << std::endl << "---- GPU memory usage : "<< memperc << "% ---" << std::endl;
 	std::cout << " used = " << used_db0 <<" MB" << std::endl;
 	std::cout << " free = " << free_db0 <<" MB" << std::endl;
 	std::cout << " total = " << total_db0 << " MB" << std::endl;
 	std::cout << "----------------------------------" << std::endl;}


 void cudaTick(cudaClock *ck)
 {
   cudaEventCreate(&(ck->start));
   cudaEventCreate(&(ck->stop));
   cudaEventRecord(ck->start,0);
 }

 double cudaTock(cudaClock *ck)
 {//modified to suit mknix
   cudaEventRecord(ck->stop, 0);
   cudaEventSynchronize(ck->stop);
   cudaEventElapsedTime(&(ck->gpu_time), ck->start, ck->stop);;
 //  std::cout << "GPU clock measured "<<  ck->gpu_time *1000.0f << " microseconds" << std::endl;
   cudaEventDestroy(ck->start); //cleaning up a bit
   cudaEventDestroy(ck->stop);
   ck->last_measure = ck->gpu_time;//not really suing it yet
   ck->gpu_time = 0.0f;
   return ck->last_measure * 1000.0;
 }

/* template bool allocate_gpu_array<double>(double *array,int size);
 template bool allocate_gpu_array<float>(float *array,int size);
 template bool allocate_gpu_array<int>(int *array,int size);*/

 template bool init_array_to_value<double>(double **array,double value,int size,int threads_per_block,cudaStream_t stream);
 template bool init_array_to_value<float>(float **array,float value,int size,int threads_per_block,cudaStream_t stream);
 template bool init_array_to_value<int>(int **array,int value,int size,int threads_per_block,cudaStream_t stream);

 template bool copy_to_gpu<double>(double *gpu_array, double *cpu_array, int size);
 template bool copy_to_gpu<float>(float *gpu_array, float *cpu_array, int size);
 template bool copy_to_gpu<int>(int *gpu_array, int *cpu_array, int size);

 template bool copy_to_gpu<double>(double *gpu_array, const double *cpu_array, int size);
 template bool copy_to_gpu<float>(float *gpu_array, const float *cpu_array, int size);
 template bool copy_to_gpu<int>(int *gpu_array, const int *cpu_array, int size);

 template bool copy_from_gpu<double>(double *cpu_array, double *gpu_array, int size);
 template bool copy_from_gpu<float>(float *cpu_array, float *gpu_array, int size);
 template bool copy_from_gpu<int>(int *cpu_array, int *gpu_array, int size);
