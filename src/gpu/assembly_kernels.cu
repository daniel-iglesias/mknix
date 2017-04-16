#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

#include "assembly_kernels.h"
#define WARPSIZE 32 //just in case
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include "cuda_helper.h"

template <typename T>
__global__ void k_assemble_global_vector(T* global_vector,
                                        int* full_map,
                                        T* local_matrices_array,
                                        int num_points,
                                        int support_node_size,
                                        int number_elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= num_points * support_node_size * support_node_size) return;
  T value = local_matrices_array[tid];
  int map_pos = full_map[tid];
  atomicAdd(&global_vector[map_pos], value);
}


template <typename T>
__global__ void k_assemble_global_matrix(T* global_matrix,
                                        int* full_map,
                                        T* local_matrices_array,
                                        int num_points,
                                        int support_node_size,
                                        int number_elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= num_points * support_node_size * support_node_size) return;
  T value = local_matrices_array[tid];
  int map_pos = full_map[tid];
  atomicAdd(&global_matrix[map_pos], value);
}
/**
 * Makes a direct assembly of the global matrix from map
 * @param  {[type]} T *global_matrix        global sparse matrix
 * @param  {[type]} int *full_map           Array mapping positions in global matrix
 * @param  {[type]} int num_points          number of gauss points
 * @param  {[type]} int support_node_size   number of nodes support for each gausspoint
 * @param  {[type]} int number_elements     number of elements in the sparse matrix
 * @param  {[type]} int threads_per_block   threads per cuda block
 * @param  {[type]} cudaStream_t stream     cuda stream executing the kernel
 */
  template <typename T>
  bool gpu_assemble_global_matrix(T *global_matrix,
                                  int* full_map,
                                  T* local_matrices_array,
                                  int num_points,
                                  int support_node_size,
                                  int number_elements,
                                  int threads_per_block,
                                  cudaStream_t stream)
  {

    int size = num_points * support_node_size * support_node_size;
    dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
    dim3 blockDim = dim3(threads_per_block,1,1);

    k_assemble_global_matrix<<<gridDim, blockDim, 0, stream>>>(global_matrix,
                                                              full_map,
                                                              local_matrices_array,
                                                              num_points,
                                                              support_node_size,
                                                              number_elements);

    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess)
    {
      std::cout << "GPU.CudaError k_assemble_global_matrix"
      "{\"cuda_error:\"" << cudaError <<
      ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
      ",\"array_size:\"" << size <<
      ",\"gridDim.x:\"" << gridDim.x <<
      ",\"threads_per_block:\"" << threads_per_block <<
      ",\"stream:\"" << stream << "}" << std::endl;

      return false;
    }
    return true;
  }

  /**
   * Makes a direct assembly of the global vector from map
   * @param  {[type]} T *global_vector        global sparse matrix
   * @param  {[type]} int *full_map           Array mapping positions in global matrix
   * @param  {[type]} int num_points          number of gauss points
   * @param  {[type]} int support_node_size   number of nodes support for each gausspoint
   * @param  {[type]} int number_elements     number of elements in the sparse matrix
   * @param  {[type]} int threads_per_block   threads per cuda block
   * @param  {[type]} cudaStream_t stream     cuda stream executing the kernel
   */
    template <typename T>
    bool gpu_assemble_global_vector(T *global_matrix,
                                    int* vector_positions,
                                    T* local_vector,
                                    int num_points,
                                    int support_node_size,
                                    int number_points,
                                    int threads_per_block,
                                    cudaStream_t stream)
    {

      int size = num_points * support_node_size * support_node_size;
      dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
      dim3 blockDim = dim3(threads_per_block,1,1);

      k_assemble_global_vector<<<gridDim, blockDim, 0, stream>>>(global_matrix,
                                                                vector_map,
                                                                local_vector_array,
                                                                num_points,
                                                                support_node_size,
                                                                number_points);

      cudaError_t cudaError = cudaGetLastError();
      if (cudaError != cudaSuccess)
      {
        std::cout << "GPU.CudaError k_assemble_global_vector"
        "{\"cuda_error:\"" << cudaError <<
        ",\"cuda_error_message\":" << "\"" << cudaGetErrorString(cudaError) << "\"" <<
        ",\"array_size:\"" << size <<
        ",\"gridDim.x:\"" << gridDim.x <<
        ",\"threads_per_block:\"" << threads_per_block <<
        ",\"stream:\"" << stream << "}" << std::endl;

        return false;
      }
      return true;
    }



template bool gpu_assemble_global_matrix <float>(float *,
                                                int*,
                                                float *,
                                                int,
                                                int,
                                                int,
                                                int,
                                                cudaStream_t stream);
template bool gpu_assemble_global_matrix <double> (double *,
                                                  int*,
                                                  double *,
                                                  int,
                                                  int,
                                                  int,
                                                  int,
                                                  cudaStream_t stream);

#endif //ASSEMBLY_KERNELS_H
