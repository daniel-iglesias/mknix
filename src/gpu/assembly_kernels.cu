#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

#include "assembly_kernels.h"
#define WARPSIZE 32 //just in case
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include "cuda_helper.h"

/**
 * Performs atomic addition on double precision.
 * @param  {[type]} double* address  memory position for the atomic operation
 * @param  {[type]} double  val      value to add atomically
 */
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

template <typename T>
__global__ void k_assemble_global_matrix(T* global_matrix,
                                        int* full_map,
                                        int num_points,
                                        int support_node_size,
                                        int number_elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= num_points * support_node_size * support_node_size) return;
  int map_pos = full_map[tid];
  atomicAdd(&global_matrix[map_pos], 1.0);
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
  * @brief Allocates and Creates a Map for a global sparse matrix in the Compressed
  * Column Storage format (CCS) as explained in http://www.netlib.org/utk/people/JackDongarra/etemplates/node373.html
  * This function also allocates all the supplementary arrays for the sparse format.
  * the allocation of the values array is not performed here as depends the data types and wether atomics will be used.
  * @param full_map Reference to the global matrix positions, accesible by [row_id * totcols + col_id ]. lhs.
  * @param row_ind. Array with the col id of each element. lhs.
  * @param col_ptr. Array with the memory position where each row starts. lhs.
  * @param presence_matrix Reference to the global matrix counter where if an element is nonzero will have a counter greater than zero. rhs.
  * @param number_rows. Integer with the total number of rows of the matrix. rhs.
  * @param number_columns. Integer with the total number of columns of the matrix. rhs.
  * @return void. lhs
  **/
    bool build_CCS_sparse_matrix_from_map(std::vector<int> &full_map,
                                          std::vector<int> &row_ind,
                                          std::vector<int> &col_ptr,
                                          int *presence_matrix,
                                          int number_rows,
                                          int number_columns)//Compressed Column Storage
    {
      //calculates total memory needed to allocate memory
        int total_elements = 0;
        for(int i = 1; i < number_rows * number_columns; i++)
            if(presence_matrix[i] > 0) total_elements++;

        row_ind.resize(total_elements);
        full_map.resize(number_rows * number_columns);
        col_ptr.resize(number_rows + 1);

        //prefix sum scan will give us the full map
        full_map[0] = 0;
        for(int i = 1; i < number_rows * number_columns; i++)
                full_map[i] = full_map[i-1] + presence_matrix[i-1];

        for(int j = 1; j < number_rows; j++){
          for(int i = 1; i < number_columns; i++){
            if(presence_matrix[j * number_rows + i] > 0){
              row_ind[full_map[j * number_rows + i]] = i;//row id of every element
            }

          }
          col_ptr[j] = full_map[j * number_rows];//pointers to start of every col
        }

        col_ptr[number_columns] = total_elements;//convention

        return true;
    }

  /**
  * @brief Creates a Map for a global sparse matrix in either the Compressed
  * Column Storage format (CCS) or Compressed Row Storage(CRS) as explained
  * in http://www.netlib.org/utk/people/JackDongarra/etemplates/node373.html
  * This function also allocates all the supplementary arrays for the sparse format.
  * the allocation of the values array is not performed here as depends the data types and wether atomics will be used.
  * @param full_map Reference to the global matrix positions, accesible by [row_id * totcols + col_id ]. lhs.
  * @param row_ind. Array with the vec_ind of each element. lhs.
  * @param cvec_ptr. Array with the memory position where each cvec starts. lhs.
  * @param presence_matrix Reference to the global matrix counter where if an element is nonzero will have a counter greater than zero. rhs.
  * @param all_point_nodes Reference to the array where each point keeps the relation of support nodes is using. rhs.
  * @param number_rows. Integer with the total number of rows of the matrix. rhs.
  * @param number_columns. Integer with the total number of columns of the matrix. rhs.
  * @return void. lhs
  **/
  bool map_global_matrix(std::vector<int> &full_map,
                        std::vector<int> &vec_ind,
                        std::vector<int> &cvec_ptr,
                        int *presence_matrix,
                        int number_rows,
                        int number_columns)
{
  /*build_CRS_sparse_matrix_from_map(full_map,
                                   vec_ind,
                                   cvec_ptr,
                                   presence_matrix,
                                   number_rows,
                                   number_columns);*/
  build_CCS_sparse_matrix_from_map(full_map,
                                  vec_ind,
                                  cvec_ptr,
                                  presence_matrix,
                                  number_rows,
                                  number_columns);

  return true;

}



  /**
   * @brief Allocates and Creates a Map for a global sparse matrix in the Compressed
   * Row Storage format (CRS) as explained in http://www.netlib.org/utk/people/JackDongarra/etemplates/node373.html
   * This function also allocates all the supplementary arrays for the sdparse format.
   * the allocation of the values array is not performed here as depends the data types and wether atomics will be used.
   * @param full_map Reference to the global matrix positions, accesible by [row_id * totcols + col_id ]. lhs.
   * @param col_ind. Array with the col id of each element. lhs.
   * @param row_ptr. Array with the memory position where each row starts. lhs.
   * @param presence_matrix Reference to the global matrix counter where if an element is nonzero will have a counter greater than zero. rhs.
   * @param number_rows. Integer with the total number of rows of the matrix. rhs.
   * @param number_columns. Integer with the total number of columns of the matrix. rhs.
   * @return void. lhs
   **/
  bool build_CRS_sparse_matrix_from_map(std::vector<int> &full_map,
                                                        std::vector<int> &col_ind,
                                                        std::vector<int> &row_ptr,
                                                        int *presence_matrix,
                                                        int number_rows,
                                                        int number_columns)//Compressed Row Storage
{
  //calculates total memory needed to allocate memory
    int total_elements = 0;
    for(int i = 1; i < number_rows * number_columns; i++)
        if(presence_matrix[i] > 0) total_elements++;

    col_ind.resize(total_elements);
    full_map.resize(number_rows * number_columns);
    row_ptr.resize(number_rows + 1);

    //prefix sum scan will give us the full map
    full_map[0] = 0;
    for(int i = 1; i < number_rows * number_columns; i++)
            full_map[i] = full_map[i-1] + presence_matrix[i-1];

    for(int i = 1; i < number_rows; i++){
      for(int j = 1; j < number_columns; j++){
        if(presence_matrix[i*number_columns + j] > 0){
          col_ind[full_map[i*number_columns + j]] = j;//col id of every element
        }

      }
      row_ptr[i] = full_map[i*number_columns];//pointers to start of every row
    }

    row_ptr[number_rows] = total_elements;//convention
    return true;
}

template bool gpu_assemble_global_matrix <float>(float *,
                                                int*,
                                                int,
                                                int,
                                                int,
                                                int,
                                                cudaStream_t stream);
template bool gpu_assemble_global_matrix <double> (double *,
                                                int*,
                                                int,
                                                int,
                                                int,
                                                int,
                                                cudaStream_t stream);

#endif //ASSEMBLY_KERNELS_H
