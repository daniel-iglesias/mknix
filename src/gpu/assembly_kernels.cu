#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

#include "assembly_kernels.h"
#define WARPSIZE 32 //just in case
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include "cuda_helper.h"
#include "cpu/structures.h"
//
//
template <typename T>
__device__ T d_interpolate1D(T query_value,
                             T *reference_values,
                             T *sought_values,
                             int init_position,
                             int counter)
{
  bool upper_bounded = false;
  int upper_index = 0;
  for(int i = 0; i < counter; i++){
    //first find the upper_bound
    T this_val = reference_values[init_position + i];
    if(query_value <= this_val){//bounded here
      upper_index = i;
      break;
    }
  }
  if(upper_index == 0) return sought_values[init_position];
  upper_bounded = true;
  if(!upper_bounded) return sought_values[init_position + counter];//not bound found so return last
  //so we have a bound
  int lower_index = upper_index - 1;
  T delta = (query_value - reference_values[init_position + lower_index]) / (reference_values[init_position + upper_index] - reference_values[init_position + lower_index]);
  return delta * sought_values[init_position + upper_index] + (1.0 - delta) * sought_values[init_position + lower_index];
}
//
__device__ double d_getMaterialKappa (int *materials_kappa_counters,
                                      int *materials_kappa_inits,
                                      double* materials_kappa_temps,
                                      double* materials_kappa_values,
                                      int material_id,
                                      double average_temperature)
{
  int n_vals =  materials_kappa_counters[material_id];
  int init_vals = materials_kappa_inits[material_id];
  //if(threadIdx.x == 0 && blockIdx.x == 0)printf("n_vals %d \n", n_vals);
  //if(threadIdx.x == 0 && blockIdx.x == 0)printf("init_vals %d \n", init_vals);
  return d_interpolate1D(average_temperature,
                         materials_kappa_temps,
                         materials_kappa_values,
                         init_vals,
                         n_vals);
}

__device__ double d_getMaterialDensity (MaterialTable *materials,
                                      int material_id)
{
  return materials->density[material_id];
}

__device__ double d_getMaterialCapacity (MaterialTable *materials,
                                      int material_id,
                                      double average_temperature)
{
  int n_vals =  materials->_capacity_counters[material_id];
  int init_vals =  materials->_capacity_inits[material_id];
  return d_interpolate1D(average_temperature,
                       materials->_capacity_temps,
                       materials->_capacity_values,
                       init_vals,
                       n_vals);

}

//////////////////////////////////////////////////////////////////////////////
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
                                        uint* full_map,
                                        uint* node_map,
                                        T* local_matrices_array,
                                        int num_points,
                                        int support_node_size,
                                        int number_elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= num_points * support_node_size * support_node_size) return;
  T value = local_matrices_array[tid];
  int node_pos = node_map[tid];
  int map_pos = full_map[node_pos];
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
                                  uint* full_map,
                                  uint* node_map,
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
                                                              node_map,
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
                                                                vector_positions,
                                                                local_vector,
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
//
template <typename T>
__global__  void k_computeSOACapacityMatrix(T *local_capacity_matrices_array,
                                            T *local_temperatures_array,
                                            T *local_weight_array,
                                            T *local_jacobian_array,
                                            T *local_shapeFun_phis,
                                            int *material_ids,
                                            MaterialTable materials,
                                            int numPoints,
                                            int supportNodeSize)
  {
  int eachPoint = threadIdx.x + blockIdx.x * blockDim.x;
  if(eachPoint >= numPoints) return;

  double avgTemp  = 0.0;
  for(int lnode = 0; lnode < supportNodeSize; lnode++){
        int nindex = eachPoint * supportNodeSize + lnode;
        avgTemp += local_temperatures_array[nindex] * local_shapeFun_phis[nindex];
  }
  int material_id = material_ids[eachPoint] - 1;
  T myDensity = materials.density[material_id];

  int counterc =  materials._capacity_counters[material_id];
  int initc =  materials._capacity_inits[material_id];

  T myCapacity = d_interpolate1D(avgTemp,
                                 materials._capacity_temps,
                                 materials._capacity_values,
                                 initc,
                                 counterc);

  T myJacobian = local_jacobian_array[eachPoint];

  T avgFactor = myDensity * myCapacity * local_weight_array[eachPoint] * myJacobian;
  for(int i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
     for(int j = 0; j < supportNodeSize; j++){
        int out_id = (eachPoint * supportNodeSize * supportNodeSize) + (i * supportNodeSize) + j;
        T value = avgFactor * local_shapeFun_phis[eachPoint * supportNodeSize + i] * local_shapeFun_phis[eachPoint * supportNodeSize + j];
        local_capacity_matrices_array[out_id] = value;
     }
  }

}
//
template <typename T>
bool gpu_computeSOACapacityMatrix(T *local_capacity_matrices_array,
                                  T *local_temperatures_array,
                                  T *local_weight_array,
                                  T *local_jacobian_array,
                                  T *local_shapeFun_phis,
                                  int *material_ids,
                                  MaterialTable materials,
                                  int numPoints,
                                  int supportNodeSize,
                                  int threads_per_block,
                                  cudaStream_t stream)
{

  int size = numPoints;
  dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
  dim3 blockDim = dim3(threads_per_block,1,1);

  k_computeSOACapacityMatrix<<<gridDim, blockDim, 0, stream>>>(local_capacity_matrices_array,
                                                               local_temperatures_array,
                                                               local_weight_array,
                                                               local_jacobian_array,
                                                               local_shapeFun_phis,
                                                               material_ids,
                                                               materials,
                                                               numPoints,
                                                               supportNodeSize);

  cudaError_t cudaError = cudaGetLastError();
  if (cudaError != cudaSuccess)
  {
    std::cout << "GPU.CudaError k_computeSOACapacityMatrix"
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
//
template <typename T>
__global__ void k_computeSOAConductivityMatrix(T *local_conductivity_matrices_array,
                                              T *local_temperatures_array,
                                              T *local_weight_array,
                                              T *local_jacobian_array,
                                              T *local_shapeFun_phis,
                                              T *local_shapeFun_phis_dim,
                                              int *material_ids,
                                              MaterialTable materials,
                                              int numPoints,
                                              int supportNodeSize)
{
  int eachPoint = threadIdx.x + blockIdx.x * blockDim.x;
  double avgTemp  = 0.0;
  for(int lnode = 0; lnode < supportNodeSize; lnode++){
      int nindex = eachPoint * supportNodeSize + lnode;
       avgTemp += local_temperatures_array[nindex] * local_shapeFun_phis[nindex];
  }
  int material_id = material_ids[eachPoint] - 1;
  int countk = materials._kappa_counters[material_id];
  int initk = materials._kappa_inits[material_id];

  T myKappa = d_interpolate1D(avgTemp,
                              materials._kappa_temps,
                              materials._kappa_values,
                              0,//initk,
                              8);//countk);

  T avgFactor = myKappa * local_weight_array[eachPoint] * local_jacobian_array[eachPoint];
  int index_p =  eachPoint * supportNodeSize * supportNodeSize;
  for(int i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
    int index_row = index_p + i * supportNodeSize;
    for(int j = 0; j < supportNodeSize; j++){
        int out_id = index_row + j;
        local_conductivity_matrices_array[out_id] = local_shapeFun_phis_dim[out_id] * avgFactor;
     }
  }

}
//
template <typename T>
bool gpu_computeSOAConductivityMatrix(T *local_conductivity_matrices_array,
                                      T *local_temperatures_array,
                                      T *local_weight_array,
                                      T *local_jacobian_array,
                                      T *local_shapeFun_phis,
                                      T *local_shapeFun_phis_dim,
                                      int *material_ids,
                                      MaterialTable materials,
                                      int numPoints,
                                      int supportNodeSize,
                                      int threads_per_block,
                                      cudaStream_t stream)
{
  int size = numPoints;
  dim3 gridDim = dim3((size+threads_per_block-1) / threads_per_block, 1, 1);
  dim3 blockDim = dim3(threads_per_block,1,1);

  k_computeSOAConductivityMatrix<<<gridDim, blockDim, 0, stream>>>(local_conductivity_matrices_array,
                                                                   local_temperatures_array,
                                                                   local_weight_array,
                                                                   local_jacobian_array,
                                                                   local_shapeFun_phis,
                                                                   local_shapeFun_phis_dim,
                                                                   material_ids,
                                                                   materials,
                                                                   numPoints,
                                                                   supportNodeSize);

  cudaError_t cudaError = cudaGetLastError();
  if (cudaError != cudaSuccess)
  {
    std::cout << "GPU.CudaError k_computeSOAConductivityMatrix"
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
//

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
/*template __device__ float getMaterialKappa <float>(MaterialTable *materials,
                                                  int material_id,
                                                  float average_temperature);
//
template __device__ double getMaterialKappa <double>(MaterialTable *materials,
                                                    int material_id,
                                                    double average_temperature);
//
template __device__ float getMaterialDensity<float>(MaterialTable *materials,
                                                    int material_id);
//
template __device__ double getMaterialDensity<double>(MaterialTable *materials,
                                                      int material_id);
//
template __device__ float getMaterialCapacity<float>(MaterialTable *materials,
                                                    int material_id,
                                                    float average_temperature);
//
template __device__ double getMaterialCapacity<double>(MaterialTable *materials,
                                                      int material_id,
                                                      double average_temperature);*/
//
template __device__ float d_interpolate1D <float>(float query_value,
                                                float *reference_values,
                                                float *sought_values,
                                                int init_position,
                                                int counter);
//
//
template __device__ double d_interpolate1D <double>(double query_value,
                                                  double *reference_values,
                                                  double *sought_values,
                                                  int init_position,
                                                  int counter);
//
template bool gpu_assemble_global_matrix <float>(float *,
                                                uint*,
                                                uint*,
                                                float *,
                                                int,
                                                int,
                                                int,
                                                int,
                                                cudaStream_t stream);
template bool gpu_assemble_global_matrix <double> (double *,
                                                  uint*,
                                                  uint*,
                                                  double *,
                                                  int,
                                                  int,
                                                  int,
                                                  int,
                                                  cudaStream_t stream);
//
template bool gpu_computeSOACapacityMatrix<float>(float *local_capacity_matrices_array,
                                                  float *local_temperatures_array,
                                                  float *local_weight_array,
                                                  float *local_jacobian_array,
                                                  float *local_shapeFun_phis,
                                                  int *material_ids,
                                                  MaterialTable materials,
                                                  int numPoints,
                                                  int supportNodeSize,
                                                  int threads_per_block,
                                                  cudaStream_t stream);
//
template bool gpu_computeSOACapacityMatrix<double>(double *local_capacity_matrices_array,
                                                  double *local_temperatures_array,
                                                  double *local_weight_array,
                                                  double *local_jacobian_array,
                                                  double *local_shapeFun_phis,
                                                  int *material_ids,
                                                  MaterialTable materials,
                                                  int numPoints,
                                                  int supportNodeSize,
                                                  int threads_per_block,
                                                  cudaStream_t stream);
//
template bool gpu_computeSOAConductivityMatrix<float>(float *local_conductivity_matrices_array,
                                                      float *local_temperatures_array,
                                                      float *local_weight_array,
                                                      float *local_jacobian_array,
                                                      float *local_shapeFun_phis,
                                                      float *local_shapeFun_phis_dim,
                                                      int *material_ids,
                                                      MaterialTable materials,
                                                      int numPoints,
                                                      int supportNodeSize,
                                                      int threads_per_block,
                                                      cudaStream_t stream);
//
template bool gpu_computeSOAConductivityMatrix<double>(double *local_conductivity_matrices_array,
                                                       double *local_temperatures_array,
                                                       double *local_weight_array,
                                                       double *local_jacobian_array,
                                                       double *local_shapeFun_phis,
                                                       double *local_shapeFun_phis_dim,
                                                       int *material_ids,
                                                       MaterialTable materials,
                                                       int numPoints,
                                                       int supportNodeSize,
                                                       int threads_per_block,
                                                       cudaStream_t stream);
//
#endif //ASSEMBLY_KERNELS_H
