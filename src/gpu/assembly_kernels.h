#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

//#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cpu/structures.h"

//#include "cuda_helper.h"
__device__ double d_getMaterialKappa (int *materials_kappa_counters,
                                      int *materials_kappa_inits,
                                      double* materials_kappa_temps,
                                      double* materials_kappa_values,
                                      int material_id,
                                      double average_temperature);

__device__ double d_getMaterialDensity (MaterialTable *materials,
                                      int material_id);

__device__ double d_getMaterialCapacity (MaterialTable *materials,
                                      int material_id,
                                      double average_temperature);
//
template <typename T>
__device__ T d_interpolate1D(T query_value,
                           T *reference_values,
                           T *sought_values,
                           int init_position,
                           int counter);
//
//namespace AssemblyKernels
template <typename T>
__global__ void k_assemble_global_vector(T* global_vector,
                                        int* full_map,
                                        T* local_matrices_array,
                                        int num_points,
                                        int support_node_size,
                                        int number_elements);
//{
  template <typename T>
  bool gpu_assemble_global_matrix(T* global_matrix,
                                  uint* full_map,
                                  uint* node_map,
                                  T* local_matrices_array,
                                  int num_cells,
                                  int support_node_size,
                                  int number_elements,
                                  int threads_per_block,
                                  cudaStream_t stream = 0);
//
template <typename T>
bool gpu_assemble_global_vector(T *global_matrix,
                                uint* vector_positions,
                                T* local_vector,
                                int num_points,
                                int support_node_size,
                                int number_points,
                                int threads_per_block,
                                cudaStream_t stream);
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
                                  cudaStream_t stream);
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
                                      cudaStream_t stream);
//}
#endif //ASSEMBLY_KERNELS_H
