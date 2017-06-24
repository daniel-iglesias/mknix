#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

//#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cpu/structures.h"

//#include "cuda_helper.h"
__device__ double d_getMaterialKappa (MaterialTable *materials,
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
                                  int* full_map,
                                  T* local_matrices_array,
                                  int num_cells,
                                  int support_node_size,
                                  int number_elements,
                                  int threads_per_block,
                                  cudaStream_t stream = 0);
//
template <typename T>
bool gpu_assemble_global_vector(T *global_matrix,
                                int* vector_positions,
                                T* local_vector,
                                int num_points,
                                int support_node_size,
                                int number_points,
                                int threads_per_block,
                                cudaStream_t stream);
//

//}
#endif //ASSEMBLY_KERNELS_H
