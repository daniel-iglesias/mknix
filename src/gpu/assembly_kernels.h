#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "cuda_helper.h"

//namespace AssemblyKernels
//{
  template <typename T>
  bool gpu_assemble_global_matrix(T* global_matrix,
                                  int* full_map,
                                  int num_cells,
                                  int support_node_size,
                                  int number_elements,
                                  int threads_per_block,
                                  cudaStream_t stream = 0);

  bool map_global_matrix(int *presence_matrix,
                         int *all_point_nodes,
                         int *position_map,
                         int *col_ind,
                         int *row_ptr);

  bool build_CRS_sparse_matrix_from_map(int *full_map,
                                        int *presence_matrix,
                                        int *col_ind,
                                        int *row_ptr,
                                        int number_rows,
                                        int number_columns);//Compressed Column Storage

  bool build_CRS_sparse_matrix_from_map(int *full_map,
                                        int *presence_matrix,
                                        int *col_ind,
                                        int *row_ptr,
                                        int number_rows,
                                        int number_columns);//Compressed Row Storage

//}
#endif //ASSEMBLY_KERNELS_H
