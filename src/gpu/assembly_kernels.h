#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H

#include <cuda.h>
#include <cuda_runtime_api.h>

//#include "cuda_helper.h"

//namespace AssemblyKernels
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

  bool map_global_matrix(std::vector<int> &full_map,
                         std::vector<int> &vec_ind,
                         std::vector<int> &cvec_ptr,
                         int *presence_matrix,
                         int number_rows,
                         int number_columns);

  bool build_CCS_sparse_matrix_from_map(std::vector<int> &full_map,
                                        std::vector<int> &col_ind,
                                        std::vector<int> &row_ptr,
                                        int *presence_matrix,
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
