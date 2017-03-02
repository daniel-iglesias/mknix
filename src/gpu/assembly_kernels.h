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



//}
#endif //ASSEMBLY_KERNELS_H
