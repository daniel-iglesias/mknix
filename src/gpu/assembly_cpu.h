#ifndef ASSEMBLY_CPU_H
#define ASSEMBLY_CPU_H

#include "parallel_helper.h"

//namespace AssemblyCPU
//{


  template <typename T>
  bool cpu_assemble_global_matrix(T* global_matrix,
                                  int* full_map,
                                  int num_cells,
                                  int support_node_size,
                                  int number_elements);

#endif //ASSEMBLY_CPU_H
