#ifndef ASSEMBLY_CPU_H
#define ASSEMBLY_CPU_H

//#include "parallel_helper.h"//for now commented out, will include multi-cpu code here

//namespace AssemblyCPU
//{
template <typename T>
void init_host_array_to_value(T *array,
                              T value,
                              int size);


  template <typename T>
  bool cpu_assemble_global_matrix(T* global_matrix,
                                  int* full_map,
                                  int num_cells,
                                  int support_node_size,
                                  int number_elements);

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


#endif //ASSEMBLY_CPU_H
