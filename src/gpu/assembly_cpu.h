#ifndef ASSEMBLY_CPU_H
#define ASSEMBLY_CPU_H
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "LMX/lmx.h"
#include "gmm/gmm_matrix.h"

#include <atomic>
#include <pthread.h>
//struct for multithreaded launch
struct p_struct{
  std::atomic<double>* globalMatrix;
  std::vector<int> *fullMap;
  double *local_matrices_array;
  int numCells;
  int supportNodeSize;
  int thread_id;
  int max_threads;
};

void atomicAssembleGlobalMatrix(std::atomic<double>* globalMatrix,
                                std::vector<int> &fullMap,
                                double *local_matrices_array,
                                int numPoints,
                                int supportNodeSize,
                                int tid,
                                int max_threads);

void* threadWrapper(void* ptr);

double inline atomic_fetch_add(std::atomic<double>* target, double value){
  double expected = target->load();
  while(!atomic_compare_exchange_weak(target, &expected, expected + value));
  return expected;
}

float inline atomic_fetch_add(std::atomic<float>* target, float value){
  float expected = target->load();
  while(!atomic_compare_exchange_weak(target, &expected, expected + value));
  return expected;
}

//#include "parallel_helper.h"//for now commented out, will include multi-cpu code here

//namespace AssemblyCPU
//{
/*template <typename T>
void cast_into_lmx_type(lmx::Matrix<data_type>& lmx_matrix,
                        T *values_array,
                        std::vector<int> &vec_ind,
                        std::vector<int> &cvec_ptr,
                        int number_rows,
                        int number_columns);*/

template <typename T>
void cast_into_gmm_type(gmm::csr_matrix<T>& gmm_matrix,
                        T *values_array,
                        std::vector<int> &vec_ind,
                        std::vector<int> &cvec_ptr,
                        int number_rows,
                        int number_columns);

template <typename T>
void init_host_array_to_value(T *array,
                              T value,
                              int size);
template <typename T>
void check_host_array_for_limits(T *array,
                                 T upper_limit,
                                 T lower_limit,
                                 int size,
                                 std::string array_name = "");

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
