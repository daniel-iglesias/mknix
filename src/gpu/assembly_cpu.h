#ifndef ASSEMBLY_CPU_H
#define ASSEMBLY_CPU_H
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "LMX/lmx.h"
#include "gmm/gmm_matrix.h"
#include "functions_cpu.h"
#include <cpu_run_type.h>
#include "../Eigen/Sparse"

#include <atomic>
#include <pthread.h>

/**
 * Structure containing all parameters for multithreaded launch of Capacity with SOA structures assembly

 */
struct p_calc_SOA_cap_struct{
  double *local_capacity_array;
  double *local_temperatures_cap_array;
  double *local_weight_cap_array;
  double *local_jacobian_cap_array;
  double *local_shapes_phis_array;
  int *material_ids;
  MaterialTable *materials;
  int number_points;
  int supportNodeSize;
  int thread_id;
  int max_threads;
};


/**
 * Structure containing all parameters for multithreaded launch of Conductivity with SOA structures assembly

 */
struct p_calc_SOA_cond_struct{
  double *local_conductivity_array;
  double *local_temperatures_cond_array;
  double *local_weight_cond_array;
  double *local_jacobian_cond_array;
  double *local_shapes_phis_array;
  double *local_shapes_phis_dim_array;
  int *material_ids;
  MaterialTable *materials;
  int number_points;
  int supportNodeSize;
  int thread_id;
  int max_threads;
};

/**
 * Structure containing all parameters for multithreaded launch of global matrix assembly
 * @param  {[type]} double* atomic array      Array of the local capacity factors
 * @param  {[type]} int vector*               Map of the global matrix
 * @param  {[type]} double* array             array of local values to be assembled
 * @param  {[type]} int size                  Number of cells
 * @param  {[type]} int size                  Number of nodes per guasspoint
 * @param  {[type]} int id                    index of thread
 * @param  {[type]} int size                  Number of threads launched
 */
struct p_struct{
  std::atomic<double>* globalMatrix;
  std::vector<uint> *fullMap;
  std::vector<uint> *nodesMap;
  double *local_matrices_array;
  int numCells;
  int supportNodeSize;
  int thread_id;
  int max_threads;
  bool use_csc;
};

//
/*void computeSOACapacityFactor(double *local_capacity_factor_array,
                              double *local_temperatures_array,
                              double *local_weight_array,
                              double *local_jacobian_array,
                              double *local_shapeFun_phis,
                              int *material_ids,
                              MaterialTable *materials,
                              int numPoints,
                              int supportNodeSize,
                              int tid);*/
//
/**
 * Computes the local SoA array for the Capacity Matrix
 * @param  {[type]} T* array       Array of the local capacity matrices
 * @param  {[type]} T* array       Array of capacity factors
 * @param  {[type]} T* array       array oflocal shapefunctions phis
 * @param  {[type]} int size       Total number of Gauss points
 * @param  {[type]} int size       number of support nodes per gausspoint
 * @param  {[type]} int id         the thread ID
 */
 template <typename T>
 void computeSOACapacityMatrix(T *local_capacity_matrices_array,
                               T *local_temperatures_array,
                               T *local_weight_array,
                               T *local_jacobian_array,
                               T *local_shapeFun_phis,
                               int *material_ids,
                               MaterialTable *materials,
                               int numPoints,
                               int supportNodeSize,
                               int tid,
                               int max_threads);
//
/*void computeSOAConductivityFactor(double *local_conductivity_factor_array,
                                  double *local_temperatures_array,
                                  double *local_weight_array,
                                  double *local_jacobian_array,
                                  double *local_shapeFun_phis,
                                  double *local_shapeFun_phis_dim,
                                  int *material_ids,
                                  MaterialTable *materials,
                                  int numPoints,
                                  int supportNodeSize,
                                  int tid);*/
//
/**
 * Computes the local SoA array for the Conductivity Matrix
 * @param  {[type]} T* array       Array of the local conductivity matrices
 * @param  {[type]} T* array       Array of conductivity factors
 * @param  {[type]} T* array       array oflocal shapefunctions phis
 * @param  {[type]} int size       Total number of Gauss points
 * @param  {[type]} int size       number of support nodes per gausspoint
 * @param  {[type]} int id         the thread ID
 */
 template <typename T>
 void computeSOAConductivityMatrix(T *local_conductivity_matrices_array,
                                   T *local_temperatures_array,
                                   T *local_weight_array,
                                   T *local_jacobian_array,
                                   T *local_shapeFun_phis,
                                   T *local_shapeFun_phis_dim,
                                   int *material_ids,
                                   MaterialTable *materials,
                                   int numPoints,
                                   int supportNodeSize,
                                   int tid,
                                   int max_threads);
//
template <typename T>
void atomicAssembleGlobalMatrix(std::atomic<T>* globalMatrix,
                                std::vector<uint> &fullMap,
                                std::vector<uint> &node_map,
                                T *local_matrices_array,
                                int numPoints,
                                int supportNodeSize,
                                int tid,
                                int max_threads,
                                bool isCSC);
//
/**
 * Assembles a Global Matrix in single threaded CPU
 * @param  {[type]} T* array       the global matrix
 * @param  {[type]} int vector     The map for direct assembly
 * @param  {[type]} T* array       array of values to be assermbled
 * @param  {[type]} int size       Total number of Gauss points
 * @param  {[type]} int size       number of support nodes per gausspoint
 */
template <typename T>
void AssembleGlobalMatrix(std::vector<T> &globalMatrix,
                          std::vector<uint> &fullMap,
                          std::vector<uint> &node_map,
                          T *local_matrices_array,
                          int numPoints,
                          int supportNodeSize,
                          bool isCSC);
//
//
/**
 * Assembles a Global Matrix in single threaded CPU
 * @param  {[type]} T* array       the global matrix
 * @param  {[type]} int vector     The map for direct assembly
 * @param  {[type]} T* array       array of values to be assermbled
 * @param  {[type]} int size       Total number of Gauss points
 * @param  {[type]} int size       number of support nodes per gausspoint
 */
template <typename T>
void AssembleGlobalMatrix(T* globalMatrix,
                          std::vector<uint> &fullMap,
                          std::vector<uint> &node_map,
                          T *local_matrices_array,
                          int numPoints,
                          int supportNodeSize,
                          bool isCSC);
//pthreads wrappers for several functions
void* threadWrapper(void* ptr);
void* computeCapacityThreadWrapper(void* ptr);
void* computeConductivityThreadWrapper(void* ptr);

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
/**
 * Cast a directly assembled matrix into GMM compatible sparse matrix
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
 /*template <typename T>
 void cast_into_eigen_type(SparseMatrix<T> &eigen_ref,
                         std::vector<T> &values_array,
                         std::vector<uint> &vec_ind,
                         std::vector<uint> &cvec_ptr,
                         int number_rows,
                         int number_columns,
                         bool use_csc);*/
//
template <typename T>
void reserve_eigen_type(SparseMatrix<T> &eigen_ref,
                        int number_elements);
//

template <typename T>
void cast_into_lmx_type(lmx::Matrix<T> &lmx_ref,
                        std::vector<T> &values_array,
                        std::vector<uint> &vec_ind,
                        std::vector<uint> &cvec_ptr,
                        int number_rows,
                        int number_columns,
                        bool use_csc);
//
template <typename T>
void cast_into_gmm_csc_type(gmm::csc_matrix<T>& gmm_matrix,
                            std::vector<T> &values_array,
                            std::vector<uint> &vec_ind,
                            std::vector<uint> &cvec_ptr,
                            int number_rows,
                            int number_columns);

//
template <typename T>
void cast_into_eigen_type(SparseMatrix<T> &eigen_ref,
                        std::vector<T> &values_array,
                        std::vector<uint> &vec_ind,
                        std::vector<uint> &cvec_ptr,
                        int number_rows,
                        int number_columns);
//
template <typename T>
void cast_into_lmx_csc_type(lmx::Matrix<T> &lmx_ref,
                            std::vector<T> &values_array,
                            std::vector<uint> &vec_ind,
                            std::vector<uint> &cvec_ptr,
                            int number_rows,
                            int number_columns);
//

template <typename T>
void cast_into_lmx_csr_type(lmx::Matrix<T> &lmx_ref,
                            std::vector<T> &values_array,
                            std::vector<uint> &vec_ind,
                            std::vector<uint> &cvec_ptr,
                            int number_rows,
                            int number_columns);

template <typename T>
void init_host_array_to_value(T *array,
                              T value,
                              int size);

/**
* initializes and array to a value in single thread
* @param  {[type]} std::vector<T> array        array
* @param  {[type]} T  value        initialization value
* @param  {[type]} int size        number of elements in the array
*/
template <typename T>
void init_host_array_to_value(std::vector<T> &array,
                              T value,
                              int size);
//
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

  bool map_global_matrix(std::vector<uint> &full_map,
                         std::vector<uint> &vec_ind,
                         std::vector<uint> &cvec_ptr,
                         int *presence_matrix,
                         int number_rows,
                         int number_columns,
                         bool isCSC = false);//CCS or CSR

  bool build_CSC_sparse_matrix_from_map(std::vector<uint> &full_map,
                                        std::vector<uint> &col_ind,
                                        std::vector<uint> &row_ptr,
                                        int *presence_matrix,
                                        int number_rows,
                                        int number_columns);//Compressed Sparse Column Storage

  bool build_CSR_sparse_matrix_from_map(std::vector<uint> &full_map,
                                        std::vector<uint> &col_ind,
                                        std::vector<uint> &row_ptr,
                                        int *presence_matrix,
                                        int number_rows,
                                        int number_columns);//Compressed Sparse Row Storage


#endif //ASSEMBLY_CPU_H
