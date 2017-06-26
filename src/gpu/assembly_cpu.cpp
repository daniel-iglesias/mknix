//assembly_cpu.cpp
#include <vector>
#include "assembly_cpu.h"
#include "functions_cpu.h"
#include "LMX/lmx_mat_matrix.h"
#include "LMX/lmx_mat_type_csc.h"
#include "LMX/lmx_mat_type_gmm_sparse1.h"

/*void computeSOACapacityFactor(double *local_capacity_factor_array,
                              double *local_temperatures_array,
                              double *local_weight_array,
                              double *local_jacobian_array,
                              double *local_shapeFun_phis,
                              int *material_ids,
                              MaterialTable *materials,
                              int numPoints,
                              int supportNodeSize,
                              int tid)
{
  for(int eachPoint = 0; eachPoint < numPoints; eachPoint++){
    double avgTemp  = 0.0;
      for(int lnode = 0; lnode < supportNodeSize; lnode++){
          int nindex = eachPoint * supportNodeSize + lnode;
           avgTemp += local_temperatures_array[nindex] * local_shapeFun_phis[nindex];
      }
      int material_id = material_ids[eachPoint] - 1;
      double myDensity = getMaterialDensity (materials,
                                             material_id);;
      double myCapacity = getMaterialCapacity(materials,
                                         material_id,
                                         avgTemp);
      //double avgFactor = myDensity * myCapacity * local_weight_array[eachPoint] * local_jacobian_array[eachPoint];
      double avgFactor =  local_jacobian_array[eachPoint];
      local_capacity_factor_array[eachPoint] = avgFactor;
  }
}*/

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
                              int max_threads)//allows multithreaded
{
  int points_thread = (numPoints + max_threads - 1)/max_threads;//points per CPU thread
  int start_Point = tid * points_thread;
  int end_Point = ((tid+1)*points_thread < numPoints)? (tid+1)*points_thread:numPoints;
  for(int eachPoint = start_Point; eachPoint < end_Point; eachPoint++){
    T avgTemp  = 0.0;
    for(int lnode = 0; lnode < supportNodeSize; lnode++){
        int nindex = eachPoint * supportNodeSize + lnode;
         avgTemp += local_temperatures_array[nindex] * local_shapeFun_phis[nindex];
    }
    int material_id = material_ids[eachPoint] - 1;
    T myDensity = getMaterialDensity (materials,
                                      material_id);;
    T myCapacity = getMaterialCapacity(materials,
                                       material_id,
                                       avgTemp);
    T myJacobian = local_jacobian_array[eachPoint];
    //if(myJacobian )
    T avgFactor = myDensity * myCapacity * local_weight_array[eachPoint] * myJacobian;
    for(int i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
       for(int j = 0; j < supportNodeSize; j++){
          int out_id = (eachPoint * supportNodeSize * supportNodeSize) + (i * supportNodeSize) + j;
          T value = avgFactor * local_shapeFun_phis[eachPoint * supportNodeSize + i] * local_shapeFun_phis[eachPoint * supportNodeSize + j];
          local_capacity_matrices_array[out_id] = value;
       }
     }
   }
}
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
                                  int tid)
{
  for(int eachPoint = 0; eachPoint < numPoints; eachPoint++){
    double avgTemp  = 0.0;
    for(int lnode = 0; lnode < supportNodeSize; lnode++){
        int nindex = eachPoint * supportNodeSize + lnode;
         avgTemp += local_temperatures_array[nindex] * local_shapeFun_phis[nindex];
    }
    int material_id = material_ids[eachPoint] - 1;
    double myKappa = getMaterialKappa(materials,
                                    material_id,
                                    avgTemp);
    double avgFactor = myKappa * local_weight_array[eachPoint] * local_jacobian_array[eachPoint];
    int index_p =  eachPoint * supportNodeSize * supportNodeSize;

    local_conductivity_factor_array[eachPoint] =  avgFactor;
   }
}*/
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
                                  int max_threads)//allows multithreaded
{
//  std::cout << "Inside computeSOAConductivityMatrix" << std::endl;
//  std::cout << "Materials pointer = " << materials << std::endl;
  //many ways to divide this,
//  std::cout << "thread "<< tid << std::endl;
  int points_thread = (numPoints + max_threads - 1)/max_threads;//points per CPU thread
//  std::cout << "points_thread "<< points_thread << std::endl;
  int start_Point = tid * points_thread;
//  std::cout << "start_Point "<< start_Point << std::endl;
  int end_Point = ((tid+1)*points_thread < numPoints)? (tid+1)*points_thread:numPoints;
//  std::cout << "end_Point "<< end_Point << std::endl;
  for(int eachPoint = start_Point; eachPoint < end_Point; eachPoint++){
    //std::cout <<"point " << eachPoint << "; ";
    T avgTemp  = 0.0;
    for(int lnode = 0; lnode < supportNodeSize; lnode++){
        int nindex = eachPoint * supportNodeSize + lnode;
         avgTemp += local_temperatures_array[nindex] * local_shapeFun_phis[nindex];
    }
    int material_id = material_ids[eachPoint] - 1;
    //std::cout << "material_id = " << material_id;
    T myKappa = getMaterialKappa(materials,
                                 material_id,
                                 avgTemp);
    //std::cout << "myKappa = " << myKappa << "; ";
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
}
/**
 * Assembles a Global Matrix with multi-threaded safe atomics
 * @param  {[type]} T* array       the global matrix
 * @param  {[type]} int vector     The map for direct assembly
 * @param  {[type]} T* array       array of values to be assermbled
 * @param  {[type]} int size       Total number of Gauss points
 * @param  {[type]} int size       number of support nodes per gausspoint
 * @param  {[type]} int id         the thread ID
 * @param  {[type]} int size       number of threads launched
 */
template <typename T>
void atomicAssembleGlobalMatrix(std::atomic<T>* globalMatrix,
                                std::vector<uint> &fullMap,
                                std::vector<uint> &node_map,
                                T *local_matrices_array,
                                int numPoints,
                                int supportNodeSize,
                                int tid,
                                int max_threads,
                                bool isCSC)
{
  //many ways to divide this,
  int points_thread = (numPoints + max_threads - 1)/max_threads;//points per CPU thread
  int start_Point = tid * points_thread;
  int end_Point = ((tid+1)*points_thread < numPoints)? (tid+1)*points_thread:numPoints;
  for(int eachPoint = start_Point; eachPoint < end_Point; eachPoint++){
    for(int i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
       for(int j = 0; j < supportNodeSize; j++){
          int pos_id = eachPoint * supportNodeSize * supportNodeSize + i * supportNodeSize + j;
          T value = local_matrices_array[pos_id];
          int node_pos = node_map[pos_id];
          int globalPos = fullMap[node_pos];
          atomic_fetch_add(&globalMatrix[globalPos], value);
       }
     }
   }
}
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
                          bool isCSC)
{
  if(isCSC){//CSC format TODO:retest and test again
    for(uint eachPoint = 0; eachPoint < numPoints; eachPoint++){
      for(uint i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
         for(uint j = 0; j < supportNodeSize; j++){
            int pos_id = (eachPoint * supportNodeSize * supportNodeSize) + (i * supportNodeSize) + j;
            T value = local_matrices_array[pos_id];
            int node_pos = node_map[pos_id];
            int globalPos = fullMap[node_pos];
            globalMatrix[globalPos] += value;
         }
       }
     }
  }else{//CSR format
    for(uint eachPoint = 0; eachPoint < numPoints; eachPoint++){
      for(uint i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
         for(uint j = 0; j < supportNodeSize; j++){
            int pos_id = eachPoint * supportNodeSize * supportNodeSize + i * supportNodeSize + j;
            T value = local_matrices_array[pos_id];
            int node_pos = node_map[pos_id];
            int globalPos = fullMap[node_pos];
            globalMatrix[globalPos] += value;
         }
       }
     }
  }

}


template <typename T>
void AssembleGlobalMatrix(T* globalMatrix,
                          std::vector<uint> &fullMap,
                          std::vector<uint> &node_map,
                          T *local_matrices_array,
                          int numPoints,
                          int supportNodeSize,
                          bool isCSC)
{
  if(isCSC){//CSC format TODO:retest and test again
    for(uint eachPoint = 0; eachPoint < numPoints; eachPoint++){
      for(uint i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
         for(uint j = 0; j < supportNodeSize; j++){
            int pos_id = (eachPoint * supportNodeSize * supportNodeSize) + (i * supportNodeSize) + j;
            T value = local_matrices_array[pos_id];
            int node_pos = node_map[pos_id];
            int globalPos = fullMap[node_pos];
            globalMatrix[globalPos] += value;
         }
       }
     }
  }else{//CSR format
    for(uint eachPoint = 0; eachPoint < numPoints; eachPoint++){
      for(uint i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
         for(uint j = 0; j < supportNodeSize; j++){
            int pos_id = eachPoint * supportNodeSize * supportNodeSize + i * supportNodeSize + j;
            T value = local_matrices_array[pos_id];
            int node_pos = node_map[pos_id];
            int globalPos = fullMap[node_pos];
            globalMatrix[globalPos] += value;
         }
       }
     }
  }

}
/**
 * Thread Wrapper to launch compute in multi threaded version.
 * @param  {[type]} void* pointer   pointer to options structure
 */
void* threadWrapper(void* ptr){
  p_struct *parameters;
  parameters = (p_struct*) ptr;
  atomicAssembleGlobalMatrix(parameters->globalMatrix,
                             *parameters->fullMap,
                             *parameters->nodesMap,
                             parameters->local_matrices_array,
			                       parameters->numCells,
                             parameters->supportNodeSize,
                             parameters->thread_id,
                             parameters->max_threads,
                             parameters->use_csc);
}

/**
 * Thread Wrapper to launch compute in multi threaded version.
 * @param  {[type]} void* pointer   pointer to options structure
 */
void* computeCapacityThreadWrapper(void* ptr)
{
  p_calc_SOA_cap_struct *parameters;
  parameters = (p_calc_SOA_cap_struct*) ptr;
  computeSOACapacityMatrix(parameters->local_capacity_array,
                           parameters->local_temperatures_cap_array,
                           parameters->local_weight_cap_array,
                           parameters->local_jacobian_cap_array,
                           parameters->local_shapes_phis_array,
                           parameters->material_ids,
                           parameters->materials,
                           parameters->number_points,
                           parameters->supportNodeSize,
                           parameters->thread_id,
                           parameters->max_threads);
}

/**
 * Thread Wrapper to launch compute in multi threaded version.
 * @param  {[type]} void* pointer   pointer to options structure
 */
void* computeConductivityThreadWrapper(void* ptr)
{
  //std::cout << "Inside computeConductivityThreadWrapper" << std::endl;
  p_calc_SOA_cond_struct *parameters;
  parameters = (p_calc_SOA_cond_struct*) ptr;
  //std::cout << "Materials pointer = " << (uint)(parameters->materials) << std::endl;
  computeSOAConductivityMatrix(parameters->local_conductivity_array,
                               parameters->local_temperatures_cond_array,
                               parameters->local_weight_cond_array,
                               parameters->local_jacobian_cond_array,
                               parameters->local_shapes_phis_array,
                               parameters->local_shapes_phis_dim_array,
                               parameters->material_ids,
                               parameters->materials,
                               parameters->number_points,
                               parameters->supportNodeSize,
                               parameters->thread_id,
                               parameters->max_threads);
}

/**
 * Cast a directly assembled matrix into GMM compatible sparse matrix
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void cast_into_gmm_csc_type(gmm::csc_matrix<T>& gmm_matrix,
                            std::vector<T> &values_array,
                            std::vector<uint> &vec_ind,
                            std::vector<uint> &cvec_ptr,
                            int number_rows,
                            int number_columns)
{
  gmm_matrix.pr = values_array.data();
  gmm_matrix.ir = vec_ind.data();
  gmm_matrix.jc = cvec_ptr.data();
  gmm_matrix.nr = number_rows;
  gmm_matrix.nc = number_columns;
}

/*template <typename T>
void cast_into_eigen_type(SparseMatrix<T> &eigen_ref,
                        std::vector<T> &values_array,
                        std::vector<uint> &vec_ind,
                        std::vector<uint> &cvec_ptr,
                        int number_rows,
                        int number_columns,
                        bool use_csc)
  {
    T* values = eigen_ref.valuePtr();
    for(int i = 0; i < values_array.size(); i++){
      values[i] = values_array[i];
    }
    int* indices = eigen_ref.innerIndexPtr();
    for(int i = 0; i < vec_ind.size(); i++){
      indices[i] = vec_ind[i];
    }
    int* vec_ptr = eigen_ref.outerIndexPtr();
    for(int i = 0; i < cvec_ptr.size(); i++){
      vec_ptr[i] = cvec_ptr[i];
    }
  }*/
//
template <typename T>
void reserve_eigen_type(SparseMatrix<T> &eigen_ref,
                        int number_elements)
{
  eigen_ref.reserve(number_elements);
  //eigen_ref.
}
//

/**
 * Cast a directly assembled matrix into GMM compatible sparse matrix
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void cast_into_lmx_type(lmx::Matrix<T> &lmx_ref,
                        std::vector<T> &values_array,
                        std::vector<uint> &vec_ind,
                        std::vector<uint> &cvec_ptr,
                        int number_rows,
                        int number_columns,
                        bool use_csc)//if false will use csr instead
{
 if(use_csc){
   cast_into_lmx_csc_type(lmx_ref,
                          values_array,
                          vec_ind,
                          cvec_ptr,
                          number_rows,
                          number_columns);
 }else{
   cast_into_lmx_csr_type(lmx_ref,
                          values_array,
                          vec_ind,
                          cvec_ptr,
                          number_rows,
                          number_columns);
 }
}

/**
 * Cast a directly assembled matrix into GMM compatible sparse matrix
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void cast_into_eigen_type(SparseMatrix<T> &eigen_ref,
                        std::vector<T> &values_array,
                        std::vector<uint> &vec_ind,
                        std::vector<uint> &cvec_ptr,
                        int number_rows,
                        int number_columns)//if false will use csr instead
{
  T* vals = eigen_ref.valuePtr();
  int* inner = eigen_ref.innerIndexPtr();
  int* outer = eigen_ref.outerIndexPtr();

  for(int i = 0; i < values_array.size(); i++){ vals[i] = values_array[i];}
  for(int i = 0; i < vec_ind.size(); i++){ vals[i] = vec_ind[i];}
  for(int i = 0; i < cvec_ptr.size(); i++){ vals[i] = cvec_ptr[i];}

  /*memcpy( vals, values_array.data(),  values_array.size() * sizeof(T));
  memcpy( inner, vec_ind.data(),  vec_ind.size() * sizeof(uint));
  memcpy( outer, cvec_ptr.data(),  cvec_ptr.size() * sizeof(uint));*/

}

/**
 * Cast a directly assembled matrix into GMM compatible sparse matrix
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void cast_into_lmx_csc_type(lmx::Matrix<T> &lmx_ref,
                            std::vector<T> &values_array,
                            std::vector<uint> &vec_ind,
                            std::vector<uint> &cvec_ptr,
                            int number_rows,
                            int number_columns)
{
  //std::cout << "inside cast_into_lmx_csc_type" << std::endl;
  gmm::csc_matrix<T> gmm_matrix;
  gmm_matrix.pr = values_array.data();
  gmm_matrix.ir = vec_ind.data();
  gmm_matrix.jc = cvec_ptr.data();
  gmm_matrix.nr = number_rows;
  gmm_matrix.nc = number_columns;
  //std::cout << "about to use gmm_csc_cast" << std::endl;
  lmx_ref.gmm_csc_cast(gmm_matrix);
  //std::cout << "after gmm_csc_cast" << std::endl;
  gmm_matrix.pr  = NULL;
  gmm_matrix.ir = NULL;
  gmm_matrix.jc = NULL;
}
/**
 * Cast a directly assembled matrix into GMM compatible sparse matrix
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void cast_into_lmx_csr_type(lmx::Matrix<T> &lmx_ref,
                            std::vector<T> &values_array,
                            std::vector<uint> &vec_ind,
                            std::vector<uint> &cvec_ptr,
                            int number_rows,
                            int number_columns)
{
  //std::cout << "inside cast_into_lmx_csr_type" << std::endl;
  gmm::csr_matrix<T> gmm_matrix;
  gmm_matrix.pr = values_array.data();
  gmm_matrix.ir = vec_ind.data();
  gmm_matrix.jc = cvec_ptr.data();
  gmm_matrix.nr = number_rows;
  gmm_matrix.nc = number_columns;
  //std::cout << "about to use gmm_csr_cast" << std::endl;
  lmx_ref.gmm_csr_cast(gmm_matrix);
  //std::cout << "after gmm_csr_cast" << std::endl;
  gmm_matrix.pr  = NULL;
  gmm_matrix.ir = NULL;
  gmm_matrix.jc = NULL;

}

/**
 * initializes and array to a value in single thread
 * @param  {[type]} T* array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void init_host_array_to_value(T *array,
                              T value,
                              int size)
{
  //std::cout<< "inside init_host_array_to_value" << std::endl;
  for(int i = 0; i < size; i++) array[i] = value;
}

/**
 * initializes and array to a value in single thread
 * @param  {[type]} std::vector<T> array        array
 * @param  {[type]} T  value        initialization value
 * @param  {[type]} int size        number of elements in the array
 */
template <typename T>
void init_host_array_to_value(std::vector<T> &array,
                              T value,
                              int size)
{
  //std::cout<< "inside init_host_array_to_value" << std::endl;
  array.assign(size,value);
  //for(int i = 0; i < size; i++) array[i] = value;
}

 /**
  * Checks that all values in are between given limits
  * @param  {[type]} T* array        array
  * @param  {[type]} T  upper_limit  upper limit value
  * @param  {[type]} T  lower_limit  lower limit value
  * @param  {[type]} int size        number of elements in the array
  * @param  {[type]} std::string name    Name of the array
  */
 template <typename T>
 void check_host_array_for_limits(T *array,
                                 T upper_limit,
                                 T lower_limit,
                                 int size,
                                 std::string array_name)
 {
   int num_errors = 0;
    for(int i = 0; i < size; i++) {
      if(array[i] < lower_limit) num_errors++;
      else if(array[i] > upper_limit) num_errors++;
    }
    if(num_errors == 0) std::cout << "Correct: Array " << array_name << " is correct between given limits " << lower_limit << " to " << upper_limit << std::endl;
    else std::cout << "Error: Array " << array_name << " has " << num_errors<< " outside given limits " << lower_limit << " to " << upper_limit << std::endl;
  }

/**
* @brief Creates a Map for a global sparse matrix in either the Compressed
* Column Storage format (CCS) or Compressed Row Storage(CRS) as explained
* in http://www.netlib.org/utk/people/JackDongarra/etemplates/node373.html
* This function also allocates all the supplementary arrays for the sparse format.
* the allocation of the values array is not performed here as depends the data types and wether atomics will be used.
* @param full_map Reference to the global matrix positions, accesible by [row_id * totcols + col_id ]. lhs.
* @param row_ind. Array with the vec_ind of each element. lhs.
* @param cvec_ptr. Array with the memory position where each cvec starts. lhs.
* @param presence_matrix Reference to the global matrix counter where if an element is nonzero will have a counter greater than zero. rhs.
* @param all_point_nodes Reference to the array where each point keeps the relation of support nodes is using. rhs.
* @param number_rows. Integer with the total number of rows of the matrix. rhs.
* @param number_columns. Integer with the total number of columns of the matrix. rhs.
* @return void. lhs
**/
bool map_global_matrix(std::vector<uint> &full_map,
                      std::vector<uint> &vec_ind,
                      std::vector<uint> &cvec_ptr,
                      int *presence_matrix,
                      int number_rows,
                      int number_columns,
                      bool isCSC)
{

if(isCSC)
build_CSC_sparse_matrix_from_map(full_map,
                                vec_ind,
                                cvec_ptr,
                                presence_matrix,
                                number_rows,
                                number_columns);
else
build_CSR_sparse_matrix_from_map(full_map,
                                vec_ind,
                                cvec_ptr,
                                presence_matrix,
                                number_rows,
                                number_columns);


return true;

}

/**
 * @brief Allocates and Creates a Map for a global sparse matrix in the Compressed
 * Row Storage format (CRS) as explained in http://www.netlib.org/utk/people/JackDongarra/etemplates/node373.html
 * This function also allocates all the supplementary arrays for the sdparse format.
 * the allocation of the values array is not performed here as depends the data types and wether atomics will be used.
 * @param full_map Reference to the global matrix positions, accesible by [row_id * totcols + col_id ]. lhs.
 * @param col_ind. Array with the col id of each element. lhs.
 * @param row_ptr. Array with the memory position where each row starts. lhs.
 * @param presence_matrix Reference to the global matrix counter where if an element is nonzero will have a counter greater than zero. rhs.
 * @param number_rows. Integer with the total number of rows of the matrix. rhs.
 * @param number_columns. Integer with the total number of columns of the matrix. rhs.
 * @return void. lhs
 **/
bool build_CSR_sparse_matrix_from_map(std::vector<uint> &full_map,
                                      std::vector<uint> &col_ind,
                                      std::vector<uint> &row_ptr,
                                      int *presence_matrix,
                                      int number_rows,
                                      int number_columns)//Compressed Row Storage
{
//calculates total memory needed to allocate memory
  uint total_elements = 0;
  for(uint i = 0; i < number_rows * number_columns; i++)
      if(presence_matrix[i] > 0) total_elements++;

  col_ind.resize(total_elements);
  full_map.resize(number_rows * number_columns);
  row_ptr.resize(number_rows + 1);

  //prefix sum scan will give us the full map
  full_map[0] = 0;
  for(uint row = 0; row < number_rows; row++){
     for(uint col = 0; col < number_columns; col++){
       if(col!=0 || row != 0){
         int index = row * number_columns + col;
         full_map[index] = full_map[index-1] + presence_matrix[index-1];//we need to avoid the first element
       }
      }
 }

  for(uint row = 0; row < number_rows; row++){
    for(uint col = 0; col < number_columns; col++){
      if(presence_matrix[row*number_columns + col] > 0){
        col_ind[full_map[row*number_columns + col]] = col;//col id of every element
      }

    }
    row_ptr[row] = full_map[row*number_columns];//pointers to start of every row
  }

  row_ptr[number_rows] = total_elements+1;//convention
//  std::cout << "Total Non Zero Elements = " << total_elements << std::endl;
  return true;
}

  /**
  * @brief Allocates and Creates a Map for a global sparse matrix in the Compressed
  * Column Storage format (CCS) as explained in http://www.netlib.org/utk/people/JackDongarra/etemplates/node373.html
  * This function also allocates all the supplementary arrays for the sparse format.
  * the allocation of the values array is not performed here as depends the data types and wether atomics will be used.
  * @param full_map Reference to the global matrix positions, accesible by [row_id * totcols + col_id ]. lhs.
  * @param row_ind. Array with the col id of each element. lhs.
  * @param col_ptr. Array with the memory position where each row starts. lhs.
  * @param presence_matrix Reference to the global matrix counter where if an element is nonzero will have a counter greater than zero. rhs.
  * @param number_rows. Integer with the total number of rows of the matrix. rhs.
  * @param number_columns. Integer with the total number of columns of the matrix. rhs.
  * @return void. lhs
  **/
    bool build_CSC_sparse_matrix_from_map(std::vector<uint> &full_map,
                                          std::vector<uint> &row_ind,
                                          std::vector<uint> &col_ptr,
                                          int *presence_matrix,
                                          int number_rows,
                                          int number_columns)//Compressed Column Storage
    {
      //calculates total memory needed to allocate memory
        uint total_elements = 0;
        for(uint i = 0; i < number_rows * number_columns; i++)
            if(presence_matrix[i] > 0) total_elements++;

        row_ind.resize(total_elements);
        full_map.resize(number_rows * number_columns);
        col_ptr.resize(number_columns + 1);

       std::vector<int> byCols_presence(number_rows * number_columns);
       for(uint row = 0; row < number_rows; row++){
          for(uint col = 0; col < number_columns; col++){
             byCols_presence[col*number_rows + row] = presence_matrix[row * number_columns + col];
          }
        }


        //prefix sum scan will give us the full map
        full_map[0] = 0;
        for(uint row = 0; row < number_rows; row++){
           for(uint col = 0; col < number_columns; col++){
             if(col + row != 0){
               int index = row * number_columns + col;
               full_map[index] = full_map[index-1] + byCols_presence[index-1];//we need to avoid the first element
               //full_map[index] = full_map[index-1] + presence_matrix[index_presence-1];//we need to avoid the first element
             }
            }
       }

          for(uint col = 0; col < number_columns; col++){
            for(uint row = 0; row < number_rows; row++){
                if(byCols_presence[col * number_rows + row] != 0){
                  row_ind[full_map[col * number_rows + row]] = row;//row id of every element
                }

          }
          col_ptr[col] = full_map[col * number_rows];//pointers to start of every col
        }
        col_ptr[number_columns] = total_elements+1;//convention
        //std::cout << "Total Non Zero Elements = " << total_elements << std::endl;

        return true;
    }


    ///////////////////////////////////////////////////////////////////////////
    //////////////// templates parts ///////////////////////////////////////////
    //
  /*  template void computeSOATemperatureAndFactors<float>(float *local_capacity_factor,//output
                                                         float *local_conductivity_factor,//output
                                                         float *local_temperatures_array,
                                                         float *local_shapeFun_phis,
                                                         float *jacobian_array,
                                                         float *weight_array,
                                                         int *material_ids,
                                                         MaterialTable *materials,
                                                         int numPoints,
                                                         int supportNodeSize);
    //
    template void computeSOATemperatureAndFactors<double>(double *local_capacity_factor,//output
                                                          double *local_conductivity_factor,//output
                                                          double *local_temperatures_array,
                                                          double *local_shapeFun_phis,
                                                          double *jacobian_array,
                                                          double *weight_array,
                                                          int *material_ids,
                                                          MaterialTable *materials,
                                                          int numPoints,
                                                          int supportNodeSize);*/
    //
    template
    void computeSOACapacityMatrix<float>(float *local_capacity_matrices_array,
                                        float *local_temperatures_array,
                                        float *local_weight_array,
                                        float *local_jacobian_array,
                                        float *local_shapeFun_phis,
                                        int *material_ids,
                                        MaterialTable *materials,
                                        int numPoints,
                                        int supportNodeSize,
                                        int tid,
                                        int max_threads);
    //
    template
    void computeSOACapacityMatrix<double>(double *local_capacity_matrices_array,
                                          double *local_temperatures_array,
                                          double *local_weight_array,
                                          double *local_jacobian_array,
                                          double *local_shapeFun_phis,
                                          int *material_ids,
                                          MaterialTable *materials,
                                          int numPoints,
                                          int supportNodeSize,
                                          int tid,
                                          int max_threads);
    //
    template
    void computeSOAConductivityMatrix<float>(float *local_conductivity_matrices_array,
                                              float *local_temperatures_array,
                                              float *local_weight_array,
                                              float *local_jacobian_array,
                                              float *local_shapeFun_phis,
                                              float *local_shapeFun_phis_dim,
                                              int *material_ids,
                                              MaterialTable *materials,
                                              int numPoints,
                                              int supportNodeSize,
                                              int tid,
                                              int max_threads);
    //
    template
    void computeSOAConductivityMatrix<double>(double *local_conductivity_matrices_array,
                                              double *local_temperatures_array,
                                              double *local_weight_array,
                                              double *local_jacobian_array,
                                              double *local_shapeFun_phis,
                                              double *local_shapeFun_phis_dim,
                                              int *material_ids,
                                              MaterialTable *materials,
                                              int numPoints,
                                              int supportNodeSize,
                                              int tid,
                                              int max_threads);
    //
    template void atomicAssembleGlobalMatrix<float>(std::atomic<float>* globalMatrix,
                                                     std::vector<uint> &fullMap,
                                                     std::vector<uint> &node_map,
                                                     float *local_matrices_array,
                                                     int numPoints,
                                                     int supportNodeSize,
                                                     int tid,
                                                     int max_threads,
                                                     bool isCSC);
    //
    template void atomicAssembleGlobalMatrix<double>(std::atomic<double>* globalMatrix,
                                                     std::vector<uint> &fullMap,
                                                     std::vector<uint> &node_map,
                                                     double *local_matrices_array,
                                                     int numPoints,
                                                     int supportNodeSize,
                                                     int tid,
                                                     int max_threads,
                                                     bool isCSC);
    //
    //
    template void AssembleGlobalMatrix<float>(float *globalMatrix,
                                              std::vector<uint> &fullMap,
                                              std::vector<uint> &node_map,
                                              float *local_matrices_array,
                                              int numPoints,
                                              int supportNodeSize,
                                              bool isCSC);
    //
    template void AssembleGlobalMatrix<double>(double* globalMatrix,
                                               std::vector<uint> &fullMap,
                                               std::vector<uint> &node_map,
                                               double *local_matrices_array,
                                               int numPoints,
                                               int supportNodeSize,
                                               bool isCSC);

    //
    template void AssembleGlobalMatrix<float>(std::vector<float> &globalMatrix,
                                              std::vector<uint> &fullMap,
                                              std::vector<uint> &node_map,
                                              float *local_matrices_array,
                                              int numPoints,
                                              int supportNodeSize,
                                              bool isCSC);
    //
    template void AssembleGlobalMatrix<double>(std::vector<double> &globalMatrix,
                                               std::vector<uint> &fullMap,
                                               std::vector<uint> &node_map,
                                               double *local_matrices_array,
                                               int numPoints,
                                               int supportNodeSize,
                                               bool isCSC);

    //
    //
    template void cast_into_eigen_type<float>(SparseMatrix<float> &eigen_ref,
                                              std::vector<float> &values_array,
                                              std::vector<uint> &vec_ind,
                                              std::vector<uint> &cvec_ptr,
                                              int number_rows,
                                              int number_columns);
    //
    template void cast_into_eigen_type<double>(SparseMatrix<double> &eigen_ref,
                                               std::vector<double> &values_array,
                                               std::vector<uint> &vec_ind,
                                               std::vector<uint> &cvec_ptr,
                                               int number_rows,
                                               int number_columns);
    //
    template void reserve_eigen_type<double>(SparseMatrix<double> &eigen_ref,
                                             int number_elements);
    //
    template void reserve_eigen_type<float>(SparseMatrix<float> &eigen_ref,
                                            int number_elements);

    //
    template void cast_into_lmx_type<float>(lmx::Matrix<float> &lmx_ref,
                                            std::vector<float> &values_array,
                                            std::vector<uint> &vec_ind,
                                            std::vector<uint> &cvec_ptr,
                                            int number_rows,
                                            int number_columns,
                                            bool use_csc);
    //
    template void cast_into_lmx_type<double>(lmx::Matrix<double> &lmx_ref,
                                            std::vector<double> &values_array,
                                            std::vector<uint> &vec_ind,
                                            std::vector<uint> &cvec_ptr,
                                            int number_rows,
                                            int number_columns,
                                            bool use_csc);
    //
    template void cast_into_gmm_csc_type<float>(gmm::csc_matrix<float>& gmm_matrix,
                                                std::vector<float> &values_array,
                                                std::vector<uint> &vec_ind,
                                                std::vector<uint> &cvec_ptr,
                                                int number_rows,
                                                int number_columns);
    //
    template void cast_into_gmm_csc_type<double>(gmm::csc_matrix<double>& gmm_matrix,
                                                 std::vector<double> &values_array,
                                                 std::vector<uint> &vec_ind,
                                                 std::vector<uint> &cvec_ptr,
                                                 int number_rows,
                                                 int number_columns);
    //
    template void cast_into_lmx_csc_type<float>(lmx::Matrix<float> &lmx_ref,
                                                std::vector<float> &values_array,
                                                std::vector<uint> &vec_ind,
                                                std::vector<uint> &cvec_ptr,
                                                int number_rows,
                                                int number_columns);
    //
    template void cast_into_lmx_csc_type<double>(lmx::Matrix<double> &lmx_ref,
                                                std::vector<double> &values_array,
                                                std::vector<uint> &vec_ind,
                                                std::vector<uint> &cvec_ptr,
                                                int number_rows,
                                                int number_columns);
    //
    template void cast_into_lmx_csr_type<float>(lmx::Matrix<float> &lmx_ref,
                                                std::vector<float> &values_array,
                                                std::vector<uint> &vec_ind,
                                                std::vector<uint> &cvec_ptr,
                                                int number_rows,
                                                int number_columns);
    //
    template void cast_into_lmx_csr_type<double>(lmx::Matrix<double> &lmx_ref,
                                                std::vector<double> &values_array,
                                                std::vector<uint> &vec_ind,
                                                std::vector<uint> &cvec_ptr,
                                                int number_rows,
                                                int number_columns);
    //
    template void init_host_array_to_value<double>(double *array,
                                                   double value,
                                                   int size);
    template void init_host_array_to_value<float>(float *array,
                                                  float value,
                                                  int size);
    template void init_host_array_to_value<int>(int *array,
                                                int value,
                                                int size);
 //
 template void init_host_array_to_value<double>(std::vector<double> &array,
                                                double value,
                                                int size);
 template void init_host_array_to_value<float>(std::vector<float> &array,
                                               float value,
                                               int size);
 template void init_host_array_to_value<int>(std::vector<int> &array,
                                             int value,
                                             int size);
//
 template void check_host_array_for_limits<double>(double *array,
                                                   double upper_limit,
                                                   double lower_limit,
                                                   int size,
                                                   std::string array_name);
  //
  template void check_host_array_for_limits<float>(float *array,
                                                   float upper_limit,
                                                   float lower_limit,
                                                   int size,
                                                   std::string array_name);
//
template void check_host_array_for_limits<int>(int *array,
                                               int upper_limit,
                                               int lower_limit,
                                               int size,
                                               std::string array_name);
