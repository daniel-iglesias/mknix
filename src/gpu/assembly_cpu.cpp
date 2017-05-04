//assembly_cpu.cpp
#include <vector>
#include "assembly_cpu.h"
#include "functions_cpu.h"
#include "LMX/lmx_mat_matrix.h"
#include "LMX/lmx_mat_type_csc.h"
#include "LMX/lmx_mat_type_gmm_sparse1.h"

template <typename T>
void computeSOATemperatureAndFactors(T *local_capacity_factor,//output
                                     T *local_conductivity_factor,//output
                                     T *local_temperatures_array,
                                     T *local_shapeFun_phis,
                                     T *jacobian_array,
                                     T *weight_array,
                                     int *material_ids,
                                     MaterialTable *materials,
                                     int numPoints,
                                     int supportNodeSize)
{

  double debug_temp = 0.0;//DEBUG LINE
  double debug_jacobian = 0.0;//DEBUG LINE
  double debug_weight = 0.0;//DEBUG LINE
  for(int eachPoint = 0; eachPoint < numPoints; eachPoint++){
    T avgTemp = 0;
    for(int i = 0; i < supportNodeSize; i++){
      //avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
      int lid = eachPoint * supportNodeSize + i;
      avgTemp += local_temperatures_array[lid] + local_shapeFun_phis[lid];
    }
    debug_temp += avgTemp;//DEBUG LINE
    int material_id = material_ids[eachPoint];
    T abs_jacobian = std::abs(jacobian_array[eachPoint]);
    debug_jacobian += abs_jacobian;//DEBUG LINE
    T weight = weight_array[eachPoint];
    debug_weight += weight;//DEBUG LINE
    T density_val = getMaterialDensity (materials,
                                            material_id);
    T cap_val = getMaterialCapacity(materials,
                                        material_id,
                                        avgTemp);
    T avgCapacityFactor = density_val * cap_val * weight_array[eachPoint] * weight * abs_jacobian;
    local_capacity_factor[eachPoint] = avgCapacityFactor;
    T kappa_val = getMaterialKappa (materials,
                                        material_id,
                                        avgTemp);
    T avgConductivityFactor = kappa_val * weight * abs_jacobian;
    local_conductivity_factor[eachPoint] = avgConductivityFactor;
  }
  std::cout << "DEBUG::computeSOATemperatureAndFactors, Sum of avrg temps = " <<  debug_temp << std::endl;  //DEBUG LINE
  std::cout << "DEBUG::computeSOATemperatureAndFactors, Sum of avrg jacobians = " <<  debug_jacobian << std::endl;  //DEBUG LINE
  std::cout << "DEBUG::computeSOATemperatureAndFactors, Sum of avrg weights = " <<  debug_weight << std::endl;  //DEBUG LINE
}

template <typename T>
void computeSOACapacityMatrix(T *local_capacity_matrices_array,
                              T *local_capacity_factor,
                              T *local_shapeFun_phis,
                              int numPoints,
                              int supportNodeSize,
                              int tid)
{
  for(int eachPoint = 0; eachPoint < numPoints; eachPoint++){
    T avgFactor  = local_capacity_factor[eachPoint];
    for(int i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
       for(int j = 0; j < supportNodeSize; j++){
          int out_id = eachPoint * supportNodeSize * supportNodeSize + i * supportNodeSize + j;
          T value = avgFactor * local_shapeFun_phis[eachPoint * supportNodeSize + i] * local_shapeFun_phis[eachPoint * supportNodeSize + j];
          local_capacity_matrices_array[out_id] = value;
       }
     }
   }
}
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
                                  T *local_conductivity_factor,
                                  T *local_shapeFun_phis,
                                  int numPoints,
                                  int supportNodeSize,
                                  int tid)
{
  for(int eachPoint = 0; eachPoint < numPoints; eachPoint++){
    double avgFactor  = local_conductivity_factor[eachPoint];
    for(int i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
       for(int j = 0; j < supportNodeSize; j++){
          int out_id = eachPoint * supportNodeSize * supportNodeSize + i * supportNodeSize + j;
          double value = avgFactor * local_shapeFun_phis[eachPoint * supportNodeSize + i] * local_shapeFun_phis[eachPoint * supportNodeSize + j];
          local_conductivity_matrices_array[out_id] = value;
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
                                T *local_matrices_array,
                                int numPoints,
                                int supportNodeSize,
                                int tid,
                                int max_threads)
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
          int globalPos = fullMap[pos_id];
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
                          T *local_matrices_array,
                          int numPoints,
                          int supportNodeSize)
{
  for(uint eachPoint = 0; eachPoint < numPoints; eachPoint++){
    for(uint i = 0; i < supportNodeSize; i++){//equivalent to obtaining thermalNumber
       for(uint j = 0; j < supportNodeSize; j++){
          int pos_id = eachPoint * supportNodeSize * supportNodeSize + i * supportNodeSize + j;
          T value = local_matrices_array[pos_id];
          int globalPos = fullMap[pos_id];
          globalMatrix[globalPos] += value;
       }
     }
   }
}
/**
 * Thread Wrapper to launch assembly in multi threaded version.
 * @param  {[type]} void* pointer   pointer to options structure
 */
void* threadWrapper(void* ptr){
  p_struct *parameters;
  parameters = (p_struct*) ptr;
  atomicAssembleGlobalMatrix(parameters->globalMatrix,
                             *parameters->fullMap,
                             parameters->local_matrices_array,
			                       parameters->numCells,
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
  std::cout << "inside cast_into_lmx_csc_type" << std::endl;
  gmm::csc_matrix<T> gmm_matrix;
  gmm_matrix.pr = values_array.data();
  //std::cout << "(uint*)vec_ind.data()" << std::endl;
  gmm_matrix.ir = vec_ind.data();
  //std::cout << "(uint*)cvec_ptr.data()" << std::endl;
  gmm_matrix.jc = cvec_ptr.data();
  gmm_matrix.nr = number_rows;
  gmm_matrix.nc = number_columns;
  std::cout << "about to use gmm_csc_cast" << std::endl;
  lmx_ref.gmm_csc_cast(gmm_matrix);
  std::cout << "after gmm_csc_cast" << std::endl;
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
void cast_into_gmm_csr_type(gmm::csr_matrix<T>& gmm_matrix,
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
  std::cout<< "inside init_host_array_to_value" << std::cout;
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
  std::cout<< "inside init_host_array_to_value" << std::cout;
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
  for(uint i = 1; i < number_rows * number_columns; i++)
      if(presence_matrix[i] > 0) total_elements++;

  col_ind.resize(total_elements);
  full_map.resize(number_rows * number_columns);
  row_ptr.resize(number_rows + 1);

  //prefix sum scan will give us the full map
  full_map[0] = 0;
  for(uint i = 1; i < number_rows * number_columns; i++)
          full_map[i] = full_map[i-1] + presence_matrix[i-1];

  for(uint i = 1; i < number_rows; i++){
    for(uint j = 1; j < number_columns; j++){
      if(presence_matrix[i*number_columns + j] > 0){
        col_ind[full_map[i*number_columns + j]] = j;//col id of every element
      }

    }
    row_ptr[i] = full_map[i*number_columns];//pointers to start of every row
  }

  row_ptr[number_rows] = total_elements;//convention
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
        for(uint i = 1; i < number_rows * number_columns; i++)
            if(presence_matrix[i] > 0) total_elements++;

        row_ind.resize(total_elements);
        full_map.resize(number_rows * number_columns);
        col_ptr.resize(number_rows + 1);

        //prefix sum scan will give us the full map
        full_map[0] = 0;
        for(uint i = 1; i < number_rows * number_columns; i++)
                full_map[i] = full_map[i-1] + presence_matrix[i-1];

        for(uint j = 1; j < number_rows; j++){
          for(uint i = 1; i < number_columns; i++){
            if(presence_matrix[j * number_rows + i] > 0){
              row_ind[full_map[j * number_rows + i]] = i;//row id of every element
            }

          }
          col_ptr[j] = full_map[j * number_rows];//pointers to start of every col
        }

        col_ptr[number_columns] = total_elements;//convention

        return true;
    }


    ///////////////////////////////////////////////////////////////////////////
    //////////////// templates parts ///////////////////////////////////////////
    //
    template void computeSOATemperatureAndFactors<float>(float *local_capacity_factor,//output
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
                                                          int supportNodeSize);
    //
    template void computeSOACapacityMatrix<float>(float *local_capacity_matrices_array,
                                                  float *local_capacity_factor,
                                                  float *local_shapeFun_phis,
                                                  int numPoints,
                                                  int supportNodeSize,
                                                  int tid);
    //
    template void computeSOACapacityMatrix<double>(double *local_capacity_matrices_array,
                                                   double *local_capacity_factor,
                                                   double *local_shapeFun_phis,
                                                   int numPoints,
                                                   int supportNodeSize,
                                                   int tid);
    //
    template void computeSOAConductivityMatrix<float>(float *local_conductivity_matrices_array,
                                                      float *local_conductivity_factor,
                                                      float *local_shapeFun_phis,
                                                      int numPoints,
                                                      int supportNodeSize,
                                                      int tid);
    //
    template void computeSOAConductivityMatrix<double>(double *local_conductivity_matrices_array,
                                                       double *local_conductivity_factor,
                                                       double *local_shapeFun_phis,
                                                       int numPoints,
                                                       int supportNodeSize,
                                                       int tid);
    //
    template void atomicAssembleGlobalMatrix<float>(std::atomic<float>* globalMatrix,
                                                     std::vector<uint> &fullMap,
                                                     float *local_matrices_array,
                                                     int numPoints,
                                                     int supportNodeSize,
                                                     int tid,
                                                     int max_threads);
    //
    template void atomicAssembleGlobalMatrix<double>(std::atomic<double>* globalMatrix,
                                                     std::vector<uint> &fullMap,
                                                     double *local_matrices_array,
                                                     int numPoints,
                                                     int supportNodeSize,
                                                     int tid,
                                                     int max_threads);
    //
    template void AssembleGlobalMatrix<float>(std::vector<float> &globalMatrix,
                                              std::vector<uint> &fullMap,
                                              float *local_matrices_array,
                                              int numPoints,
                                              int supportNodeSize);
    //
    template void AssembleGlobalMatrix<double>(std::vector<double> &globalMatrix,
                                               std::vector<uint> &fullMap,
                                               double *local_matrices_array,
                                               int numPoints,
                                               int supportNodeSize);
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
    template void cast_into_gmm_csr_type<float>(gmm::csr_matrix<float>& gmm_matrix,
                                                std::vector<float> &values_array,
                                                std::vector<uint> &vec_ind,
                                                std::vector<uint> &cvec_ptr,
                                                int number_rows,
                                                int number_columns);
    //
    template void cast_into_gmm_csr_type<double>(gmm::csr_matrix<double>& gmm_matrix,
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
