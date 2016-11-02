//assembly_cpu.cpp


/**
 * Makes a direct assembly of the global matrix from map in single thread
 * @param  {[type]} T *global_matrix        global sparse matrix
 * @param  {[type]} int *full_map           Array mapping positions in global matrix
 * @param  {[type]} int num_points          number of gauss points
 * @param  {[type]} int support_node_size   number of nodes support for each gausspoint
 * @param  {[type]} int number_elements     number of elements in the sparse matrix
 */

template <typename T>
bool cpu_assemble_global_matrix(T* global_matrix,
                                int* full_map,
                                int num_cells,
                                int support_node_size,
                                int number_elements)
{
  return true;                  
}
