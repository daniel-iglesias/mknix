#ifndef CUDA_HELPER_H
#define CUDA_HELPER_H
#include <iostream>
#include <iomanip>
#include <cuda.h>
//#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

namespace CudaHelper{

    template <typename T>
    bool divide_array_by_value(T*  array,
                               T value,
                               int size,
                               int threads_per_block,
                               cudaStream_t stream = 0);

    bool exclusive_scan(int *input_first,
                        int *input_last,
                        int *output,
                        int init);

    size_t free_memory_gpu();

    template <typename T>
    bool init_array_to_value(T *array,
                             T value,
                             int size,
                             int threads_per_block,
                             cudaStream_t stream = 0);

    bool init_index_array(int *index_array,
                          int size,
                          int threads_per_block,
                          cudaStream_t stream = 0);

    template <typename T, typename T2>
    bool multiply_arrays_by_element(T *array,
                                    T2 *array2,
                                    int size,
                                    int threads_per_block,
                                    cudaStream_t stream = 0);

    void print_gpu_memory();

}


#endif //CUDA_HELPER_H
