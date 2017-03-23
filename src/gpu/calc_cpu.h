#ifndef CALC_CPU_H
#define CALC_CPU_H

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "LMX/lmx.h"
#include "gmm/gmm_matrix.h"

#include <atomic>
#include <pthread.h>
//struct for multithreaded launch
struct p_struct_c{
  /*std::atomic<double>* globalMatrix;
  std::vector<int> *fullMap;
  double *local_matrices_array;
  int numCells;
  int supportNodeSize;
  int thread_id;
  int max_threads;*/
};
/*
void atomicAssembleGlobalMatrix(std::atomic<double>* globalMatrix,
                                std::vector<int> &fullMap,
                                double *local_matrices_array,
                                int numPoints,
                                int supportNodeSize,
                                int tid,
                                int max_threads);*/

void* threadWrapper(void* ptr);

calculateCapacityMatrix

calculateConductivityMatrix

//templating part


#endif //CALC_CPU_H
