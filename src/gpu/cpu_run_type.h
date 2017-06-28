#ifndef CPU_RUN_TYPE_H
#define CPU_RUN_TYPE_H
//temporary switches
#define OLD_CODE false
#define NEWCPU true
#define MULTICPU false
#define GPU false
//
#define MAX_THREADS 4
//
#define USECSC false

#include <common.h>

#include <../Eigen/Sparse>
using namespace Eigen;

/*template <typename T>
using SparseMatrix = SparseMatrix<T>;*/
template <typename T>
using VectorX = Matrix<T, Dynamic, 1>;

#endif //CPU_RUN_TYPE_H
