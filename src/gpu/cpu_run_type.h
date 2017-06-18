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

 #include<../Eigen/Sparse> //eigen
//template <typename T>

using namespace Eigen;

template <typename T>
using SparseMatrix = SparseMatrix<T>;
template <typename T>
using VectorX = Matrix<T, Dynamic, 1>;
/*namespace lmx {
  template <typename T>
using Matrix = Eigen::SparseMatrix<T>;
//using SparseMatrix<float> Matrix<float>;
 //typedef Matrix<data_type, 1, Dynamic> lmx::Vector;//row vector
//template <typename T>
//using Matrix<double, Eigen::Dynamic, 1> Vector<double>;//column vector
//using Matrix<float, Eigen::Dynamic, 1> Vector<float>;//column vector
}*/

#endif //CPU_RUN_TYPE_H
