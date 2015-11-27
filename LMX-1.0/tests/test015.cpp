// #define HAVE_GMM

#include "LMX/lmx.h"

using namespace lmx;

int main(int argc, char* argv[])
{
  setMatrixType(1);
  setVectorType(0);

  { //Example of "sparsePattern" in CSC matrix:
    Matrix<double> A(3,3);
    Vector<size_type> ja(4);
    Vector<size_type> ia(2);

    ja(0) = 1;
    ja(1) = 1; //no elements in first col
    ja(2) = 2; //one element in second col
    ja(3) = 3; //one element in third col

    ia(0) = 1;//first element in first row
    ia(1) = 2;//second element in second row

    // shapping...
    A.sparsePattern(ia, ja);

    // filling...
    A(0,1) = 3.0;
    A(1,2) = 4.0;
    
    cout << A << endl /*<< b << endl << c*/ << endl;
  }
} 
