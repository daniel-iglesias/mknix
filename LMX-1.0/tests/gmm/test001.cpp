// Tests harwellBoeing save and load

#include "LMX/lmx.h"

using namespace lmx;

int main(int argc, char* argv[]){

  setMatrixType(0);
  setVectorType(0);

  Matrix<double> aMatrix(4,4);
  
  aMatrix(0,0) = 10.;
  aMatrix(0,1) = 1.;
  aMatrix(1,0) = 1.;
  aMatrix(1,1) = 20.;
  aMatrix(2,1) = 2.;
  aMatrix(1,2) = 2.;
  aMatrix(2,2) = 30.;
  aMatrix(2,3) = 3.;
  aMatrix(3,2) = 3.;
  aMatrix(3,3) = 40.;
  
  cout << aMatrix;
  
  aMatrix.harwellBoeingSave("test001.mat");
  
  Matrix<double> loadedMatrix;
  loadedMatrix.harwellBoeingLoad("test001.mat");
  
  cout << loadedMatrix;
  
  return 1;
} 
