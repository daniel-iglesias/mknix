// #define HAVE_GMM

#include "LMX/lmx.h"

using namespace lmx;

int main(int argc, char* argv[])
{
  setMatrixType(0);
  setVectorType(0);
  setLinSolverType(0);

  Vector<double> a(3);
  Vector<double> b(3);
  Vector<double> c(3);
  double norm;
  
  a.writeElement(1.,0);
  a.writeElement(2.,0);
  b.writeElement(3.5,1);
  
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;

  c.add(a,b);
  std::cout << "a+b = " << c << std::endl;
  std::cout << "a+b = " << a+b << std::endl;

  c.subs(a,b);
  std::cout << "a-b = " << c << std::endl;
  std::cout << "a-b = " << a-b << std::endl;

  c.multElements(a,b);
  std::cout << "a(i)*b(i) = " << c << std::endl;
  std::cout << "a*b = " << a*b << std::endl;


  norm = a.norm1();
  cout << a << endl << "norm(a) = " << norm << endl;
  
} 
