
#include "LMX/lmx.h"

using namespace lmx;

class SolveTest
{
  public:
    SolveTest()
    {}
    ~SolveTest()
    {}
    void solve( int dim )
    {
      A.resize(0,0);
      A.resize(dim,dim);
      cout << "Solving system of " << dim << " equations... ";
      cout.flush();
      b.resize(dim);

      b.fillIdentity();
      A.fillRandom();

      LinearSystem<double> x(A,b);

      {
        ExactStopwatch sw;
        sw.setQuiet();
        x.solveYourself();
        cout << sw.getTime() << endl;
      }

    }

  private:
    Matrix<double> A;
    Vector<double> b;
};

int main(int argc, char* argv){

  setMatrixType(0);
  setVectorType(0);
  setLinSolverType(0);

  int i;
  SolveTest tester;

  for(i=10; i<=1000; i=10*i){
    tester.solve( i );
  }

  return 1;
} 
