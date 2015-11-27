// #define HAVE_GMM

#include "LMX/lmx.h"
#include "LMX/lmx_nlsolvers.h"

using namespace lmx;

class MyExternalSystem{
public:
  MyExternalSystem(){}

  ~MyExternalSystem(){}

  void residuo(lmx::Vector<double>& res_in, lmx::Vector<double>& q_in){
    // x_i^3 = lhs_i
    if ( lhs.size() == 0 ){
      temp.resize(res_in.size());
      lhs.resize(res_in.size());
      for (int i=0; i<lhs.size(); ++i)
        lhs.writeElement( std::pow(i+1.,6), i); // lhs_i = i^6
    }
    temp.multElements( q_in, q_in ); //(x_i)^2
    res_in.multElements( temp, q_in ); //(x_i)^3
    res_in -= lhs;
  }

  void jacobiano(lmx::Matrix<double>& jac_in, lmx::Vector<double>& q_in){
    for (int i=0; i<q_in.size() ; ++i)
      jac_in.writeElement( 3.*q_in.readElement(i)*q_in.readElement(i), i, i);
  }

  bool convergence( lmx::Vector<double>& res_in ){
    if (res_in.norm1() < 0.001) return 1;
    else return 0;
  }

  private:
    lmx::Vector<double> lhs;
    lmx::Vector<double> temp;
};

int main(int argc, char* argv[])
{
  setMatrixType(0);
  setVectorType(0);
  setLinSolverType(0);

  Vector<double> b(3);

  b.fillIdentity();

  MyExternalSystem theSystem;
  lmx::NLSolver<MyExternalSystem> theSolver;
  theSolver.setInitialConfiguration( b );

  theSolver.setSystem( theSystem );
  theSolver.setResidue( &MyExternalSystem::residuo );
  theSolver.setJacobian( &MyExternalSystem::jacobiano );
  theSolver.setConvergence( &MyExternalSystem::convergence );
  theSolver.solve( 100 );
  cout << "Resultado: " << theSolver.getSolution() << endl;
}


