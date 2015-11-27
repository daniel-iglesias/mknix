/***************************************************************************
 *   Copyright (C) 2014 by Daniel Iglesias                                 *
 *   daniel@extremo                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// #define HAVE_GMM
// #define HAVE_SUPERLU

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include"LMX/lmx.h"
#include"LMX/lmx_diff_problem_first_second.h"

using namespace std;
//// Weakly coupled systems of different order:
// System1:            q1dot +  H*q1 = g(q2,t)
// System2: q2ddot +(0*q2dot)+ K*q2 = f(q1,t)
class MyDiffSystem{
  public:
    MyDiffSystem(lmx::Vector<double> q1_in, lmx::Vector<double> q2_in)
    {
      localQ1.resize(q1_in.size());
      localQ1 = q1_in;
      localQ2.resize(q2_in.size());
      localQ2 = q2_in;
      
      H.resize(3,3);
      H(0,0) = 2.;
      H(1,1) = 3.;
      H(2,2) = 4.;

      K.resize(2,2);
      K(0,0) = 2.;
      K(1,0) = 1.;
      K(0,1) = 3.;
      K(1,1) = 1.;
    }

    ~MyDiffSystem(){}

    void myEvaluation1( const lmx::Vector<double>& q1,
                        lmx::Vector<double>& q1dot,
                        double time
                      )
    {
      q1dot -= H*q1;
    }

    void myEvaluation2( const lmx::Vector<double>& q2,
                        const lmx::Vector<double>& q2dot,
                        lmx::Vector<double>& q2ddot,
                        double time
                      )
    {
      q2ddot -= K*q2;
    }

    void myResidue1( lmx::Vector<double>& residue,
                    const lmx::Vector<double>& q1,
                    const lmx::Vector<double>& q1dot,
                    double time
                  )
    {
      localQ1 = q1;
      residue = q1dot + H*q1;
      // g(q2, t) = { q2[1]-q2[0]+t, q2[0], q2[1] }
      residue(0) -= ( localQ2.readElement(1)
		     -localQ2.readElement(0)
		     +time);
      residue(1) -= localQ2(0); // same as .readElement(0) but occasionaly less efficient
      residue(1) -= localQ2(1);
    }

    void myResidue2( lmx::Vector<double>& residue,
                    const lmx::Vector<double>& q2,
                    const lmx::Vector<double>& q2dot,
                    const lmx::Vector<double>& q2ddot,
                    double time
                  )
    {
      localQ2 = q2;
      residue = q2ddot + K*q2;
      // f(q1, t) = { q1[0]+t, q1[1]+2*t }
      residue(0) -= localQ1(0)+time;
      residue(1) -= localQ1(1)+2*time;
    }

    void myTangent1( lmx::Matrix<double>& tangent,
                    const lmx::Vector<double>& q1,
                    double partial_q1dot,
		    double time
                  )
    {
      tangent.fillIdentity( partial_q1dot );
      tangent += H;
    }
    void myTangent2( lmx::Matrix<double>& tangent,
                    const lmx::Vector<double>& q2,
                    const lmx::Vector<double>& q2dot,
                    double partial_q2dot,
                    double partial_q2ddot,
		    double time
                  )
    {
      tangent.fillIdentity( partial_q2ddot );
      tangent += K;
    }

  private:
    lmx::Vector<double> localQ1;
    lmx::Vector<double> localQ2;
    lmx::Matrix<double> H;
    lmx::Matrix<double> K;    
};

int main(int argc, char* argv[])
{

  lmx::setMatrixType( 0 );
  lmx::setVectorType( 0 );
  lmx::setLinSolverType( 0 );

  lmx::Vector<double> q1_0(3);
  q1_0(0)=0.2;
  q1_0(1)=0.0;
  q1_0(2)=0.1;
  lmx::Vector<double> q2_0(2);
  q2_0(0)=0.1;
  q2_0(1)=0.0;
  lmx::Vector<double> q2_dot0(2);
  q2_0(0)=0.0;
  q2_0(1)=0.0;
  lmx::DiffProblemFirstSecond< MyDiffSystem > theProblem;
  MyDiffSystem theSystem(q1_0, q2_0);
  
  theProblem.setDiffSystem( theSystem );
  theProblem.setIntegrator1( "BDF-2" );
  theProblem.setIntegrator2( "BDF-1" );
  theProblem.setInitialConfiguration1( q1_0 );
  theProblem.setInitialConfiguration2( q2_0, q2_dot0 );
  theProblem.setTimeParameters( 0, 0.2, 0.01 );
//  theProblem.setOutputFile("res.dat", 0);
  theProblem.setEvaluation1( &MyDiffSystem::myEvaluation1 );
  theProblem.setEvaluation2( &MyDiffSystem::myEvaluation2 );
  theProblem.setResidue1( &MyDiffSystem::myResidue1 );
  theProblem.setResidue2( &MyDiffSystem::myResidue2 );
  theProblem.setJacobian1( &MyDiffSystem::myTangent1 );
  theProblem.setJacobian2( &MyDiffSystem::myTangent2 );
  theProblem.setConvergence( 1E-5 );
  theProblem.solve();
  
  return EXIT_SUCCESS;
}
