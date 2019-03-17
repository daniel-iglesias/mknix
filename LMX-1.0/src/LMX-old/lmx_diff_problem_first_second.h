/***************************************************************************
 *   Copyright (C) 2005 by Daniel Iglesias                                 *
 *   http://code.google.com/p/lmx                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef LMXDIFF_PROBLEM_FIRST_SECOND_H
#define LMXDIFF_PROBLEM_FIRTS_SECOND_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_problem_first_second.h

      \brief DiffProblemFirstSecond class implementation

      Describes an initial value for a partitioned dynamic system with an ODE or DAE 
      description. The system is composed by two subsystem of orders 1 and 2, respectively.

      This is the base file of lmx_diff systems' manipulation and solution.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "lmx_diff_problem_double.h"
#include "lmx_diff_integrator_newmark.h"

namespace lmx {

    /**
    \class DiffProblemFirstSecond 
    \brief Template class DiffProblemFirstSecond.
    Implementation for weakly coupled First-Second Order ODE system solvers.

    This class implements methods for defining and solving initial value 
    problems described by a TotalDiff class' derivided object, and initial 
    conditions in the form \f$ q(t_o) = q_o \f$.

    @author Daniel Iglesias.
    */
template <typename Sys, typename T=double> 
class DiffProblemFirstSecond
 : public DiffProblemDouble<Sys, T>{

  public:

    /** Empty constructor. */
    DiffProblemFirstSecond()
     : b_solveInitialEquilibrium(1)
//        , b_alpha(0)
       , b_convergence1(0)
       , b_convergence2(0)
       , b_sparse1(0)
       , b_sparse2(0)
    {}

    /** Destructor. */
    ~DiffProblemFirstSecond()
    {
      this->p_delta_q1 = 0;
      this->p_delta_q2 = 0;
    }

    void setResidue1( void (Sys::* residue_in)( lmx::Vector<T>& residue1,
                                               const lmx::Vector<T>& q1,
                                               const lmx::Vector<T>& qdot1,
                                               double time
                                             )
                    );
    void setResidue2
        ( void (Sys::* residue_in)( lmx::Vector<T>& residue2,
                                    const lmx::Vector<T>& q22,
                                    const lmx::Vector<T>& qdot2,
                                    const lmx::Vector<T>& qddot2,
                                    double time
                                  )
        );

    void setJacobian1( void (Sys::* jacobian_in)( lmx::Matrix<T>& tangent,
                                                 const lmx::Vector<T>& q,
                                                 double partial_qdot,
                                                 double time
                                               )
                     );

    void setJacobian2
        ( void (Sys::* jacobian_in)(
                                     lmx::Matrix<T>& jacobian,
                                     const lmx::Vector<T>& q2,
                                     const lmx::Vector<T>& qdot2,
                                     double partial_qdot2,
                                     double partial_qddot2,
                                     double time
                                   )
        );

    void setEvaluation1( void (Sys::* eval_in)( const lmx::Vector<T>& q1,
                                               lmx::Vector<T>& qdot1,
                                               double time
                                             )
                       );
    void setEvaluation2
     ( void (Sys::* eval_in)( const lmx::Vector<T>& q2,
                              const lmx::Vector<T>& qdot2,
                              lmx::Vector<T>& qddot2,
                              double time
                            )
     );

    /**
     * Defines the integrator that will be used for second order system configuration advance & updating.
     * @param type Key of integrator family to use.
     * @param opt1 Optional value for some integrators (usually specifies the order).
     * @param opt2 Optional value for some integrators.
     */
    void setIntegrator2( int type, int opt1=0, int opt2=0 )
    { DiffProblemDouble<Sys, T>::setIntegrator2( type, opt1, opt2 ); }

    /**
     * Defines the integrator that will be used for configuration advance & updating.
     * @param type Key of integrator family to use.
     * @param opt2 Optional value for some integrators.
     */
    void setIntegrator2( const char* type, int opt2=0 )
    { DiffProblemDouble<Sys, T>::setIntegrator2( type, opt2 ); }

//     void setIntegrator( const char* type, double alpha_in );


    void setIntegrator2( const char* type,
                        double beta,
                        double gamma,
                        double alpha = 0
                      );

    /**
     * Backward resolution of mother function.
     * @param L2 norm maximum residual.
     */
    void setConvergence( double eps_in )
    { DiffProblemDouble<Sys, T>::setConvergence( eps_in ); }

    void setConvergence1
        ( bool (Sys::* convergence1)( const lmx::Vector<T>& q,
                                     const lmx::Vector<T>& qdot,
                                     double time
                                   )

        );
    void setConvergence2
        ( bool (Sys::* convergence2)( const lmx::Vector<T>& q,
                                     const lmx::Vector<T>& qdot,
                                     const lmx::Vector<T>& qddot,
                                     double time
                                   )

        );

    // Needs documentation
    void setSparsePatternJacobian1( lmx::DenseMatrix<T>& mat_sparse1 )
    { mat_sparse1.writeSparsePattern( v_rows1, v_cols1 );
      b_sparse1 = 1;
    }
    void setSparsePatternJacobian2( lmx::DenseMatrix<T>& mat_sparse2 )
    { mat_sparse2.writeSparsePattern( v_rows2, v_cols2 );
      b_sparse2 = 1;
    }

    void iterationResidue1( lmx::Vector<T>& residue, lmx::Vector<T>& q_current );

    void iterationJacobian1( lmx::Matrix<T>& jacobian, lmx::Vector<T>& q_current );

    bool iterationConvergence1( lmx::Vector<T>& q_current );

    void iterationResidue2( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q );

    void iterationJacobian2( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q );

    bool iterationConvergence2( lmx::Vector<T>& delta_q );

    void initialize( );
    
    void solve( );

    void stepSolve( );
    
  private:
    void stepSolveExplicit( );
    void stepSolveImplicit( );

  private:
    bool b_solveInitialEquilibrium; ///< default TRUE.
//     bool b_alpha; ///< 1 if HHT-alpha integrator is set.
    bool b_convergence1; ///< 1 if external convergence function is set for 1st order system.
    bool b_convergence2; ///< 1 if external convergence function is set for 2nd order system.
    bool b_sparse1; ///< 1 if sparse pattern is defined for for 1st order jacobian matrix. Only for implicit integrators.
    bool b_sparse2; ///< 1 if sparse pattern is defined for for 2nd order jacobian matrix. Only for implicit integrators.
//     double alpha;
    std::vector<size_type> v_rows1, v_cols1; // vectors for sparse pattern of 1st order jacobian matrix.
    std::vector<size_type> v_rows2, v_cols2; // vectors for sparse pattern of 2nd order jacobian matrix.
    NLSolverDouble< DiffProblemFirstSecond<Sys, T>, T > theNLSolver;
    void (Sys::* res1)( lmx::Vector<T>& residue,
                       const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       double time
                      );
    void (Sys::* res2)( lmx::Vector<T>& residue,
                       const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       const lmx::Vector<T>& qddot,
                       double time
                      );
    void (Sys::* jac1)( lmx::Matrix<T>& jacobian,
                       const lmx::Vector<T>& q,
                       double partial_qdot,
		       double time
                      );
    void (Sys::* jac2)( lmx::Matrix<T>& jacobian,
                       const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       double partial_qdot,
                       double partial_qddot,
		       double time
                      );
    void (Sys::* eval1)( const lmx::Vector<T>& q,
                         lmx::Vector<T>& qdot,
                         double time
                        );
    void (Sys::* eval2)( const lmx::Vector<T>& q,
                         const lmx::Vector<T>& qdot,
                         lmx::Vector<T>& qddot,
                         double time
                       );
    bool (Sys::* conv1)( const lmx::Vector<T>& q,
                        const lmx::Vector<T>& qdot,
                        double time
                       );
    bool (Sys::* conv2)( const lmx::Vector<T>& q,
                        const lmx::Vector<T>& qdot,
                        const lmx::Vector<T>& qddot,
                        double time
                       );

};


/////////////////////////////// Implementation of the methods defined previously


/**
 * Sets the external function for residue evaluation. Must be a Sys member function.
 * @param residue_in Residue function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::
        setResidue1( void (Sys::* residue_in)( lmx::Vector<T>& residue,
                                             const lmx::Vector<T>& q,
                                             const lmx::Vector<T>& qdot,
                                             double time
                                           )
                  )
{
  this->res1 = residue_in;
}

/**
 * Sets the external function for residue evaluation. Must be a Sys member function.
 * @param residue_in Residue function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::setResidue2
        ( void (Sys::* residue_in)(
                                    lmx::Vector<T>& residue,
                                    const lmx::Vector<T>& q,
                                    const lmx::Vector<T>& qdot,
                                    const lmx::Vector<T>& qddot,
                                    double time
                                  )
        )
{
  this->res2 = residue_in;
}

/**
 * Sets the external function for tangent to q. Must be a Sys member function.
 * :::change documentation:::
 * @param jacobian_in Tangent function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::
        setJacobian1( void (Sys::* jacobian_in)(  lmx::Matrix<T>& tangent,
                                                 const lmx::Vector<T>& q,
                                                 double partial_qdot,
                                                 double time
                                              )
                   )
{
  this->jac1 = jacobian_in;
}


/**
 * Sets the external function for tangent to q. Must be a Sys member function.
 * ::: change documentation :::
 * @param jacobian_in Tangent function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::setJacobian2
        ( void (Sys::* jacobian_in)(
                                     lmx::Matrix<T>& jacobian,
                                     const lmx::Vector<T>& q,
                                     const lmx::Vector<T>& qdot,
                                     double partial_qdot,
                                     double partial_qddot,
                                     double time
                                   )
        )
{
  this->jac2 = jacobian_in;
}

/**
 * Sets the external (first order configuration) evaluation function. Must be a Sys member function.
 * @param eval_in The acceleration function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::
        setEvaluation1( void (Sys::* eval_in)( const lmx::Vector<T>& q,
                                              lmx::Vector<T>& qdot,
                                              double time
                                            )
                     )
{
  this->eval1 = eval_in;
}


/**
 * Sets the external (second order configuration) evaluation function. Must be a Sys member function.
 * @param eval_in The acceleration function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::setEvaluation2
        ( void (Sys::* eval_in)(
                                 const lmx::Vector<T>& q,
                                 const lmx::Vector<T>& qdot,
                                 lmx::Vector<T>& qddot,
                                 double time
                               )
        )
{
  this->eval2 = eval_in;
}

// /**
//  * Defines the integrator that will be used for configuration advance & updating (Alpha version).
//  * @param type Key of integrator family to use.
//  * @param alpha_in HHT's alpha.
//  */
// template <typename Sys, typename T>
//     void DiffProblemFirstSecond<Sys,T>::setIntegrator
//         ( const char* type,
//           double alpha_in
//         )
// {
//   if (!strcmp(type, "ALPHA")){
//     this->b_alpha = 1;
//     this->alpha = alpha_in;
//     this->theIntegrator = new IntegratorNEWMARK<T>( .25*std::pow(1+alpha,2), .5+alpha );
//   }
// }

/**
 * Defines the integrator that will be used for configuration advance & updating (Newmark and Alpha version).
 * @param type Key of integrator family to use.
 * @param beta_in Newmark's beta.
 * @param gamma_in Newmark's gamma.
 * @param alpha_in HHT's alpha.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::setIntegrator2
        ( const char* type,
          double beta_in,
          double gamma_in,
          double alpha_in
        )
{
//   if (!strcmp(type, "ALPHA")){
//     this->b_alpha = 1;
//     this->alpha = alpha_in;
//     this->theIntegrator = new IntegratorNEWMARK<T>( beta_in*std::pow(1+alpha,2), gamma_in+alpha );
//   }
//   else if (!strcmp(type, "NEWMARK"))
    this->theIntegrator2 = new IntegratorNEWMARK<T>( beta_in, gamma_in );
}

/**
 * Sets the external function that implements a different convergence criteria from those available in LMX.
 * Must be a Sys member function.
 * 
 * @param conv_in The convergence evaluation function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::setConvergence1
        ( bool (Sys::* conv_in)( const lmx::Vector<T>& q,
                                 const lmx::Vector<T>& qdot,
                                 double time
                               )

        )
{
  this->conv1 = conv_in;
  b_convergence1 = 1;
}

/**
 * Sets the external function that implements a different convergence criteria from those available in LMX.
 * Must be a Sys member function.
 *
 * @param conv_in The convergence evaluation function.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::setConvergence2
        ( bool (Sys::* conv_in)( const lmx::Vector<T>& q,
                                 const lmx::Vector<T>& qdot,
                                 const lmx::Vector<T>& qddot,
                                 double time
                               )

        )
{
  this->conv2 = conv_in;
  b_convergence2 = 1;
}


/**
 * Function for NLSolver residue computation of 1st order system.
 * @param residue Residue vector.
 * @param q_current Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::iterationResidue1( lmx::Vector<T>& residue, lmx::Vector<T>& q_current )
{
  static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator1)->integratorUpdate( q_current );

  (this->theSystem->*res1)( residue,
              this->theConfiguration1->getConf(0),
              this->theConfiguration1->getConf(1),
              this->theConfiguration1->getTime( )
            );

}

/**
 * Function for NLSolver residue computation of 2nd order system.
 * @param residue Residue vector.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::
        iterationResidue2( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q )
{
  static_cast< IntegratorBaseImplicit<T>* >
      (this->theIntegrator2)->integratorUpdate( delta_q );
  (this->theSystem->*res2)( residue,
                           this->theConfiguration2->getConf(0),
                           this->theConfiguration2->getConf(1),
                           this->theConfiguration2->getConf(2),
                           this->theConfiguration2->getTime( )
                         );
}


/**
 * Function for NLSolver jacobian computation for 1st order system.
 * @param jacobian Tangent matrix.
 * @param q_current Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::iterationJacobian1( lmx::Matrix<T>& jacobian, lmx::Vector<T>& q_current )
{
  (this->theSystem->*jac1)( jacobian,
              this->theConfiguration1->getConf(0),
              static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator1)->getPartialQdot( ),
              this->theConfiguration1->getTime( )
            );
}

/**
 * Function for NLSolver jacobian computation for 2nd order system.
 * @param jacobian Tangent matrix.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::iterationJacobian2( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q )
{
  (this->theSystem->*jac2)( jacobian,
              this->theConfiguration2->getConf(0),
              this->theConfiguration2->getConf(1),
              static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator2)->getPartialQdot( ),
              static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator2)->getPartialQddot( ),
              this->theConfiguration2->getTime( )
            );
}


/**
 * Function for NLSolver convergence evaluation of 1st order system.
 * @param q_current Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    bool DiffProblemFirstSecond<Sys,T>::iterationConvergence1( lmx::Vector<T>& q_current )
{
  return (this->theSystem->*conv1)( this->theConfiguration1->getConf(0),
                                    this->theConfiguration1->getConf(1),
                                    this->theConfiguration1->getTime( )
                                  );
}

/**
 * Function for NLSolver convergence evaluation of 2nd order system.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    bool DiffProblemFirstSecond<Sys,T>::iterationConvergence2( lmx::Vector<T>& delta_q )
{
  return (this->theSystem->*conv2)( this->theConfiguration2->getConf(0),
                                    this->theConfiguration2->getConf(1),
                                    this->theConfiguration2->getConf(2),
                                    this->theConfiguration2->getTime( )
                                  );
}


/**
 * Initialize solving function
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::initialize( )
{
  this->theConfiguration1->setTime( this->to );
  this->theConfiguration2->setTime( this->to );
  this->theIntegrator1->initialize( this->theConfiguration1 );
  this->theIntegrator2->initialize( this->theConfiguration2 );
  std::stringstream message;
  message << "ERROR : Implicit-Explicit integrators defined in DiffProblemFirstSecond." << endl;
  if ( ! this->theIntegrator1->isExplicit() ){
    if ( ! this->theIntegrator2->isExplicit() ){
      if (b_solveInitialEquilibrium){ // default TRUE
	(this->theSystem->*eval1)( this->theConfiguration1->getConf(0),
				  this->theConfiguration1->setConf(1),
				  this->theConfiguration1->getTime( )
				);
	(this->theSystem->*eval2)( this->theConfiguration2->getConf(0),
				  this->theConfiguration2->getConf(1),
				  this->theConfiguration2->setConf(2),
				  this->theConfiguration2->getTime( )
				);
      }
      if (b_sparse1) theNLSolver.setSparse1( v_rows1, v_cols1 );
      if (b_sparse2) theNLSolver.setSparse2( v_rows2, v_cols2 );
      theNLSolver.setInitialConfiguration1( this->theConfiguration1->getConf(0) );
      theNLSolver.setInitialConfiguration2( this->theConfiguration2->getConf(0) );
      theNLSolver.setDeltaInResidues(  );
      theNLSolver.setSystem( *this );
      if( b_convergence1 )
	theNLSolver.setConvergence1( &DiffProblemFirstSecond<Sys,T>::iterationConvergence1 );
      if( b_convergence2 )
	theNLSolver.setConvergence2( &DiffProblemFirstSecond<Sys,T>::iterationConvergence2 );
      
      theNLSolver.setResidue1( &DiffProblemFirstSecond<Sys,T>::iterationResidue1 ); // Also advances the integrator
      theNLSolver.setResidue2( &DiffProblemFirstSecond<Sys,T>::iterationResidue2 ); // Also advances the integrator
      theNLSolver.setJacobian1( &DiffProblemFirstSecond<Sys,T>::iterationJacobian1 );
      theNLSolver.setJacobian2( &DiffProblemFirstSecond<Sys,T>::iterationJacobian2 );
      this->p_delta_q1 = &( theNLSolver.getSolution1() );
      this->p_delta_q2 = &( theNLSolver.getSolution2() );
    }
    else LMX_THROW(lmx::failure_error, message.str() );
  }
  else if ( ! this->theIntegrator2->isExplicit() ){
    LMX_THROW(lmx::failure_error, message.str() );
  }
  this->writeStepFiles();
}


/**
 * Solve main function
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::solve( )
{
  this->initialize();
  int max = (int)( (this->tf - this->to) / this->stepSize );
  for ( int i=0; i<max; ++i){
    this->stepSolve();
  }
}

/**
 * Solve only one step 
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::stepSolve( )
{
  std::stringstream message;
  message << "ERROR : Implicit-Explicit integrators defined in DiffProblemFirstSecond." << endl;
  if ( this->theIntegrator1->isExplicit() ){
    if ( this->theIntegrator2->isExplicit() ){
      this->stepSolveExplicit();
    }
    else 
      LMX_THROW(lmx::failure_error, message.str() );
  }
  else if( this->theIntegrator2->isExplicit() ){
    LMX_THROW(lmx::failure_error, message.str() );
  }
  else this->stepSolveImplicit();
}

/**
 * Explicit time scheme solver.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::stepSolveExplicit( )
{
  (this->theSystem->*eval1)( this->theConfiguration1->getConf(0),
			      this->theConfiguration1->setConf(1),
			      this->theConfiguration1->getTime( )
		      );
  (this->theSystem->*eval2)( this->theConfiguration2->getConf(0),
			      this->theConfiguration2->getConf(1),
			      this->theConfiguration2->setConf(2),
			      this->theConfiguration2->getTime( )+this->stepSize
		      );
  this->writeStepFiles();
  this->theConfiguration1->nextStep( this->stepSize );
  this->theConfiguration2->nextStep( this->stepSize );
  this->theIntegrator1->advance( );
  this->theIntegrator2->advance( );
  if(this->b_steptriggered) (this->theSystem->*(this->stepTriggered))( );
  this->writeStepFiles();
}

/**
 * Implicit time scheme solver.
 */
template <typename Sys, typename T>
    void DiffProblemFirstSecond<Sys,T>::stepSolveImplicit( )
{
  this->writeStepFiles();
  this->theConfiguration1->nextStep( this->stepSize );
  this->theConfiguration2->nextStep( this->stepSize );
  this->theIntegrator1->advance( );
  this->theIntegrator2->advance( );
  theNLSolver.solve( );
  if(this->b_steptriggered) (this->theSystem->*(this->stepTriggered))( );
  this->writeStepFiles();
}

}; // namespace lmx


#endif
