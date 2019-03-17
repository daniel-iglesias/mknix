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

#ifndef LMXDIFF_PROBLEM_SECOND_H
#define LMXDIFF_PROBLEM_SECOND_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_problem_second.h

      \brief DiffProblemSecond class implementation

      Describes an initial value for a dynamic system with an ODE or DAE 
      description.

      This is the base file of lmx_diff systems' manipulation and solution.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "lmx_diff_problem.h"
#include "lmx_diff_integrator_newmark.h"

namespace lmx {

    /**
    \class DiffProblemSecond 
    \brief Template class DiffProblemSecond.
    Implementation for Fist Order ODE system solvers.

    This class implements methods for defining and solving initial value 
    problems described by a TotalDiff class' derivided object, and initial 
    conditions in the form \f$ q(t_o) = q_o \f$.

    @author Daniel Iglesias.
    */
template <typename Sys, typename T=double> 
class DiffProblemSecond
 : public DiffProblem<Sys, T>{

  public:

    /** Empty constructor. */
    DiffProblemSecond()
     : b_solveInitialEquilibrium(1)
       , b_residueByParts(0)
       , b_jacobianByParts(0)
       , b_alpha(0)
       , b_convergence(0)
       , b_sparse(0)
    {}

    /** Destructor. */
    ~DiffProblemSecond()
    {
      for( unsigned int i=0; i<residueParts.size(); ++i){
        delete residueParts[i];
        residueParts[i]=0;
      }
      for( unsigned int i=0; i<jacobianParts.size(); ++i){
        delete jacobianParts[i];
        jacobianParts[i]=0;
      }
    this->p_delta_q = 0;
    }

    void setResidue
        ( void (Sys::* residue_in)( lmx::Vector<T>& residue,
                                    const lmx::Vector<T>& q,
                                    const lmx::Vector<T>& qdot,
                                    const lmx::Vector<T>& qddot,
                                    double time
                                  )
        );

    void setResidueByParts
        ( void (Sys::* residue_q_qdot)( lmx::Vector<T>& res_q_qdot,
                                        const lmx::Vector<T>& q,
                                        const lmx::Vector<T>& qdot
                                      ),
        void (Sys::* residue_qddot)( lmx::Vector<T>& res_qddot,
                                     const lmx::Vector<T>& qddot
                                   ),
        void (Sys::* residue_time)( lmx::Vector<T>& res_time,
                                    double time
                                  )
        );

    void setJacobian
        ( void (Sys::* jacobian_in)(
                                     lmx::Matrix<T>& jacobian,
                                     const lmx::Vector<T>& q,
                                     const lmx::Vector<T>& qdot,
                                     double partial_qdot,
                                     double partial_qddot,
                                     double time
                                   )
        );

    void setJacobianByParts
     ( void (Sys::* jacobian_q_qdot)( lmx::Matrix<T>& jac_q_qdot,
                                      const lmx::Vector<T>& q,
                                      const lmx::Vector<T>& qdot,
                                      double partial_qdot,
                                      double time
                                    ),
       void (Sys::* jacobian_qddot)( lmx::Matrix<T>& jac_qddot,
                                     double partial_qddot,
                                     double time
                                   )
     );

    void setEvaluation
     ( void (Sys::* eval_in)( const lmx::Vector<T>& q,
                              const lmx::Vector<T>& qdot,
                              lmx::Vector<T>& qddot,
                              double time
                            )
     );

    /**
     * Defines the integrator that will be used for configuration advance & actualization.
     * @param type Key of integrator family to use.
     * @param opt1 Optional value for some integrators (usually specifies the order).
     * @param opt2 Optional value for some integrators.
     */
    void setIntegrator( int type, int opt1=0, int opt2=0 )
    { DiffProblem<Sys, T>::setIntegrator( type, opt1, opt2 ); }

    /**
     * Defines the integrator that will be used for configuration advance & actualization.
     * @param type Key of integrator family to use.
     * @param opt2 Optional value for some integrators.
     */
    void setIntegrator( const char* type, int opt2=0 )
    { DiffProblem<Sys, T>::setIntegrator( type, opt2 ); }

    void setIntegrator( const char* type, double alpha_in );


    void setIntegrator( const char* type,
                        double beta,
                        double gamma,
                        double alpha = 0
                      );

    /**
     * Backward resolution of mother function.
     * @param L2 norm maximum residual.
     */
    void setConvergence( double eps_in )
    { DiffProblem<Sys, T>::setConvergence( eps_in ); }

    void setConvergence
        ( bool (Sys::* convergence)( const lmx::Vector<T>& q,
                                     const lmx::Vector<T>& qdot,
                                     const lmx::Vector<T>& qddot,
                                     double time
                                   )

        );

    // Needs documentation
    void setSparsePatternJacobian( lmx::DenseMatrix<T>& mat_sparse )
    { mat_sparse.writeSparsePattern( v_rows, v_cols );
      b_sparse = 1;
    }

    void iterationResidue( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q );

    void iterationResidueByParts( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q );

    void iterationResidueForAlpha( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q );

    void iterationJacobian( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q );

    void iterationJacobianByParts( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q );

    void iterationJacobianForAlpha( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q );

    bool iterationConvergence( lmx::Vector<T>& delta_q );

    void initialize( );
    
    void solve( );

    void stepSolve( );
    
  private:
    void stepSolveExplicit( );
    void stepSolveImplicit( );

  private:
    bool b_solveInitialEquilibrium; ///< default TRUE.
    bool b_residueByParts; ///< 0 if setResidue is called, 1 if setResidueByParts is called.
    bool b_jacobianByParts; ///< 0 if setJacobian is called, 1 if setJacobianByParts is called.
    bool b_alpha; ///< 1 if HHT-alpha integrator is set.
    bool b_convergence; ///< 1 if external convergence function is set.
    bool b_sparse; ///< 1 if sparse pattern is defined for jacobian matrix. Only for implicit integrators.
    double alpha;
    std::vector< lmx::Vector<T>* > residueParts;
    std::vector< lmx::Matrix<T>* > jacobianParts;
    std::vector<size_type> v_rows, v_cols; // vectors for sparse pattern of jacobian matrix.
    NLSolver< DiffProblemSecond<Sys, T>, T > theNLSolver;
    void (Sys::* res)( lmx::Vector<T>& residue,
                       const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       const lmx::Vector<T>& qddot,
                       double time
                     );
    void (Sys::* res_q_qdot)( lmx::Vector<T>& res_q_qdot,
                              const lmx::Vector<T>& q,
                              const lmx::Vector<T>& qdot
                            );
    void (Sys::* res_qddot)( lmx::Vector<T>& res_qddot,
                             const lmx::Vector<T>& qddot
                           );
    void (Sys::* res_time)( lmx::Vector<T>& res_time,
                            double time
                          );
    void (Sys::* jac)( lmx::Matrix<T>& jacobian,
                       const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       double partial_qdot,
                       double partial_qddot,
					   double time
                     );
   void (Sys::* jac_q_qdot)( lmx::Matrix<T>& jac_q_qdot,
                             const lmx::Vector<T>& q,
                             const lmx::Vector<T>& qdot,
                             double partial_qdot,
							 double time
                           );
   void (Sys::* jac_qddot)( lmx::Matrix<T>& jac_qddot,
                            double partial_qddot,
							double time
                          );
   void (Sys::* eval)( const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       lmx::Vector<T>& qddot,
                       double time
                     );
   bool (Sys::* conv)( const lmx::Vector<T>& q,
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
    void DiffProblemSecond<Sys,T>::setResidue
        ( void (Sys::* residue_in)(
                                    lmx::Vector<T>& residue,
                                    const lmx::Vector<T>& q,
                                    const lmx::Vector<T>& qdot,
                                    const lmx::Vector<T>& qddot,
                                    double time
                                  )
        )
{
  this->b_residueByParts = 0;
  this->res = residue_in;
}

/**
 * Sets the external function for residue evaluation. Must be a Sys member function.
 * @param residue_q_qdot Residue member function depending on \f$ q \f$ and \f$ \frac{dq}{dt} \f$ .
 * @param residue_qddot Residue member function depending on \f$ \frac{d^2q}{dt^2} \f$ .
 * @param residue_time Residue member function depending on time, \f$ t \f$ .
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setResidueByParts
        (
          void (Sys::* residue_q_qdot)( lmx::Vector<T>& res_q_qdot,
                                        const lmx::Vector<T>& q,
                                        const lmx::Vector<T>& qdot
                                      ),
          void (Sys::* residue_qddot)( lmx::Vector<T>& res_qddot,
                                       const lmx::Vector<T>& qddot
                                     ),
          void (Sys::* residue_time)( lmx::Vector<T>& res_time,
                                     double time
                                    )
        )
{
  this->b_residueByParts = 1;
  this->res_q_qdot = residue_q_qdot;
  this->res_qddot = residue_qddot;
  this->res_time = residue_time;
}

/**
 * Sets the external function for tangent to q. Must be a Sys member function.
 * ::: change documentation :::
 * @param jacobian_in Tangent function.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setJacobian
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
  this->b_jacobianByParts = 0;
  this->jac = jacobian_in;
}

/**
 * Sets the external function for residue evaluation. Must be a Sys member function.
 * ::: change documentation :::
 * @param jacobian_q_qdot Tangent member function depending on \f$ q \f$ and \f$ \frac{dq}{dt} \f$ .
 * @param jacobian_qddot Tangent member function depending on \f$ \frac{d^2q}{dt^2} \f$ .
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setJacobianByParts
     ( void (Sys::* jacobian_q_qdot)( lmx::Matrix<T>& jac_q_qdot,
                                      const lmx::Vector<T>& q,
                                      const lmx::Vector<T>& qdot,
                                      double partial_qdot,
									  double time
                                    ),
       void (Sys::* jacobian_qddot)( lmx::Matrix<T>& jac_qddot,
                                     double partial_qddot,
									 double time
                                   )
     )
{
  this->b_jacobianByParts = 1;
  this->jac_q_qdot = jacobian_q_qdot;
  this->jac_qddot  = jacobian_qddot;
}

/**
 * Sets the external (first order configuration) evaluation function. Must be a Sys member function.
 * @param eval_in The acceleration function.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setEvaluation
        ( void (Sys::* eval_in)(
                                 const lmx::Vector<T>& q,
                                 const lmx::Vector<T>& qdot,
                                 lmx::Vector<T>& qddot,
                                 double time
                               )
        )
{
  this->eval = eval_in;
}

/**
 * Defines the integrator that will be used for configuration advance & actualization (Alpha version).
 * @param type Key of integrator family to use.
 * @param alpha_in HHT's alpha.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setIntegrator
        ( const char* type,
          double alpha_in
        )
{
  if (!strcmp(type, "ALPHA")){
    this->b_alpha = 1;
    this->alpha = alpha_in;
    this->theIntegrator = new IntegratorNEWMARK<T>( .25*std::pow(1+alpha,2), .5+alpha );
  }
}

/**
 * Defines the integrator that will be used for configuration advance & actualization (Newmark and Alpha version).
 * @param type Key of integrator family to use.
 * @param beta_in Newmark's beta.
 * @param gamma_in Newmark's gamma.
 * @param alpha_in HHT's alpha.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setIntegrator
        ( const char* type,
          double beta_in,
          double gamma_in,
          double alpha_in
        )
{
  if (!strcmp(type, "ALPHA")){
    this->b_alpha = 1;
    this->alpha = alpha_in;
    this->theIntegrator = new IntegratorNEWMARK<T>( beta_in*std::pow(1+alpha,2), gamma_in+alpha );
  }
  else if (!strcmp(type, "NEWMARK"))
    this->theIntegrator = new IntegratorNEWMARK<T>( beta_in/**std::pow(1+alpha,2)*/, gamma_in/*+alpha*/ );
}

/**
 * Sets the external function that implements a different convergence criteria from those available in LMX.
 * Must be a Sys member function.
 *
 * @param conv_in The convergence evaluation function.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::setConvergence
        ( bool (Sys::* conv_in)( const lmx::Vector<T>& q,
                                 const lmx::Vector<T>& qdot,
                                 const lmx::Vector<T>& qddot,
                                 double time
                               )

        )
{
  this->conv = conv_in;
  b_convergence = 1;
}


/**
 * Function for NLSolver residue computation.
 * @param residue Residue vector.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::
        iterationResidue( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q )
{
  static_cast< IntegratorBaseImplicit<T>* >
      (this->theIntegrator)->integratorUpdate( delta_q );
  (this->theSystem->*res)( residue,
                           this->theConfiguration->getConf(0),
                           this->theConfiguration->getConf(1),
                           this->theConfiguration->getConf(2),
                           this->theConfiguration->getTime( )
                         );

}

/**
 * Function for NLSolver residue computation. Used when ResidueByParts is set.
 * @param residue Residue vector.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::
        iterationResidueByParts( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q )
{
  static_cast< IntegratorBaseImplicit<T>* >
      (this->theIntegrator)->integratorUpdate( delta_q );

  (this->theSystem->*res_q_qdot)( *(residueParts[0]),
                                  this->theConfiguration->getConf(0),
                                  this->theConfiguration->getConf(1)
                                );
  (this->theSystem->*res_qddot)( *(residueParts[1]),
                                 this->theConfiguration->getConf(2)
                               );
  (this->theSystem->*res_time)( *(residueParts[2]),
                                this->theConfiguration->getTime( )
                              );
  residue = *residueParts[0] + *residueParts[1] + *residueParts[2];
}

/**
 * Function for NLSolver residue computation. Used when ALPHA integrator is set.
 * @param residue Residue vector.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::
        iterationResidueForAlpha( lmx::Vector<T>& residue, lmx::Vector<T>& delta_q )
{
  static_cast< IntegratorBaseImplicit<T>* >
      (this->theIntegrator)->integratorUpdate( delta_q );

  (this->theSystem->*res_q_qdot)( *(residueParts[0]),
                                  this->theConfiguration->getConf(0),
                                  this->theConfiguration->getConf(1)
                                );
  (this->theSystem->*res_qddot)( *(residueParts[1]),
                                 this->theConfiguration->getConf(2)
                               );
  (this->theSystem->*res_time)( *(residueParts[2]),
                                this->theConfiguration->getTime( )
                                    - alpha * this->theConfiguration->getLastStepSize()
                              );
  residue = (T)(1.-alpha) * (*residueParts[0])
      +     (T)   (alpha) * (*residueParts[3])
      +                  (*residueParts[1])
      +                   *residueParts[2];
}

/**
 * Function for NLSolver jacobian computation.
 * @param jacobian Tangent matrix.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::iterationJacobian( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q )
{
  (this->theSystem->*jac)( jacobian,
              this->theConfiguration->getConf(0),
              this->theConfiguration->getConf(1),
              static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQdot( ),
              static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQddot( ),
              this->theConfiguration->getTime( )
            );
}

/**
 * Function for NLSolver jacobian computation. Used when ResidueByParts is set.
 * @param jacobian Tangent matrix.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::iterationJacobianByParts( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q )
{
  (this->theSystem->*jac_q_qdot)( *(jacobianParts[0]),
                                  this->theConfiguration->getConf(0),
                                  this->theConfiguration->getConf(1),
                                  static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQdot( ),
                                  this->theConfiguration->getTime( )
                                );
  (this->theSystem->*jac_qddot)( *(jacobianParts[1]),
                                 static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQddot( ),
                                 this->theConfiguration->getTime( )
                               );
  jacobian = *jacobianParts[0] + *jacobianParts[1];
}

/**
 * Function for NLSolver jacobian computation. Used when ALPHA integrator is set.
 * @param jacobian Tangent matrix.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::iterationJacobianForAlpha( lmx::Matrix<T>& jacobian, lmx::Vector<T>& delta_q )
{
  (this->theSystem->*jac_q_qdot)( *(jacobianParts[0]),
                                  this->theConfiguration->getConf(0),
                                  this->theConfiguration->getConf(1),
                                  static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQdot( ),
                                  this->theConfiguration->getTime( )
                                );
  (this->theSystem->*jac_qddot)( *(jacobianParts[1]),
                                 static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQddot( ),
								 this->theConfiguration->getTime( )
                               );
  jacobian = (T)(1-alpha)* (*jacobianParts[0]) + *jacobianParts[1];
}

/**
 * Function for NLSolver convergence evaluation.
 * @param delta_q Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    bool DiffProblemSecond<Sys,T>::iterationConvergence( lmx::Vector<T>& delta_q )
{
  return (this->theSystem->*conv)( this->theConfiguration->getConf(0),
                                   this->theConfiguration->getConf(1),
                                   this->theConfiguration->getConf(2),
                                   this->theConfiguration->getTime( )
                                 );
}


/**
 * Initialize solving function
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::initialize( )
{
  this->theConfiguration->setTime( this->to );
  this->theIntegrator->initialize( this->theConfiguration );
  if ( ! this->theIntegrator->isExplicit() ){
    if (b_solveInitialEquilibrium) // default TRUE
      (this->theSystem->*eval)( this->theConfiguration->getConf(0),
				this->theConfiguration->getConf(1),
				this->theConfiguration->setConf(2),
				this->theConfiguration->getTime( )+this->stepSize
			      );
    if( b_residueByParts ){
      residueParts.push_back( new lmx::Vector<T>(this->theConfiguration->getConf(0).size() ) );
      residueParts.push_back( new lmx::Vector<T>(this->theConfiguration->getConf(0).size() ) );
      residueParts.push_back( new lmx::Vector<T>(this->theConfiguration->getConf(0).size() ) );
      if(b_alpha){
	// Must create the residueParts[3] for storing the residue_q_qdot_{n-1}
	residueParts.push_back( new lmx::Vector<T>(this->theConfiguration->getConf(0).size() ) );
	// First evaluation of the residue in t=to as the first step needs residue_q_qdot_{1-1}
	(this->theSystem->*res_q_qdot)( *(residueParts[3]),
					this->theConfiguration->getConf(0),
					this->theConfiguration->getConf(1)
				      );
      }
    }
    if( b_jacobianByParts ){
      jacobianParts.push_back( new lmx::Matrix<T>
			      (this->theConfiguration->getConf(0).size(),
				this->theConfiguration->getConf(0).size()
			      )
			    );
      jacobianParts.push_back( new lmx::Matrix<T>
			      (this->theConfiguration->getConf(0).size(),
				this->theConfiguration->getConf(0).size()
			      )
			    );
    }
    if (b_sparse) theNLSolver.setSparse( v_rows, v_cols );
    theNLSolver.setInitialConfiguration( this->theConfiguration->getConf(0) );
    theNLSolver.setDeltaInResidue(  );
    theNLSolver.setSystem( *this );
    if( b_convergence ){
      theNLSolver.setConvergence( &DiffProblemSecond<Sys,T>::iterationConvergence );
    }
    if( b_residueByParts ){
      theNLSolver.setResidue( &DiffProblemSecond<Sys,T>::iterationResidueByParts ); // Also advances the integrator
      if( b_alpha )
	theNLSolver.setResidue( &DiffProblemSecond<Sys,T>::iterationResidueForAlpha ); // Also advances the integrator
    }
    else
      theNLSolver.setResidue( &DiffProblemSecond<Sys,T>::iterationResidue ); // Also advances the integrator
    if( b_jacobianByParts ){
      theNLSolver.setJacobian( &DiffProblemSecond<Sys,T>::iterationJacobianByParts );
      if( b_alpha )
	theNLSolver.setJacobian( &DiffProblemSecond<Sys,T>::iterationJacobianForAlpha );
    }
    else
      theNLSolver.setJacobian( &DiffProblemSecond<Sys,T>::iterationJacobian );
    
    this->p_delta_q = &( theNLSolver.getSolution() );
  } 
  this->writeStepFiles();
}


/**
 * Solve main function
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::solve( )
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
    void DiffProblemSecond<Sys,T>::stepSolve( )
{
  if ( this->theIntegrator->isExplicit() )
    this->stepSolveExplicit();
  else this->stepSolveImplicit();
}

/**
 * Explicit time scheme solver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::stepSolveExplicit( )
{
  (this->theSystem->*eval)( this->theConfiguration->getConf(0),
			    this->theConfiguration->getConf(1),
			    this->theConfiguration->setConf(2),
			    this->theConfiguration->getTime( )
		    );
  this->theConfiguration->nextStep( this->stepSize );
  this->theIntegrator->advance( );
  if(this->b_steptriggered) (this->theSystem->*(this->stepTriggered))( );
  this->writeStepFiles();
}

/**
 * Implicit time scheme solver.
 */
template <typename Sys, typename T>
    void DiffProblemSecond<Sys,T>::stepSolveImplicit( )
{
  this->theConfiguration->nextStep( this->stepSize );
  this->theIntegrator->advance( );
  theNLSolver.solve( );
  if(b_alpha) *residueParts[3] = *residueParts[0];
  if(this->b_steptriggered) (this->theSystem->*(this->stepTriggered))( );
  this->writeStepFiles();
}

}; // namespace lmx


#endif
