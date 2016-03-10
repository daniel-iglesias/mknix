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

#ifndef LMXDIFF_PROBLEM_FIRST_H
#define LMXDIFF_PROBLEM_FIRST_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_problem_first.h

      \brief DiffProblem class implementation

      Describes an initial value for a dynamic system with an ODE or DAE description.

      This is the base file of lmx_diff systems' manipulation and solution.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "lmx_diff_problem.h"

namespace lmx {

    /**
    \class DiffProblemFirst 
    \brief Template class DiffProblemFirst.
    Implementation for First Order ODE system solvers.

    This class implements methods for defining and solving initial value
    problems described by a TotalDiff class' derivided object, and initial
    conditions in the form \f$ q(t_o) = q_o \f$.

    @author Daniel Iglesias .
    */
template <typename Sys, typename T=double> 
class DiffProblemFirst
 : public DiffProblem<Sys, T>{

  public:

    /** Empty constructor. */
    DiffProblemFirst()
     : solveInitialEquilibrium(1)
       , b_convergence(0)
    {}

    /** Destructor. */
    ~DiffProblemFirst()
    {}

    void setResidue( void (Sys::* residue_in)( lmx::Vector<T>& residue,
                                               const lmx::Vector<T>& q,
                                               const lmx::Vector<T>& qdot,
                                               double time
                                             )
                   );
    void setJacobian( void (Sys::* jacobian_in)( lmx::Matrix<T>& tangent,
                                                 const lmx::Vector<T>& q,
                                                 double partial_qdot,
                                                 double time
                                               )
                    );
    void setEvaluation( void (Sys::* eval_in)( const lmx::Vector<T>& q,
                                               lmx::Vector<T>& qdot,
                                               double time
                                             )
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
                                     double time
                                   )

        );

    void setStepTriggered( void (Sys::* stepTriggered_in)() );

    void iterationResidue( lmx::Vector<T>& residue, lmx::Vector<T>& q_current );

    void iterationJacobian( lmx::Matrix<T>& jacobian, lmx::Vector<T>& q_current );

    bool iterationConvergence( lmx::Vector<T>& q_current );

    void initialize( );

    void solve( );

    void stepSolve( );
    
  private:
    void stepSolveExplicit( );
    void stepSolveImplicit( );

  private:
    bool solveInitialEquilibrium; ///<default TRUE
    bool b_convergence; ///< 1 if external convergence function is set.
    NLSolver< DiffProblemFirst<Sys, T> > theNLSolver;
    void (Sys::* stepTriggered)(); ///< function called at the end of each time step
    void (Sys::* res)( lmx::Vector<T>& residue,
                       const lmx::Vector<T>& q,
                       const lmx::Vector<T>& qdot,
                       double time
                     );
    void (Sys::* jac)( lmx::Matrix<T>& tangent,
                       const lmx::Vector<T>& q,
                       double partial_qdot,
                       double time
                     );
    void (Sys::* eval)( const lmx::Vector<T>& q,
                        lmx::Vector<T>& qdot,
                        double time
                      );
    bool (Sys::* conv)( const lmx::Vector<T>& q,
                        const lmx::Vector<T>& qdot,
                        double time
                      );

};


/////////////////////////////// Implementation of the methods defined previously


/**
 * Sets the external function for residue evaluation. Must be a Sys member function.
 * @param residue_in Residue function.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::
        setResidue( void (Sys::* residue_in)( lmx::Vector<T>& residue,
                                             const lmx::Vector<T>& q,
                                             const lmx::Vector<T>& qdot,
                                             double time
                                           )
                  )
{
  this->res = residue_in;
}

/**
 * Sets the external function for tangent to q. Must be a Sys member function.
 * :::change documentation:::
 * @param jacobian_in Tangent function.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::
        setJacobian( void (Sys::* jacobian_in)(  lmx::Matrix<T>& tangent,
                                                 const lmx::Vector<T>& q,
                                                 double partial_qdot,
                                                 double time
                                              )
                   )
{
  this->jac = jacobian_in;
}

/**
 * Sets the external (first order configuration) evaluation function. Must be a Sys member function.
 * @param eval_in The acceleration function.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::
        setEvaluation( void (Sys::* eval_in)( const lmx::Vector<T>& q,
                                              lmx::Vector<T>& qdot,
                                              double time
                                            )
                     )
{
  this->eval = eval_in;
}

/**
 * Sets the external function that implements a different convergence criteria from those available in LMX.
 * Must be a Sys member function.
 * 
 * @param conv_in The convergence evaluation function.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::setConvergence
        ( bool (Sys::* conv_in)( const lmx::Vector<T>& q,
                                 const lmx::Vector<T>& qdot,
                                 double time
                               )

        )
{
  this->conv = conv_in;
  b_convergence = 1;
}


  /**
   * Defines a function call between time steps.
   *
   */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::setStepTriggered( void (Sys::* stepTriggered_in)() )
{
  this->stepTriggered = stepTriggered_in;
  this->b_steptriggered = 1;
}


/**
 * Function for NLSolver residue computation.
 * @param residue Residue vector.
 * @param q_current Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::iterationResidue( lmx::Vector<T>& residue, lmx::Vector<T>& q_current )
{
  static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->integratorUpdate( q_current );

  (this->theSystem->*res)( residue,
              this->theConfiguration->getConf(0),
              this->theConfiguration->getConf(1),
              this->theConfiguration->getTime( )
            );

}

/**
 * Function for NLSolver jacobian computation.
 * @param jacobian Tangent matrix.
 * @param q_current Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::iterationJacobian( lmx::Matrix<T>& jacobian, lmx::Vector<T>& q_current )
{
  (this->theSystem->*jac)( jacobian,
              this->theConfiguration->getConf(0),
              static_cast< IntegratorBaseImplicit<T>* >(this->theIntegrator)->getPartialQdot( ),
              this->theConfiguration->getTime( )
            );
}

/**
 * Function for NLSolver convergence evaluation.
 * @param q_current Configuration computed by the NLSolver.
 */
template <typename Sys, typename T>
    bool DiffProblemFirst<Sys,T>::iterationConvergence( lmx::Vector<T>& q_current )
{
  return (this->theSystem->*conv)( this->theConfiguration->getConf(0),
                                   this->theConfiguration->getConf(1),
                                   this->theConfiguration->getTime( )
                                 );
}


/**
 * Initialize solving function
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::initialize( )
{
  this->theConfiguration->setTime( this->to );
  this->theIntegrator->initialize( this->theConfiguration );
  if ( ! this->theIntegrator->isExplicit() ){
    if (solveInitialEquilibrium) // default TRUE
    (this->theSystem->*eval)( this->theConfiguration->getConf(0),
                              this->theConfiguration->setConf(1),
                              this->theConfiguration->getTime( )
                            );
    if(this->vervosity<2) theNLSolver.setInfo(0);
    theNLSolver.setInitialConfiguration( this->theConfiguration->getConf(0) );
    theNLSolver.setDeltaInResidue(  );
    theNLSolver.setSystem( *this );
    if( b_convergence ){
      theNLSolver.setConvergence( &DiffProblemFirst<Sys,T>::iterationConvergence );
    }
    theNLSolver.setResidue( &DiffProblemFirst<Sys,T>::iterationResidue ); // Also advances the integrator
    theNLSolver.setJacobian( &DiffProblemFirst<Sys,T>::iterationJacobian );
  } 
  this->writeStepFiles();

}

/**
 * Solve main function
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::solve( )
{
  this->initialize();
  int max = (int)( (this->tf - this->to) / this->stepSize );
  for ( int i=0; i<max; ++i)
    this->stepSolve();
}

/**
 * Solve only one step 
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::stepSolve( )
{
  if ( this->theIntegrator->isExplicit() )
    this->stepSolveExplicit();
  else this->stepSolveImplicit();
}

/**
 * Explicit time scheme solver.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::stepSolveExplicit( )
{
  (this->theSystem->*eval)( this->theConfiguration->getConf(0),
			    this->theConfiguration->setConf(1),
			    this->theConfiguration->getTime( )+this->stepSize
		    );
  this->theConfiguration->nextStep( this->stepSize );
  this->theIntegrator->advance( );
  if(this->b_steptriggered) (this->theSystem->*stepTriggered)( );
  this->writeStepFiles();
}

/**
 * Implicit time scheme solver.
 */
template <typename Sys, typename T>
    void DiffProblemFirst<Sys,T>::stepSolveImplicit( )
{
  this->theConfiguration->nextStep( this->stepSize );
  this->theIntegrator->advance( );
  theNLSolver.solve( 10000 );
  if(this->b_steptriggered) (this->theSystem->*stepTriggered)( );
  this->writeStepFiles();
}

}; // namespace lmx


#endif
