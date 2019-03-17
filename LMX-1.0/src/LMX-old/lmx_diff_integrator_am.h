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

#ifndef LMXINTEGRATOR_AM_H
#define LMXINTEGRATOR_AM_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_am.h
      
      \brief IntegratorAM class implementation

      Implements ADAMS-MOULTON integrator class for solving dynamic systems.

      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include<cmath>

#include"lmx_diff_integrator_base_implicit.h"

namespace lmx {

    /**
    \class IntegratorAM
    \brief Template class IntegratorAM.
    Adams-Moulton integrator implementation for ODE systems.
    
    @author Daniel Iglesias.
    */
template <class T> class IntegratorAM : public IntegratorBaseImplicit<T>
{
private:
  int order;
  T b[5][5];  /**< Array (of arrays) of method's coefficients.*/
  Configuration<T>* q;

public:

  /** Empty constructor. */
  IntegratorAM(){}

  /** Standard constructor. */
  IntegratorAM(int ord);

  /** Destructor. */
  ~IntegratorAM(){}

  /** Initialize integration function. */
  void initialize( Configuration<T>* );

  /** Advance to next time-step function. */
  void advance();

  /** Actualize with delta in actual time-step. */
  void integratorUpdate( lmx::Vector<T> delta );

  /** Calculates the factor \f$ \frac{\partial qdot_n}{\partial q_n} \f$. */
  double getPartialQdot( )
  { return (q->getTimeSize() > this->order-1) ?
      1./ ( q->getLastStepSize()*b[order-1][0] ) :
      1./ ( q->getLastStepSize()*b[q->getTimeSize()-2][0] ); }

  /** Calculates the factor \f$ \frac{\partial qddot_n}{\partial q_n} \f$. */
  double getPartialQddot( )
  { return (q->getTimeSize() > this->order-1) ?
      1./ ( std::pow( q->getLastStepSize()*b[order-1][0], 2 ) ) :
      1./ ( std::pow( q->getLastStepSize()*b[q->getTimeSize()-2][0], 2 ) ); }
};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

  template <class T> IntegratorAM<T>::IntegratorAM(int ord) : order(ord)
  {
        // order == 1 -> Euler Implicit
          b[0][0] = (1.); // b[0] = 1
        
        // order == 2 -> Trapezoidal rule
          b[1][0] = (0.5); // b[0] = 1
          b[1][1] = (0.5); // b[1] = 1
        
        // order == 3
          b[2][0] = (5./12.); // b[0]
          b[2][1] = (8./12.); // b[1]
          b[2][2] = (-1./12.); // b[2]
        
        // order == 4
          b[3][0] = (9./24.); // b[0]
          b[3][1] = (19./24.); // b[1]
          b[3][2] = (-5./24.); // b[2]
          b[3][3] = (1./24.); // b[3]
        
        // order == 5
          b[4][0] = (251./720.); // b[0]
          b[4][1] = (646./720.); // b[1]
          b[4][2] = (-264./720.); // b[2]
          b[4][3] = (106./720.); // b[3]
          b[4][4] = (-19./720.); // b[4]
  }

  template <class T>
      void IntegratorAM<T>::initialize( Configuration<T>* configuration_in )
  {
#undef max

    q = configuration_in;
    q->setStoredSteps( 2, std::max(order,2), order );
  }


  template <class T>
      void IntegratorAM<T>::advance( )
  {
    int i, j;
    if (q->getTimeSize() < order/*+1*/){ //orden de integrador-1 > steps
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        // f{n} = (q{n} - q{n-1}) / ( h*b_0 ) + ...
        q->setConf( i+1,( q->getConf(i,0) - q->getConf(i,1) )
            * ( 1. / (q->getLastStepSize() * b[q->getTimeSize()-2][0])) );
        for ( j=1; j < q->getTimeSize()-1; ++j){
        // ... - b_i/b_0 * f{n-i}
          q->setConf( i+1, q->getConf(i+1,0) - ( b[q->getTimeSize()-2][j]* q->getConf(i+1,j)  )
              * ( 1. / ( b[q->getTimeSize()-2][0])) );
        }
      }
    }
    else{//se aplica el orden real de integrador
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        // f{n} = (q{n} - q{n-1}) / ( h*b_0 ) + ...
        q->setConf( i+1,( q->getConf(i,0) - q->getConf(i,1) )
            * ( 1. / (q->getLastStepSize() * b[order-1][0])) );
        for ( j=1; j < order; ++j){
        // ... - b_i/b_0 * f{n-i}
          q->setConf( i+1, q->getConf(i+1,0) - ( b[order-1][j]* q->getConf(i+1,j)  )
              * ( 1. / ( b[order-1][0])) );
        }
      }
    }
  }

  template <class T>
      void IntegratorAM<T>::integratorUpdate( lmx::Vector<T> delta )
  {
    int i;
    if (q->getTimeSize() < order/*+1*/){ //orden de integrador-1 > steps
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        q->setConf( i+1 ) += delta
          * ( 1. /
             std::pow(q->getLastStepSize() * b[q->getTimeSize()-2][0], i+1) );
      }
    }
    else{
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        q->setConf( i+1 ) += delta
          * ( 1. /
            std::pow(q->getLastStepSize() * b[order-1][0], i+1) );
      }
    }
    q->setConf( 0, q->getConf(0, 0) + delta );
  }

}; // namespace lmx


#endif
