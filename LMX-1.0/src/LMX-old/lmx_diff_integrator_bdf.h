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

#ifndef LMXINTEGRATOR_BDF_H
#define LMXINTEGRATOR_BDF_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_bdf.h
      
      \brief IntegratorBDF class implementation

      Implements BDF integrator class for solving dynamic systems.

      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include<cmath>

#include"lmx_diff_integrator_base_implicit.h"

namespace lmx {

    /**
    \class IntegratorBDF
    \brief Template class IntegratorBDF.
    Gear's BDF integrator implementation for ODE systems.
    
    @author Daniel Iglesias.
    */
template <class T> class IntegratorBDF : public IntegratorBaseImplicit<T>
{
  public:
  
    /** Empty constructor. */
    IntegratorBDF(){}
  
    /** Standard constructor. */
    IntegratorBDF(int ord);
  
    /** Destructor. */
    ~IntegratorBDF(){}
  
    /** Initialize integration function. */
    void initialize( Configuration<T>* );
  
    /** Advance to next time-step function. */
    void advance();
  
    /** Actualize with delta in actual time-step. */
    void integratorUpdate( lmx::Vector<T> delta );

    /** Calculates the factor \f$ \frac{\partial qdot_n}{\partial q_n} \f$. */
    double getPartialQdot( )
    {
      return (q->getTimeSize() > this->order) ?
          1./ ( q->getLastStepSize()*b[order-1] ) :
          1./ ( q->getLastStepSize()*b[q->getTimeSize()-2] );
    }

    /** Calculates the factor \f$ \frac{\partial qddot_n}{\partial q_n} \f$. */
    double getPartialQddot( )
    {
      return (q->getTimeSize() > this->order) ?
          1./ ( std::pow( q->getLastStepSize()*b[order-1],2 ) ) :
          1./ ( std::pow( q->getLastStepSize()*b[q->getTimeSize()-2],2 ) );
    }
  
  private:
    int order;
    T b[5];  /**< Array of method's RHS coefficients.*/
    T a[5][5];  /**< Array (of arrays) of method's coefficients.*/
    Configuration<T>* q;

};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

  template <class T> IntegratorBDF<T>::IntegratorBDF(int ord) : order(ord)
  {
        // order == 1 -> Euler Implicit
          b[0] = (1.); // b = 1
          a[0][0] = (1.); // a[0] = 1
        
        // order == 2
          b[1] = (2./3.); // b
          a[1][0] = (4./3.); // a[0]
          a[1][1] = (-1./3.); // a[1]
        
        // order == 3
          b[2] = (6./11.); // b
          a[2][0] = (18./11.); // a[0]
          a[2][1] = (-9./11.); // a[1]
          a[2][2] = (2./11.); // a[2]
        
        // order == 4
          b[3] = (12./25.); // b
          a[3][0] = (48./25.); // a[0]
          a[3][1] = (-36./25.); // a[1]
          a[3][2] = (16./25.); // a[2]
          a[3][3] = (-3./25.); // a[3]
        
        // order == 5
          b[4] = (60./137.); // b
          a[4][0] = (300./137.); // a[0]
          a[4][1] = (-300./137.); // a[1]
          a[4][2] = (200./137.); // a[2]
          a[4][3] = (-75./137.); // a[3]
          a[4][4] = (-12./137.); // a[4]
  }

  template <class T>
      void IntegratorBDF<T>::initialize( Configuration<T>* configuration_in )
  {
    q = configuration_in;
    q->setStoredSteps( order+1, order+1, 1 );
  }


  template <class T>
      void IntegratorBDF<T>::advance( )
  { 
    int i, j;
    if (q->getTimeSize() < order+1){ //orden de integrador > steps
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        q->setConf( i+1,( q->getConf(i,0) )
            * ( 1. / (q->getLastStepSize() * b[q->getTimeSize()-2])) );
        for ( j=0; j < q->getTimeSize()-1; ++j){
          q->setConf( i+1, q->getConf(i+1,0) - ( a[q->getTimeSize()-2][j]* q->getConf(i,j+1)  )
              * ( 1. / (q->getLastStepSize() * b[q->getTimeSize()-2])) );
        }
      }
    }
     else{//se aplica el orden real de integrador
       for ( i = 0; i < q->getDiffOrder(); ++i ){
//          q->setConf( i+1,0).fillIdentity(0);
         q->setConf( i+1,( q->getConf(i,0) )
           * ( 1. / (q->getLastStepSize() * b[order-1])) );
         for ( j=0; j < order; ++j){
           q->setConf( i+1, q->getConf(i+1,0) - ( a[order-1][j]* q->getConf(i,j+1)  )
               * ( 1. / (q->getLastStepSize() * b[order-1])) );
         }
       }
     } 
  }

  template <class T>
      void IntegratorBDF<T>::integratorUpdate( lmx::Vector<T> delta )
  {
    int i;
    if (q->getTimeSize() < order+1){ //orden de integrador > steps
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        q->setConf( i+1) +=  delta
            * std::pow( 1. / (q->getLastStepSize() * b[q->getTimeSize()-2]), i+1);
      }
    }
    else{//se aplica el orden real de integrador
      for ( i = 0; i < q->getDiffOrder(); ++i ){
        q->setConf( i+1) += delta
            * std::pow( 1. / (q->getLastStepSize() * b[order-1]), i+1);
      }
    }    
    
    q->setConf( 0, q->getConf(0, 0) + delta );
  }

}; // namespace lmx


#endif

