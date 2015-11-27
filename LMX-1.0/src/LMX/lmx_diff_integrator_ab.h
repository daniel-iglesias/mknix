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

#ifndef LMXINTEGRATOR_AB_H
#define LMXINTEGRATOR_AB_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_ab.h
      
      \brief IntegratorAB class implementation

      Implements ADAMS-BASHFORD integrator class for solving dynamic systems.

      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include"lmx_diff_integrator_base.h"

namespace lmx {

    /**
    \class IntegratorAB
    \brief Template class IntegratorAB.
    Adams-Bashford integrator implementation for ODE systems.
    
    @author Daniel Iglesias .
    */
template <class T> class IntegratorAB : public IntegratorBase<T>
{

public:

  /** Empty constructor. */
  IntegratorAB(){}

  IntegratorAB(int ord);

  /** Destructor. */
  ~IntegratorAB(){}

  void initialize( Configuration<T>* );

  /** Returns 1 (TRUE) if it is an explicit-scheme integrator. */
  bool isExplicit()
  { return 1; }

  void advance( );


  private:
    int order;
    T b[5][5];  /**< Array (of arrays) of method's coefficients.*/
    Configuration<T>* theConfiguration;

};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

  template <class T> IntegratorAB<T>::IntegratorAB(int ord) : order(ord)
      /**
       * Standard constructor.
       * @param ord Differential order of system.
       */
  {
        // order == 1 -> Euler Explicit
          b[0][0] = (1.); // b[0] = 1
        
        // order == 2
          b[1][0] = (1.5); // b[0]
          b[1][1] = (-0.5); // b[1]
        
        // order == 3
          b[2][0] = (23./12.); // b[0]
          b[2][1] = (-16./12.); // b[1]
          b[2][2] = (5./12.); // b[2]
        
        // order == 4
          b[3][0] = (55./24.); // b[0]
          b[3][1] = (-59./24.); // b[1]
          b[3][2] = (37./24.); // b[2]
          b[3][3] = (-9./24.); // b[3]
        
        // order == 5
          b[4][0] = (1901./720.); // b[0]
          b[4][1] = (-2774./720.); // b[1]
          b[4][2] = (2616./720.); // b[2]
          b[4][3] = (-1274./720.); // b[3]
          b[4][4] = (251./720.); // b[4]

  }


  template <class T>
      void IntegratorAB<T>::initialize( Configuration<T>* configuration_in )
    /**
     * Miscellaneus setup operations.
     * @param configuration_in
     */
  {
    theConfiguration = configuration_in;
    theConfiguration->setStoredSteps( 2, order+1, order+1 );
  }

  template <class T>
      void IntegratorAB<T>::advance( )
      /**
       * Advances to next time-step.
       */
  { int i, j;
    for ( i = theConfiguration->getDiffOrder()-1; i>=0; --i ){
      theConfiguration->setConf( i, theConfiguration->getConf(i, 1) );//  *q[i][0] = *q[i][1]; // q[n] = q[n-1] + ...
      if( theConfiguration->getTimeSize() > this->order)
        for ( j=0; j < this->order; ++j){
        theConfiguration->setConf( i,
                                    theConfiguration->getConf( i, 0 ) + theConfiguration->getLastStepSize()*b[order-1][j] * theConfiguration->getConf( i+1, j+1 ) );
                           //        *q[i][0] += h * b[order-1][j] * (*q[i+1][j+1]); // ... + b_i*f_{n-i}
        }
      // Si se est� arrancando no se puede aplicar el orden completo de integraci�n:
      else
      {
        for ( j=0; j < theConfiguration->getTimeSize()-1; j++){
          theConfiguration->setConf( i,
                                      theConfiguration->getConf( i, 0 ) + theConfiguration->getLastStepSize() *
                                      b[theConfiguration->getTimeSize()-2][j] * theConfiguration->getConf( i+1, j+1 ) );
                          //            *q[i][0] += h * (b[step][j] * (*q[i+1][j+1]));
        }
      }
//       cout << "Derivative order = " << i << ", q_t = " << theConfiguration->getConf( i, 0 ) << endl;
    }

  }



}; // namespace lmx


#endif
