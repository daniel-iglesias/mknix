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

#ifndef LMXINTEGRATOR_NEWMARK_H
#define LMXINTEGRATOR_NEWMARK_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_newmark.h
      
      \brief IntegratorNEWMARK class implementation

      Implements Beta-Newmark family integrators for solving dynamic systems.

      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

// #include<cmath>

#include"lmx_diff_integrator_base_implicit.h"

namespace lmx {

    /**
  \class IntegratorNEWMARK
  \brief Template class IntegratorNEWMARK.
  Gear's BDF integrator implementation for ODE systems.
    
  @author Daniel Iglesias.
     */
  template <class T> class IntegratorNEWMARK : public IntegratorBaseImplicit<T>
  {
    public:
  
      /** Empty constructor. */
      IntegratorNEWMARK(){}
  
      /** Standard constructor. */
      IntegratorNEWMARK(double beta, double gamma);
  
      /** Destructor. */
      ~IntegratorNEWMARK(){}
  
      /** Initialize integration function. */
      void initialize( Configuration<T>* );
  
      /** Advance to next time-step function. */
      void advance();
  
      /** Actualize with delta in actual time-step. */
      void integratorUpdate( lmx::Vector<T> delta );

      /** Calculates the factor \f$ \frac{\partial qdot_n}{\partial q_n} \f$. */
      double getPartialQdot( )
      {
        return gamma / ( q->getLastStepSize()*beta );
      }

      /** Calculates the factor \f$ \frac{\partial qddot_n}{\partial q_n} \f$. */
      double getPartialQddot( )
      {
        return 1./ ( std::pow( q->getLastStepSize(), 2. )*beta );
      }
  
    private:
      double beta, gamma;
      Configuration<T>* q;

  };

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

  template <class T> IntegratorNEWMARK<T>::
      IntegratorNEWMARK( double beta_in,
                         double gamma_in)
  : beta(beta_in), gamma(gamma_in)
  {
  }

  template <class T>
      void IntegratorNEWMARK<T>::initialize( Configuration<T>* configuration_in )
  {
    q = configuration_in;
    q->setStoredSteps( 2, 2, 2 );
  }


  template <class T>
      void IntegratorNEWMARK<T>::advance( )
  {
    q->setConf( 1,
                q->getConf( 1, 1 ) +
                    (T)(1. - gamma) * q->getLastStepSize()
                    * q->getConf( 2 , 1 ) );
    q->setConf( 0,
                q->getConf( 0, 1 ) +
                    (T)q->getLastStepSize() * q->getConf( 1, 1 ) +
                    (T)(0.5-beta) * std::pow( q->getLastStepSize(), 2 )
                    * q->getConf( 2, 1 ) );
    q->setConf( 2 ).fillIdentity( 0 );
  }

  template <class T>
      void IntegratorNEWMARK<T>::integratorUpdate( lmx::Vector<T> delta )
  {
    q->setConf( 0 ) += delta;
    q->setConf( 1 ) += delta * ( gamma/(beta*q->getLastStepSize() ) );
    q->setConf( 2 ) += delta * ( 1. / ( beta * std::pow(q->getLastStepSize(), 2 ) ) );
  }

}; // namespace lmx


#endif
