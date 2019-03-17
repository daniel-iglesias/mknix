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

#ifndef LMXINTEGRATOR_CENTRALDIFF_H
#define LMXINTEGRATOR_CENTRALDIFF_H


//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_centraldiff.h
      
      \brief IntegratorCentralDiff class implementation

      Implements central differences explicit integrator class for solving dynamic systems.

      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include"lmx_diff_integrator_base.h"

namespace lmx {

    /**
    \class IntegratorCentralDifference
    \brief Template class IntegratorCentralDifference.
    Central difference integrator implementation for ODE systems.
    
    @author Daniel Iglesias.
    */
template <class T> class IntegratorCentralDifference : public IntegratorBase<T>
{
  public:

    /** Empty constructor. */
    IntegratorCentralDifference() : firstIteration(1)
    {}

    /** Destructor. */
    ~IntegratorCentralDifference(){}

    void initialize( Configuration<T>* );

    /** Returns 1 (TRUE) if it is an explicit-scheme integrator. */
    bool isExplicit()
    { return 1; }

    /** Advance to next time-step function. */
    void advance( );

  private:
    Configuration<T>* q;
    bool firstIteration;

};

}; // namespace lmx

/////////////////////////////// Implementation of the methods defined previously

namespace lmx {

  template <class T>
      void IntegratorCentralDifference<T>::initialize( Configuration<T>* configuration_in )
  {
    q = configuration_in;
    q->setStoredSteps( 2, 2, 2 );
  }

  template <class T>
      void IntegratorCentralDifference<T>::advance( )
  {
    if( q->getDiffOrder() == 2 ){
      if ( firstIteration )
        q->setConf( 1 ) += (T)(q->getLastStepSize() / 2) * q->getConf(2,1);
      else
        q->setConf( 1 ) += (T)q->getLastStepSize() * q->getConf(2,1);
      q->setConf( 0 ) += (T)q->getLastStepSize() * q->getConf(1,0);
//       cout << "Derivative order = 1" << ", q_t = " << q->getConf( 1, 0 ) << endl;
//       cout << "Derivative order = 0" << ", q_t = " << q->getConf( 0, 0 ) << endl;
    }
    else{
      std::stringstream message;
      message << "Differential system must be a second order to apply this method." << endl;
      LMX_THROW(failure_error, message.str() );
    }
  }

}; // namespace lmx


#endif
