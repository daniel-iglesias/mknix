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

#ifndef LMXINTEGRATOR_BASE_IMPLICIT_H
#define LMXINTEGRATOR_BASE_IMPLICIT_H

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_base_implicit.h
      
      \brief IntegratorBaseImplicit abstract class implementation

      Implements the basic structure that will have the implicit integrators.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "lmx_diff_integrator_base.h"

namespace lmx {

    /**
    \class IntegratorBaseImplicit
    \brief Template abstract class IntegratorBaseImplicit.
    
    Implicit integrator basic structure.
    
    @author Daniel Iglesias .
    */
template <class T> class IntegratorBaseImplicit : public IntegratorBase<T>{

public:

  /** Empty constructor. */
  IntegratorBaseImplicit(){}

  /** Destructor. */
  virtual ~IntegratorBaseImplicit( ){}

  /** Returns 1 (TRUE) if it is an explicit-scheme integrator. */
  bool isExplicit()
  { return 0; }

  /** Calculates the factor \f$ \frac{\partial qdot_n}{\partial q_n} \f$. */
  virtual double getPartialQdot( ) = 0;

  /** Calculates the factor \f$ \frac{\partial qddot_n}{\partial q_n} \f$. */
  virtual double getPartialQddot( ) = 0;

  /** Actualizes variables applying an increment. */
  virtual void integratorUpdate( lmx::Vector<T> delta ) = 0;
};

}; // namespace lmx


#endif
