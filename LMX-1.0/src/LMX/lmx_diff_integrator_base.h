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

#ifndef LMXINTEGRATOR_BASE_H
#define LMXINTEGRATOR_BASE_H

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_diff_integrator_base.h
      
      \brief IntegratorBase abstract class implementation

      Implements the basic structure that will have the integrators that solve initial value problems.

      \author Daniel Iglesias

    */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "lmx_diff_configuration.h"

namespace lmx {

    /**
    \class IntegratorBase
    \brief Template abstract class IntegratorBase.
    
    Integrator basic structure.
    
    @author Daniel Iglesias .
    */
template <class T> class IntegratorBase{

public:

  /** Empty constructor. */
  IntegratorBase(){}

  /** Destructor. */
  virtual ~IntegratorBase( ){}

  /** Setup the object. */
  virtual void initialize( Configuration<T>* ) = 0;

  /** Returns 1 (TRUE) if it is an explicit-scheme integrator. */
  virtual bool isExplicit() = 0;

  /** Actualization of the variables with non-actual configuration terms. */
  virtual void advance( ) = 0;
};

}; // namespace lmx


#endif
