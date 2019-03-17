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

#ifndef LMX_H
#define LMX_H

//////////////////////////////////////////// Doxygen file documentation entry:
    /**
     * \file lmx.h
     *
     * \brief This is the main file that has to be included from an extern program so that the basic features of the library can be used.
     *
     * This basic facilities include all the linear algebra operations over Vectors and Matrices, Tensor algebra and linear systems solvers.
     *
     * \author Daniel Iglesias
     * 
     */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "lmx_linsolvers.h"
#include "cofe_fmc.h"
#include "lmx_except.h"
#include "lmx_base_stopwatch.h"
#include "lmx_base_selector.h"

/** \namespace lmx{}
 *
 *  \brief Linked MatriX methods library uses lmx namespace
 *
 *  This is the namespace that contains LMX-lite.
 *
 */

namespace lmx{

}

#endif
