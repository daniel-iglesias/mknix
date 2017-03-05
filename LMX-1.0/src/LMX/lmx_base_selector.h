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
#ifndef LMXSELECTOR_H
#define LMXSELECTOR_H

#include <iostream>
#include <sstream>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file lmx_base_selector.h

  \brief This file contains set-get functions for switching between the Matrix, Vector and Linear Solvers types.

  Selector class is not usually instantiated. Instead, its friend functions are used by Matrix, Vector and LinearSystem objects for getting the type that is being used.

  \author Daniel Iglesias

*/
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {

/** Function that changes the type of matrix container that will be used.
 */
inline int setMatrixType(int type)
{
    static int matrix_type = 0;
    static bool type_locked = 0;
    if (type < 0) {
        type_locked = 1;
        return matrix_type;
    }
    else if (type_locked == 0) {
        return matrix_type = type;
    } else {
        std::stringstream message;
        message << "\nERROR: MatrixType can't be set after a Matrix object is created." << endl;
        LMX_THROW(failure_error, message.str());
    }
}

/** Function that changes the type of vector container that will be used.
 */
inline int setVectorType(int type)
{
    static int vector_type = 0;
    static bool type_locked = 0;
    if (type < 0) {
        type_locked = 1;
        return vector_type;
    }
    else if (type_locked == 0) {
        return vector_type = type;
    } else {
        std::stringstream message;
        message << "\nERROR: VectorType can't be set after a Vector object is created." << endl;
        LMX_THROW(failure_error, message.str());
    }
}

/** Function that changes the type of linear solver that will be used.
 */
inline int setLinSolverType(int type)
{
    static int lin_solver_type = 0;
    if (type < 0) {
        return lin_solver_type;
    } else { return lin_solver_type = type; }
}

/** Function reads the type of Matrix container that is used.
 */
inline int getMatrixType() { return setMatrixType(-1); }

/** Function reads the type of Vector container that is used.
 */
inline int getVectorType() { return setVectorType(-1); }

/** Function reads the type of LinearSystem that is used.
 */
inline int getLinSolverType() { return setLinSolverType(-1); }


constexpr int nl_solver_type = 0; /**< This variable switches between different types of non-linear solvers. */


} // namespace lmx 


#endif
