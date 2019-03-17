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

#ifndef LINEAR_SOLVERS_H
#define LINEAR_SOLVERS_H

#include"lmx_linsolvers_system.h"

//////////////////////////////////////////// Doxygen file documentation entry:
    /*!
      \file lmx_linsolvers.h
      
      \brief Linear_solvers collection

      Implements linear methods for a typical "A*x = b" system

      \author Daniel Iglesias
      
    */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {

template <typename T>
  /**
   * \brief Automatized driver for solving a linear system.
   *
   * Tries to guess the solver that better fits the system. It's use is intended for fast developing of non-optimized applications. For better performance and control of application, use the LinearSystem class.
   *
   * @param A LHS Matrix.
   * @param x LHS Vector of unknowns.
   * @param b RHS Vector.
   */
  void solveLinear(Matrix<T>& A, Vector<T>& x, Vector<T>& b)
{
  int old_solver_type = getLinSolverType(); //< used to restore the setting at the end of function.

  if (getMatrixType() == 2 || getMatrixType() == 3){
#ifdef HAVE_GMM
      cout << "Hay superLU..." <<endl;
    if (A.rows() < 1000)
      setLinSolverType(1);
    else setLinSolverType(3);
#endif
  }
  else if (A.rows() < 1000){
    if (getMatrixType() == 1){
#ifdef HAVE_SUPERLU
      setLinSolverType(1);
#else
      setLinSolverType(2);
    cout << "::WARNING::" << endl
         << "::Auto selected a solver for symmetric matrices,"
         << "::if there is no convergence force the use of other solver." << endl
         << "::END WARNING::" << endl;
#endif
    }
    else setLinSolverType(0);
  }
  else{
    setLinSolverType(2);
    cout << "::WARNING::" << endl
         << "::Auto selected a solver for symmetric matrices,"
         << "::if there is no convergence force the use of other solver." << endl
         << "::END WARNING::" << endl;
  }

  LinearSystem<T> system(A,x,b);
  cout << "Matrix type: " << getMatrixType() << endl;
  cout << "Solver selected: " << getLinSolverType() << endl;
  system.solveYourself();

  setLinSolverType(old_solver_type);
}


}; // namespace lmx


#endif
