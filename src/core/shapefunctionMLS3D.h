/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of Nemesis.                                             *
 *                                                                            *
 *  Nemesis is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  Nemesis is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with Nemesis.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#ifndef SHAPEFUNCTIONMLS_H
#define SHAPEFUNCTIONMLS_H

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file shapefunctionMLS.h

  \brief Function for aproximation of basic variables by Moving Least
          Squares fit.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "shapefunction.h"

namespace mknix {

/**
@author Daniel Iglesias
*/
class ShapeFunctionMLS : public ShapeFunction {

private:
//    int dim; /**< Dimension of problem's space.*/
    int nn; /**< Number of support nodes */
    int mm; /**< Order of monomials(mm). */
    int m; /**< Size of the monomial vectors(depends on mm). */
    int weightType; /**< Selector of radial basis functions type.
                      * Can be == 0 for Gaussian (EXP)
                      *        == 1 Cubic spline
                      *        == 2 Quartic spline              */
    double s_max; /**< function's parameters. */
    lmx::DenseMatrix<double> w; /**< Weights and derivatives */
    VectorX<double> p; /**< Shifted polynomials of interest point */
    lmx::DenseMatrix<double> P; /**< Shifted Polynomials of support points */
    std::vector< lmx::DenseMatrix<double> > A; /**< vector of moment matrices*/
    std::vector< std::vector< VectorX<double> > > B; /** Container of B_I */

public:
    ShapeFunctionMLS();

    ShapeFunctionMLS(int, int, int, double&, double&, Point* );

    ~ShapeFunctionMLS();

    void calc();

private:
    void computeWeights();
    void computeMomentMatrix();

    void computePhi(double, double, double z=0);


};

} //Namespace mknix

#endif
