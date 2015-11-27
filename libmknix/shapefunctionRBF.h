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

#ifndef SHAPEFUNCTIONRBF_H
#define SHAPEFUNCTIONRBF_H

//////////////////////////////////////////// Doxygen file documentation entry:
/**
 * \file shapefunctionRBF.h
 *
 * \brief Function for interpolation of basic variables
 *  by means of Radial basis Functions.
 *
 * \author Daniel Iglesias
 *
 */
//////////////////////////////////////////// Doxygen file documentation (end)

#include "shapefunction.h"

namespace mknix {

/**
@author Daniel Iglesias
*/
class ShapeFunctionRBF : public ShapeFunction {

private:
    int nn, mm; /**< Number of radial basis functions(nn) and monomials(mm). */
    int rbfType; /**< Selector of radial basis functions type.
                  * Can be == 0 for Multi-Quadrics (MQ)
                  *        == 1 for Gaussian (EXP)
                  *        == 2 Thin plate spline (TPS)
                  *        == 10 Wu-C2                  */
    double alpha_c, d_c, q; /**< function's parameters. */
    lmx::DenseMatrix<double> g_o;

public:
    ShapeFunctionRBF();

    ShapeFunctionRBF(int, int, int,
                     double&, double&, double&, Point* );

    ~ShapeFunctionRBF();

    void calc();

private:
    void computeMomentMatrix();

    void computePhi(double, double, double);


};

} //Namespace mknix

#endif
