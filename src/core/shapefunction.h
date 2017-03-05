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

#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H

#include "LMX/lmx.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file shapefunction.h

  \brief Function for interpolation of basic variables.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

class Point;

/**
@author Daniel Iglesias
*/
class ShapeFunction {

protected:
    size_t dim;
    lmx::DenseMatrix<double> phi;
    Point* gp;

public:
    ShapeFunction();

    ShapeFunction( const ShapeFunction* );

    ShapeFunction( Point* );

    virtual ~ShapeFunction();

    virtual void calc() = 0;

    virtual double getPhi(size_t i, size_t j) // i = derivative order, j = node
    {
        return phi.readElement(i, j);
    }

    virtual void setPhi(double value_in, size_t i, size_t j) // i = derivative order, j = node
    {
        phi.writeElement(value_in, i, j);
    }

    virtual void outputValues();

    virtual void gnuplotOut();


};

} //Namespace mknix

#endif
