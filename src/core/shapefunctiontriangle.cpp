/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                             *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#include "shapefunctiontriangle.h"
#include "point.h"
#include "node.h"

#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor for ShapeFunctionTriangle.
 */
ShapeFunctionTriangle::ShapeFunctionTriangle()
    : ShapeFunction()
{
//   this->phi.resize(6, 3);
}


/**
 * @brief Constructs a triangular FEM shape function for the given Gauss point.
 * @param gp_in Pointer to the evaluation Point (expects exactly 3 support nodes).
 */
ShapeFunctionTriangle::ShapeFunctionTriangle( Point* gp_in )
    : ShapeFunction(gp_in)
{
    this->phi.resize(3, 3);
}


/**
 * @brief Destructor for ShapeFunctionTriangle.
 */
ShapeFunctionTriangle::~ShapeFunctionTriangle()
{
}


/**
 * @brief Computes linear triangular (area-coordinate) FEM shape functions and first spatial derivatives
 *        using the standard Zienkiewicz & Taylor formula.
 */
void ShapeFunctionTriangle::calc()
{
    // phi(0,0) = ( x1*y2 - x2*y1 + (y1-y2)*x_gp + (x2-x1)*y_gp ) / (2*J)
    phi.writeElement
    (
        (gp->supportNodes[1]->getX() * gp->supportNodes[2]->getY()
         - gp->supportNodes[1]->getY() * gp->supportNodes[2]->getX()
         + ( gp->supportNodes[1]->getY() - gp->supportNodes[2]->getY() ) * gp->X
         + ( gp->supportNodes[2]->getX() - gp->supportNodes[1]->getX() ) * gp->Y
        ) / ( 2*gp->jacobian )
        , 0, 0
    );
    //////////////////////////////////////////////////////////////////
    // phi(0,1) = ( x2*y0 - x0*y2 + (y2-y0)*x_gp + (x0-x2)*y_gp ) / (2*J)
    phi.writeElement
    (
        (gp->supportNodes[2]->getX() * gp->supportNodes[0]->getY()
         - gp->supportNodes[2]->getY() * gp->supportNodes[0]->getX()
         + ( gp->supportNodes[2]->getY() - gp->supportNodes[0]->getY() ) * gp->X
         + ( gp->supportNodes[0]->getX() - gp->supportNodes[2]->getX() ) * gp->Y
        ) / ( 2*gp->jacobian )
        , 0, 1
    );
    //////////////////////////////////////////////////////////////////
    // phi(0,2) = ( x0*y1 - x1*y0 + (y0-y1)*x_gp + (x1-x0)*y_gp ) / (2*J)
    phi.writeElement
    (
        (gp->supportNodes[0]->getX() * gp->supportNodes[1]->getY()
         - gp->supportNodes[0]->getY() * gp->supportNodes[1]->getX()
         + ( gp->supportNodes[0]->getY() - gp->supportNodes[1]->getY() ) * gp->X
         + ( gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX() ) * gp->Y
        ) / ( 2*gp->jacobian )
        , 0, 2
    );
    //////////////////////////////////////////////////////////////////
    // FIRST DERIVATIVES:
    //////////////////////////////////////////////////////////////////
    phi.writeElement(
        ( gp->supportNodes[1]->getY() - gp->supportNodes[2]->getY() )
        / ( 2*gp->jacobian ), 1, 0 );
    phi.writeElement(
        ( gp->supportNodes[2]->getX() - gp->supportNodes[1]->getX() )
        / ( 2*gp->jacobian ), 2, 0 );
    //////////////////////////////////////////////////////////////////
    phi.writeElement(
        ( gp->supportNodes[2]->getY() - gp->supportNodes[0]->getY() )
        / ( 2*gp->jacobian ), 1, 1 );
    phi.writeElement(
        ( gp->supportNodes[0]->getX() - gp->supportNodes[2]->getX() )
        / ( 2*gp->jacobian ), 2, 1 );
    //////////////////////////////////////////////////////////////////
    phi.writeElement(
        ( gp->supportNodes[0]->getY() - gp->supportNodes[1]->getY() )
        / ( 2*gp->jacobian ), 1, 2 );
    phi.writeElement(
        ( gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX() )
        / ( 2*gp->jacobian ), 2, 2 );

//   cout << "phi = " << phi << endl;
}


}
