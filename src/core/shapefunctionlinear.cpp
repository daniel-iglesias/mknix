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

#include "shapefunctionlinear.h"
#include "point.h"
#include "node.h"

namespace mknix
{

ShapeFunctionLinear::ShapeFunctionLinear()
    : ShapeFunction()
{
//   this->phi.resize(6, 3);
}


ShapeFunctionLinear::ShapeFunctionLinear( Point* gp_in )
    : ShapeFunction(gp_in)
{
    this->phi.resize(4, 2); // deriv 0 and dx dy for two nodes support
}


ShapeFunctionLinear::~ShapeFunctionLinear()
{
}


void ShapeFunctionLinear::calc()
{
    // Signed calculation as expressed in http://stackoverflow.com/questions/552916/how-to-find-sign-of-directed-distance
    double distance_x1x0 = gp->supportNodes[1]->distance( *(gp->supportNodes[0]) ); // should use the jacobian

    cout << "DISTANCE_X1X0 = " << distance_x1x0 << endl;

    phi.writeElement( - ( ( gp->getX() - gp->supportNodes[1]->getX() ) *  ( gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX() )
                          +( gp->getY() - gp->supportNodes[1]->getY() ) *  ( gp->supportNodes[1]->getY() - gp->supportNodes[0]->getY() )
                          +( gp->getZ() - gp->supportNodes[1]->getZ() ) *  ( gp->supportNodes[1]->getZ() - gp->supportNodes[0]->getZ() )
                        ) / (distance_x1x0 * distance_x1x0)
                      , 0, 0);
    phi.writeElement(  (  ( gp->getX() - gp->supportNodes[0]->getX() ) *  ( gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX() )
                          +( gp->getY() - gp->supportNodes[0]->getY() ) *  ( gp->supportNodes[1]->getY() - gp->supportNodes[0]->getY() )
                          +( gp->getZ() - gp->supportNodes[0]->getZ() ) *  ( gp->supportNodes[1]->getZ() - gp->supportNodes[0]->getZ() )
                       ) / (distance_x1x0 * distance_x1x0)
                       , 0, 1);

    cout << "support nodes = ("
         << gp->supportNodes[0]->getX() << " ," << gp->supportNodes[0]->getY() << ")("
         << gp->supportNodes[1]->getX() << " ," << gp->supportNodes[1]->getY() << ")" << endl;

    cout << "PHI-1D (" << gp->getX() <<", " << gp->getY() <<" ) = " << phi(0,0) << ", " << phi(0,1) << endl;
    //////////////////////////////////////////////////////////////////
    // FIRST DERIVATIVES:
    //////////////////////////////////////////////////////////////////
    phi.writeElement( -1./ distance_x1x0, 1, 0 );
    phi.writeElement( -1./ distance_x1x0, 2, 0 );
    phi.writeElement( -1./ distance_x1x0, 3, 0 );
    //////////////////////////////////////////////////////////////////
    phi.writeElement( 1. / distance_x1x0, 1, 1 );
    phi.writeElement( 1. / distance_x1x0, 2, 1 );
    phi.writeElement( 1. / distance_x1x0, 3, 1 );
//   cout << "phi = " << phi << endl;
}


}
