/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   http://code.google.com/p/mknix                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "shapefunctionlinear-x.h"
#include "point.h"
#include "node.h"

namespace mknix {

ShapeFunctionLinearX::ShapeFunctionLinearX()
    : ShapeFunction()
{
//   this->phi.resize(6, 3);
}


ShapeFunctionLinearX::ShapeFunctionLinearX( Point* gp_in )
    : ShapeFunction(gp_in)
{
    this->phi.resize(4, 2); // deriv 0 and dx dy for two nodes support
}


ShapeFunctionLinearX::~ShapeFunctionLinearX()
{
}


void ShapeFunctionLinearX::calc()
{   // Signed calculation using only the x components of the points
    double distance_x1x0 = gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX(); // should use the jacobian

    cout << "DISTANCE_X1X0 = " << distance_x1x0 << endl;

    phi.writeElement( - ( ( gp->getX() - gp->supportNodes[1]->getX() ) *  ( gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX() )
                        ) / (distance_x1x0 * distance_x1x0)
                      , 0, 0);
    phi.writeElement(  (  ( gp->getX() - gp->supportNodes[0]->getX() ) *  ( gp->supportNodes[1]->getX() - gp->supportNodes[0]->getX() )
                       ) / (distance_x1x0 * distance_x1x0)
                       , 0, 1);

    cout << "support nodes = ("
         << gp->supportNodes[0]->getX() << ")("
         << gp->supportNodes[1]->getX() << ")" << endl;

    cout << "PHI-1D-X (" << gp->getX() <<" ) = " << phi(0,0) << ", " << phi(0,1) << endl;
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
