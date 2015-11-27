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
#include "shapefunctiontriangle3D.h"
#include "point.h"
#include "node.h"
#include "simulation.h"

namespace mknix {

ShapeFunctionTriangleSigned::ShapeFunctionTriangleSigned()
    : ShapeFunction()
{
//   this->phi.resize(6, 3);
}


ShapeFunctionTriangleSigned::ShapeFunctionTriangleSigned( Point* gp_in )
    : ShapeFunction(gp_in)
{
    this->phi.resize(3, 3);
}


ShapeFunctionTriangleSigned::~ShapeFunctionTriangleSigned()
{
}


void ShapeFunctionTriangleSigned::calc()
{
    // based on the answer in http://answers.unity3d.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
    // check out also: http://www.had2know.com/academics/triangle-area-perimeter-angle-3-coordinates.html
    lmx::Vector<double> f(dim), p1(dim), p2(dim), p3(dim);
    f.writeElement(gp->getX(),0);
    f.writeElement(gp->getY(),1);
    f.writeElement(gp->getZ(),2);
    p1.writeElement(gp->supportNodes[0]->getX(), 0);
    p1.writeElement(gp->supportNodes[0]->getY(), 1);
    p1.writeElement(gp->supportNodes[0]->getZ(), 2);
    p2.writeElement(gp->supportNodes[1]->getX(), 0);
    p1.writeElement(gp->supportNodes[1]->getY(), 1);
    p1.writeElement(gp->supportNodes[1]->getZ(), 2);
    p3.writeElement(gp->supportNodes[2]->getX(), 0);
    p1.writeElement(gp->supportNodes[2]->getY(), 1);
    p1.writeElement(gp->supportNodes[2]->getZ(), 2);

    lmx::Vector<double> f1(dim), f2(dim), f3(dim); // calculate vectors from point f to vertices p1, p2 and p3:
    f1.subs( p1, f);
    f2.subs( p2, f);
    f3.subs( p3, f);
    lmx::Vector<double> va(dim), va1(dim), va2(dim), va3(dim);
    double a, a1, a2, a3;
    va.mult(p1-p2, p1-p3);
    va1.mult(f2, f3);
    va2.mult(f3, f1);
    va3.mult(f1, f2);
    lmx::Vector<double> vaa1(dim), vaa2(dim), vaa3(dim);
    a = va.norm2();
    a1 = std::copysign( va1.norm2()/a, va*va1 );
    a2 = std::copysign( va2.norm2()/a, va*va2 );
    a3 = std::copysign( va3.norm2()/a, va*va3 );

    phi.writeElement( a1, 0, 0 );
    phi.writeElement( a2, 0, 0 );
    phi.writeElement( a3, 0, 0 );
    //////////////////////////////////////////////////////////////////
    // FIRST DERIVATIVES:
    //////////////////////////////////////////////////////////////////
//   phi.writeElement(
//       ( gp->supportNodes[1]->gety() - gp->supportNodes[2]->gety() )
//       / ( 2*gp->jacobian ), 1, 0 );
//   phi.writeElement(
//       ( gp->supportNodes[2]->getx() - gp->supportNodes[1]->getx() )
//       / ( 2*gp->jacobian ), 2, 0 );
//   //////////////////////////////////////////////////////////////////
//   phi.writeElement(
//       ( gp->supportNodes[2]->gety() - gp->supportNodes[0]->gety() )
//       / ( 2*gp->jacobian ), 1, 1 );
//   phi.writeElement(
//       ( gp->supportNodes[0]->getx() - gp->supportNodes[2]->getx() )
//       / ( 2*gp->jacobian ), 2, 1 );
//   //////////////////////////////////////////////////////////////////
//   phi.writeElement(
//       ( gp->supportNodes[0]->gety() - gp->supportNodes[1]->gety() )
//       / ( 2*gp->jacobian ), 1, 2 );
//   phi.writeElement(
//       ( gp->supportNodes[1]->getx() - gp->supportNodes[0]->getx() )
//       / ( 2*gp->jacobian ), 2, 2 );

//   cout << "phi = " << phi << endl;
}


}
