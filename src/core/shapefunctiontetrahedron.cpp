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
#include "shapefunctiontetrahedron.h"
#include "point.h"
#include "node.h"

namespace mknix
{

ShapeFunctionTetrahedron::ShapeFunctionTetrahedron()
    : ShapeFunction()
{
}


ShapeFunctionTetrahedron::ShapeFunctionTetrahedron( Point* gp_in )
    : ShapeFunction(gp_in)
{
    // rows: shape function + 3 first derivatives (x,y,z)
    // columns: 4 points
    this->phi.resize(4, 4);
}


ShapeFunctionTetrahedron::~ShapeFunctionTetrahedron()
{
}


void ShapeFunctionTetrahedron::calc()
{
    //////////////////////////////////////////////////////////////////
    //  PRELIMINARY CALCULATIONS:
    //   (From Zien and Taylor, pg 128)
    this->compute_abcd(1,2,3); //node 0
    this->compute_abcd(2,3,0); //node 1
    this->compute_abcd(3,0,1); //node 2
    this->compute_abcd(0,1,2); //node 3

    //////////////////////////////////////////////////////////////////
    // phi(0,i) = (a_i + b_i*x + c_i*y + d_i*z) / (6*J)
    for( int i=0; i<4; ++i )
    {
        phi.writeElement
        (
            ( a[i] + b[i]*gp->X + c[i]*gp->Y + d[i]*gp->Z )
            / ( 6*gp->jacobian )
            , 0, i
        );
    }
    //////////////////////////////////////////////////////////////////
    // FIRST DERIVATIVES:
    //////////////////////////////////////////////////////////////////
    for( int i=0; i<4; ++i )
    {
        phi.writeElement( b[i] / ( 6*gp->jacobian ), 1, i ); // ,x
        phi.writeElement( c[i] / ( 6*gp->jacobian ), 2, i ); // ,y
        phi.writeElement( d[i] / ( 6*gp->jacobian ), 3, i ); // ,z
    }

    //     cout << "jacobian = " << gp->jacobian << endl;
    //     cout << "gp = (" << gp->X << ", " << gp->Y << ", " << gp->Z << ");\n";
    //     cout << "a = (" << a[0] << ", " << a[1] << ", " << a[2] << ", " << a[3] << ");\n";
    //     cout << "b = (" << b[0] << ", " << b[1] << ", " << b[2] << ", " << b[3] << ");\n";
    //     cout << "c = (" << c[0] << ", " << c[1] << ", " << c[2] << ", " << c[3] << ");\n";
    //     cout << "d = (" << d[0] << ", " << d[1] << ", " << d[2] << ", " << d[3] << ");\n";
    //     cout << "phi = " << phi << endl;
}

void ShapeFunctionTetrahedron::compute_abcd( int i1, int i2, int i3 )
{
    //   (From Zien and Taylor, pg 128) does not work...
    // Trying from www.iue.tuwien.ac.at/phd/nentchev/node30.html
    static int sign[4] = {1, -1, 1, -1};

    cofe::TensorRank2< 3, double > det_a, det_bcd;

    det_a(0,0) = gp->supportNodes[i1]->getX();
    det_a(0,1) = gp->supportNodes[i1]->getY();
    det_a(0,2) = gp->supportNodes[i1]->getZ();
    det_a(1,0) = gp->supportNodes[i2]->getX();
    det_a(1,1) = gp->supportNodes[i2]->getY();
    det_a(1,2) = gp->supportNodes[i2]->getZ();
    det_a(2,0) = gp->supportNodes[i3]->getX();
    det_a(2,1) = gp->supportNodes[i3]->getY();
    det_a(2,2) = gp->supportNodes[i3]->getZ();

//     cout << det_a << endl;

    a.push_back(sign[i2]*det_a.determinant());

    int i;

    det_bcd = det_a;
    for(i=0; i<3; ++i)
    {
        det_bcd(i,0) = 1.;
    }
    b.push_back( -sign[i2]*(det_bcd.determinant()) );

    det_bcd = det_a;
    for(i=0; i<3; ++i)
    {
        det_bcd(i,0) = 1.;
        det_bcd(i,1) = det_a(i,0);
    }
    c.push_back( sign[i2]*(det_bcd.determinant()) );

    for(i=0; i<3; ++i)
    {
        det_bcd(i,0) = 1.;
        det_bcd(i,1) = det_a(i,0);
        det_bcd(i,2) = det_a(i,1);
    }
    d.push_back( -sign[i2]*(det_bcd.determinant()) );

}


}
