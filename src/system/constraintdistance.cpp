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
#include "constraintdistance.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix {

ConstraintDistance::ConstraintDistance()
    : Constraint()
{
}

ConstraintDistance::ConstraintDistance( Node* a_in, Node* b_in, double& alpha_in, std::string& method_in )
    : Constraint(alpha_in, method_in)
{
    nodes.push_back( a_in );
    nodes.push_back( b_in );
    size_type total_support_nodes = a_in->getSupportSize(0) + b_in->getSupportSize(0);
//   ro =  std::pow( nodes[1]->getx() - nodes[0]->getx(), 2 )
//       + std::pow( nodes[1]->gety() - nodes[0]->gety(), 2 )
//       + std::pow( nodes[1]->getz() - nodes[0]->getz(), 2 );

    cout << "node A: " << nodes[0]->getNumber() << endl;
    cout << "node A support size: " << a_in->getSupportSize(0) << endl;
    cout << "node B: " << nodes[1]->getNumber() << endl;
    cout << "node B support size: " << b_in->getSupportSize(0) << endl;

    calcRo();

    this->stiffnessMatrix.resize(total_support_nodes*dim,total_support_nodes*dim);
    this->internalForces.resize(total_support_nodes*dim);
    this->lambda.resize(1);
    this->lambda[0]=0.0;
    this->phi.resize(1);
    this->phi_q.resize(1);
    this->phi_q[0].resize(total_support_nodes*dim);
    this->phi_qq.resize(1);
    this->phi_qq[0].resize(total_support_nodes*dim,total_support_nodes*dim);
}

ConstraintDistance::~ConstraintDistance()
{
}

void ConstraintDistance::calcRo() {
    ro = std::sqrt( std::pow( nodes[1]->getConf(0) - nodes[0]->getConf(0), 2 )
                    +std::pow( nodes[1]->getConf(1) - nodes[0]->getConf(1), 2 )
                    +std::pow( nodes[1]->getConf(2) - nodes[0]->getConf(2), 2 ) );
//     cout<< "\nro: " << ro
// 	<< " = (" <<  nodes[1]->getConf(0) << " - " << nodes[0]->getConf(0) << ")^2"
// 	<< " + (" <<  nodes[1]->getConf(1) << " - " << nodes[0]->getConf(1) << ")^2"
// 	<< " + (" <<  nodes[1]->getConf(2) << " - " << nodes[0]->getConf(2) << ")^2"
// 	<< endl;
}


void ConstraintDistance::calcPhi()
{
//   rt =  std::pow( nodes[1]->getx() - nodes[0]->getx(), 2 )
//       + std::pow( nodes[1]->gety() - nodes[0]->gety(), 2 )
//       + std::pow( nodes[1]->getz() - nodes[0]->getz(), 2 );
//
//   this->phi[0] = rt - ro ;

    rt = std::sqrt( std::pow( nodes[1]->getConf(0) - nodes[0]->getConf(0), 2 )
                    +std::pow( nodes[1]->getConf(1) - nodes[0]->getConf(1), 2 )
                    +std::pow( nodes[1]->getConf(2) - nodes[0]->getConf(2), 2 ) );

//     cout<< "\nrt: " << rt
// 	<< " = (" <<  nodes[1]->getConf(0) << " - " << nodes[0]->getConf(0) << ")^2"
// 	<< " + (" <<  nodes[1]->getConf(1) << " - " << nodes[0]->getConf(1) << ")^2"
// 	<< " + (" <<  nodes[1]->getConf(2) << " - " << nodes[0]->getConf(2) << ")^2"
// 	<< endl;
//     cout << "ro = " << ro << ", rt = " << rt << endl;
    this->phi[0] = rt*rt - ro*ro ;

//   cout<< "\tro: " << ro << endl;
//   cout<< "\trt: " << rt << endl;
//   cout<< "\tphi: " << phi[0] << endl;
//       << " = (" <<  nodes[1]->getx() << " - " << nodes[0]->getx() << ")^2"
//       << " + (" <<  nodes[1]->gety() << " - " << nodes[0]->gety() << ")^2"
//       << " + (" <<  nodes[1]->getz() << " - " << nodes[0]->getz() << ")^2"
//       << endl;

}

void ConstraintDistance::calcPhiq()
{
    this->phi_q[0](0) = -2.0*( nodes[1]->getConf(0) - nodes[0]->getConf(0) ) ;
    this->phi_q[0](1) = -2.0*( nodes[1]->getConf(1) - nodes[0]->getConf(1) ) ;
    if (dim == 3)
        this->phi_q[0](2) = -2.0*( nodes[1]->getConf(2) - nodes[0]->getConf(2) ) ;
    this->phi_q[0](dim) = +2.0*( nodes[1]->getConf(0) - nodes[0]->getConf(0) ) ;
    this->phi_q[0](dim+1) = +2.0*( nodes[1]->getConf(1) - nodes[0]->getConf(1) ) ;
    if (dim == 3)
        this->phi_q[0](dim+2) = +2.0*( nodes[1]->getConf(2) - nodes[0]->getConf(2) ) ;
}

void ConstraintDistance::calcPhiqq()
{
    this->phi_qq[0](0,0) =  2.0 ;
    this->phi_qq[0](1,1) =  2.0 ;
    this->phi_qq[0](dim,dim) =  2.0 ;
    this->phi_qq[0](dim+1,dim+1) =  2.0 ;
    this->phi_qq[0](0,dim) = -2.0 ;
    this->phi_qq[0](dim,0) = -2.0 ;
    this->phi_qq[0](1,dim+1) = -2.0 ;
    this->phi_qq[0](dim+1,1) = -2.0 ;
    if (dim == 3) {
        this->phi_qq[0](2,2) =  2.0 ;
        this->phi_qq[0](2,dim+2) = -2.0 ;
        this->phi_qq[0](dim+2,2) = -2.0 ;
        this->phi_qq[0](dim+2,dim+2) =  2.0 ;
    }
}

}
