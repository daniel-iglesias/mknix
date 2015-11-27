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
#include "constraintfixedcoordinates.h"
#include "simulation.h"
#include "node.h"

namespace mknix {

ConstraintFixedCoordinates::ConstraintFixedCoordinates()
    : Constraint()
{
}


ConstraintFixedCoordinates::ConstraintFixedCoordinates( Node* a_in, Node* b_in, double& alpha_in, std::string& method_in )
    : Constraint(alpha_in, method_in)
{
    nodes.push_back( a_in );
    nodes.push_back( b_in );
    size_type total_support_nodes = a_in->getSupportSize(0) + b_in->getSupportSize(0);

//   cout << "node A: " << nodes[0]->getNumber() << endl;
//   cout << "node A support size: " << a_in->getSupportSize(0) << endl;
//   cout << "node B: " << nodes[1]->getNumber() << endl;
//   cout << "node B support size: " << b_in->getSupportSize(0) << endl;

    rxo = nodes[1]->getConf(0) - nodes[0]->getConf(0) ;
    ryo = nodes[1]->getConf(1) - nodes[0]->getConf(1) ;
    rzo = nodes[1]->getConf(2) - nodes[0]->getConf(2) ;

    this->stiffnessMatrix.resize(total_support_nodes*dim,total_support_nodes*dim);
    this->internalForces.resize(total_support_nodes*dim);

    this->lambda.resize(dim);
    this->lambda[0]=0.0;
    this->lambda[1]=0.0;
    if (dim == 3)
        this->lambda[2]=0.0;

    this->phi.resize(dim);

    this->phi_q.resize(dim);
    this->phi_q[0].resize(total_support_nodes*dim);
    this->phi_q[1].resize(total_support_nodes*dim);
    if (dim == 3)
        this->phi_q[2].resize(total_support_nodes*dim);

    this->phi_qq.resize(dim);
    this->phi_qq[0].resize(total_support_nodes*dim,total_support_nodes*dim);
    this->phi_qq[1].resize(total_support_nodes*dim,total_support_nodes*dim);
    if (dim == 3)
        this->phi_qq[2].resize(total_support_nodes*dim,total_support_nodes*dim);
}


ConstraintFixedCoordinates::~ConstraintFixedCoordinates()
{
}

void ConstraintFixedCoordinates::calcPhi()
{
    rxt = nodes[1]->getConf(0) - nodes[0]->getConf(0) ;
    ryt = nodes[1]->getConf(1) - nodes[0]->getConf(1) ;
    if (dim == 3)
        rzt = nodes[1]->getConf(2) - nodes[0]->getConf(2);

    this->phi[0] = rxt - rxo ;
    this->phi[1] = ryt - ryo ;
    if (dim == 3)
        this->phi[2] = rzt - rzo ;

//  cout << endl << "rxt = " << nodes[1]->getConf(0) << " - " << nodes[0]->getConf(0)<< " = " << rxt << endl;
//  cout << endl << "Phi in fixedcoordinates = " << phi[0] << endl;
}

void ConstraintFixedCoordinates::calcPhiq()
{   size_type i,j;
    for(i=0; i<nodes[0]->getSupportSize(0); ++i) {
        this->phi_q[0](dim*i+0) = -nodes[0]->getShapeFunValue(0,i);
        this->phi_q[1](dim*i+1) = -nodes[0]->getShapeFunValue(0,i);
        if (dim == 3)
            this->phi_q[2](dim*i+2) = -nodes[0]->getShapeFunValue(0,i);
    }
    // Now i=supportNodes[0].size()
    for(j=0; j<nodes[1]->getSupportSize(0); ++j) {
        this->phi_q[0](dim*(i+j)) = nodes[1]->getShapeFunValue(0,j);
        this->phi_q[1](dim*(i+j)+1) = nodes[1]->getShapeFunValue(0,j);
        if (dim == 3)
            this->phi_q[2](dim*(i+j)+2) = nodes[1]->getShapeFunValue(0,j);
    }
//  this->phi_q[0](0) = -1.0;
//  this->phi_q[0](dim) = +1.0;
//
//  this->phi_q[1](1) = -1.0;
//  this->phi_q[1](dim+1) = +1.0;
//
//  if (dim == 3){
//    this->phi_q[2](2) = -1.0;
//    this->phi_q[2](dim+2) = +1.0;
//  }
}

void ConstraintFixedCoordinates::calcPhiqq()
{
}

}
