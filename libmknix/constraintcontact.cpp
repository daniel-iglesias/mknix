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
#include "constraintcontact.h"
#include "simulation.h"
#include "node.h"

namespace mknix {

ConstraintContact::ConstraintContact()
    : Constraint()
{
}

ConstraintContact::ConstraintContact( Node* q1_in, Node* q2_in, Node* p_in, double& alpha_in, std::string& method_in )
    : Constraint(alpha_in, method_in)
{
    nodes.push_back( q1_in );
    nodes.push_back( q2_in );
    nodes.push_back( p_in );

    this->stiffnessMatrix.resize(3*dim, 3*dim);
    this->internalForces.resize(3*dim);
    this->lambda.resize(1);
    this->lambda[0]=0.0;
    this->phi.resize(1);
    this->phi_q.resize(1);
    this->phi_q[0].resize(3*dim);
    this->phi_qq.resize(1);
    this->phi_qq[0].resize(3*dim, 3*dim);
    normal.resize(dim);
}

ConstraintContact::~ConstraintContact()
{
}

void ConstraintContact::calcPhi()
{
    normal[0] = +( nodes[1]->getConf(1) - nodes[0]->getConf(1) );
    normal[1] = -( nodes[1]->getConf(0) - nodes[0]->getConf(0) );
    double normal_len = std::sqrt( normal[0]*normal[0] + normal[1]*normal[1] );
    normal[0] /= normal_len;
    normal[1] /= normal_len;
    //gap = rt = ||(P-Q1) · normal|| = ||(P-Q2) · normal||
    //         = ||0.5*( ((P-Q1) + (P-Q2)) · normal)||
    rh = (nodes[2]->getConf(0) - nodes[1]->getConf(0) ) * normal[0]
         +(nodes[2]->getConf(1) - nodes[1]->getConf(1) ) * normal[1];
    rt = (nodes[2]->getConf(0) - nodes[0]->getConf(0) ) * normal[0]
         + (nodes[2]->getConf(1) - nodes[0]->getConf(1) ) * normal[1];

    this->phi[0] = rt*rt ;
//  if( rt < 0. ) cout << endl << "NEGATIVE GAP: " << rt << ", " << rh << endl;
}

void ConstraintContact::calcPhiq()
{
    if( rt < 0. )
    {
        this->phi_q[0](0) = -normal[0]*rt;
        this->phi_q[0](1) = -normal[1]*rt;
//    if (dim == 3)
//      this->phi_q[0](2) = -2.0*( nodes[1]->getz() - nodes[0]->getz() ) ;

        this->phi_q[0](dim) = -normal[0]*rt;
        this->phi_q[0](dim+1) = -normal[1]*rt;
//    if (dim == 3)
//      this->phi_q[0](2) = -2.0*( nodes[1]->getz() - nodes[0]->getz() ) ;

        this->phi_q[0](2*dim) = 2*normal[0]*rt;
        this->phi_q[0](2*dim+1) = 2*normal[1]*rt;
//    if (dim == 3)
//      this->phi_q[0](dim+2) = +2.0*( nodes[1]->getz() - nodes[0]->getz() ) ;
    }
    else if( rt > 0. )
    {
        //this->phi_q[0].fillIdentity( 0.0 )
        this->phi_q[0](0) = 0.0 ;
        this->phi_q[0](1) = 0.0 ;
        if (dim == 3)
            this->phi_q[0](2) = 0.0 ;
        this->phi_q[0](dim) = 0.0 ;
        this->phi_q[0](dim+1) = 0.0 ;
        if (dim == 3)
            this->phi_q[0](dim+2) = 0.0 ;
        this->phi_q[0](2*dim) = 0.0 ;
        this->phi_q[0](2*dim+1) = 0.0 ;
        if (dim == 3)
            this->phi_q[0](2*dim+2) = 0.0 ;
    }
}

void ConstraintContact::calcPhiqq()
{
    if( rt < 0. )
    {
        this->phi_qq[0](0,0) =  pow(normal[0],2)/2. ;
        this->phi_qq[0](0,1) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](0,dim) =  pow(normal[0],2)/2. ;
        this->phi_qq[0](0,dim+1) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](0,2*dim) =  -pow(normal[0],2) ;
        this->phi_qq[0](0,2*dim+1) =  -normal[0]*normal[1] ;

        this->phi_qq[0](1,0) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](1,1) =  pow(normal[1],2)/2. ;
        this->phi_qq[0](1,dim) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](1,dim+1) =  pow(normal[1],2)/2. ;
        this->phi_qq[0](1,2*dim) =  -normal[0]*normal[1] ;
        this->phi_qq[0](1,2*dim+1) =  -pow(normal[1],2) ;

        this->phi_qq[0](2,0) =  pow(normal[0],2)/2. ;
        this->phi_qq[0](2,1) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](2,dim) =  pow(normal[0],2)/2. ;
        this->phi_qq[0](2,dim+1) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](2,2*dim) =  -pow(normal[0],2) ;
        this->phi_qq[0](2,2*dim+1) =  -normal[0]*normal[1] ;

        this->phi_qq[0](3,0) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](3,1) =  pow(normal[1],2)/2. ;
        this->phi_qq[0](3,dim) =  normal[0]*normal[1]/2. ;
        this->phi_qq[0](3,dim+1) =  pow(normal[1],2)/2. ;
        this->phi_qq[0](3,2*dim) =  -normal[0]*normal[1] ;
        this->phi_qq[0](3,2*dim+1) =  -pow(normal[1],2) ;

        this->phi_qq[0](4,0) =  -pow(normal[0],2) ;
        this->phi_qq[0](4,1) =  -normal[0]*normal[1] ;
        this->phi_qq[0](4,dim) =  -pow(normal[0],2) ;
        this->phi_qq[0](4,dim+1) =  -normal[0]*normal[1] ;
        this->phi_qq[0](4,2*dim) =  2*pow(normal[0],2) ;
        this->phi_qq[0](4,2*dim+1) =  2*normal[0]*normal[1] ;

        this->phi_qq[0](5,0) =  -normal[0]*normal[1] ;
        this->phi_qq[0](5,1) =  -pow(normal[1],2) ;
        this->phi_qq[0](5,dim) =  -normal[0]*normal[1] ;
        this->phi_qq[0](5,dim+1) =  -pow(normal[1],2) ;
        this->phi_qq[0](5,2*dim) =  2*normal[0]*normal[1] ;
        this->phi_qq[0](5,2*dim+1) =  2*pow(normal[1],2) ;


        if (dim == 3) {
//      this->phi_qq[0](2,2) =  2.0 ;
//      this->phi_qq[0](dim+2,dim+2) =  2.0 ;
//      this->phi_qq[0](2,dim+2) = -2.0 ;
//      this->phi_qq[0](dim+2,2) = -2.0 ;
        }
    }
    else if( rt > 0. )
    {
        //this->phi_qq[0].fillIdentity( 0.0 )
        this->phi_qq[0](0,0) = 0.0 ;
        this->phi_qq[0](1,1) = 0.0 ;
        this->phi_qq[0](dim,dim) = 0.0 ;
        this->phi_qq[0](dim+1,dim+1) = 0.0 ;
        this->phi_qq[0](2*dim,2*dim) =  0.0 ;
        this->phi_qq[0](2*dim+1,2*dim+1) =  0.0 ;
        this->phi_qq[0](0,2*dim) = 0.0 ;
        this->phi_qq[0](2*dim,0) = 0.0 ;
        this->phi_qq[0](1,2*dim+1) = 0.0 ;
        this->phi_qq[0](2*dim+1,1) = 0.0 ;
        this->phi_qq[0](dim,2*dim) = 0.0 ;
        this->phi_qq[0](2*dim,dim) = 0.0 ;
        this->phi_qq[0](dim+1,2*dim+1) = 0.0 ;
        this->phi_qq[0](2*dim+1,dim+1) = 0.0 ;
        if (dim == 3) {
//      this->phi_qq[0](2,2) = 0.0 ;
//      this->phi_qq[0](2,dim+2) = 0.0 ;
//      this->phi_qq[0](dim+2,2) = 0.0 ;
//      this->phi_qq[0](dim+2,dim+2) = 0.0 ;
        }
    }
}

}
