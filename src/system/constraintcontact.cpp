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

#include "constraintcontact.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor.
 */
ConstraintContact::ConstraintContact()
    : Constraint()
{
}

/**
 * @brief Constructor for a contact constraint between a surface segment and a contact point.
 * @param q1_in     Pointer to the first node defining the contact surface.
 * @param q2_in     Pointer to the second node defining the contact surface.
 * @param p_in      Pointer to the contacting node.
 * @param alpha_in  Penalty parameter.
 * @param method_in Enforcement method keyword.
 */
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

/**
 * @brief Destructor.
 */
ConstraintContact::~ConstraintContact()
{
}

/**
 * @brief Evaluates the contact gap function phi using the surface normal and signed gap distance.
 */
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

/**
 * @brief Computes the gradient of phi with respect to the nodal coordinates.
 *        Returns zero when the contact gap is open (rt > 0).
 */
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

/**
 * @brief Computes the Hessian of phi with respect to the nodal coordinates.
 */
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


        if (dim == 3)
        {
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
        if (dim == 3)
        {
//      this->phi_qq[0](2,2) = 0.0 ;
//      this->phi_qq[0](2,dim+2) = 0.0 ;
//      this->phi_qq[0](dim+2,2) = 0.0 ;
//      this->phi_qq[0](dim+2,dim+2) = 0.0 ;
        }
    }
}

}
