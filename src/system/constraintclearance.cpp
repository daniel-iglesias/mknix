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

#include "constraintclearance.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

ConstraintClearance::ConstraintClearance()
    : Constraint()
{
}

ConstraintClearance::ConstraintClearance( Node* a_in, Node* b_in, double& rh_in, double& alpha_in, std::string& method_in )
    : Constraint(alpha_in, method_in)
    , rh(rh_in)
{
    nodes.push_back( a_in );
    nodes.push_back( b_in );

    this->stiffnessMatrix.resize(2*dim, 2*dim);
    this->internalForces.resize(2*dim);
    this->lambda.resize(1);
    this->lambda[0]=0.0;
    this->phi.resize(1);
    this->phi_q.resize(1);
    this->phi_q[0].resize(2*dim);
    this->phi_qq.resize(1);
    this->phi_qq[0].resize(2*dim, 2*dim);
}

ConstraintClearance::~ConstraintClearance()
{
}

void ConstraintClearance::calcPhi()
{
    rt = std::sqrt( std::pow( nodes[1]->getConf(0) - nodes[0]->getConf(0), 2 )
                    + std::pow( nodes[1]->getConf(1) - nodes[0]->getConf(1), 2 )
                    + std::pow( nodes[1]->getConf(2) - nodes[0]->getConf(2), 2 ) );

    this->phi[0] = rt*rt - rh*rh ;
}

void ConstraintClearance::calcPhiq()
{
    if( rt > rh )
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
    else if( rt <= rh )
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
    }
}

void ConstraintClearance::calcPhiqq()
{
    if( rt > rh )
    {
        this->phi_qq[0](0,0) =  2.0 ;
        this->phi_qq[0](1,1) =  2.0 ;
        this->phi_qq[0](dim,dim) =  2.0 ;
        this->phi_qq[0](dim+1,dim+1) =  2.0 ;
        this->phi_qq[0](0,dim) = -2.0 ;
        this->phi_qq[0](dim,0) = -2.0 ;
        this->phi_qq[0](1,dim+1) = -2.0 ;
        this->phi_qq[0](dim+1,1) = -2.0 ;
        if (dim == 3)
        {
            this->phi_qq[0](2,2) =  2.0 ;
            this->phi_qq[0](dim+2,dim+2) =  2.0 ;
            this->phi_qq[0](2,dim+2) = -2.0 ;
            this->phi_qq[0](dim+2,2) = -2.0 ;
        }
    }
    else if( rt <= rh )
    {
        //this->phi_qq[0].fillIdentity( 0.0 )
        this->phi_qq[0](0,0) = 0.0 ;
        this->phi_qq[0](1,1) = 0.0 ;
        this->phi_qq[0](dim,dim) = 0.0 ;
        this->phi_qq[0](dim+1,dim+1) = 0.0 ;
        this->phi_qq[0](0,dim) = 0.0 ;
        this->phi_qq[0](dim,0) = 0.0 ;
        this->phi_qq[0](1,dim+1) = 0.0 ;
        this->phi_qq[0](dim+1,1) = 0.0 ;
        if (dim == 3)
        {
            this->phi_qq[0](2,2) = 0.0 ;
            this->phi_qq[0](2,dim+2) = 0.0 ;
            this->phi_qq[0](dim+2,2) = 0.0 ;
            this->phi_qq[0](dim+2,dim+2) = 0.0 ;
        }
    }
}

}
