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

#include "constraintfixedaxis.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor.
 */
ConstraintFixedAxis::ConstraintFixedAxis()
    : Constraint()
{
}


/**
 * @brief Constructor for a constraint that fixes the relative displacement along one axis.
 * @param a_in        Pointer to node A.
 * @param b_in        Pointer to node B.
 * @param axisName_in Name of the locked axis: "x", "y", or "z".
 * @param alpha_in    Penalty parameter.
 * @param method_in   Enforcement method keyword.
 */
ConstraintFixedAxis::ConstraintFixedAxis( Node* a_in, Node* b_in, std::string& axisName_in, double& alpha_in, std::string& method_in )
    : Constraint(alpha_in, method_in)
    , axisName(axisName_in)
{
    nodes.push_back( a_in );
    nodes.push_back( b_in );
    size_type total_support_nodes = a_in->getSupportSize(0) + b_in->getSupportSize(0);

//   cout << "node A: " << nodes[0]->getNumber() << endl;
//   cout << "node B: " << nodes[1]->getNumber() << endl;

    if( axisName == "x" || axisName == "X" )
        ro = nodes[1]->getConf(0) - nodes[0]->getConf(0) ;
    if( axisName == "y" || axisName == "Y" )
        ro = nodes[1]->getConf(1) - nodes[0]->getConf(1) ;
    if( axisName == "z" || axisName == "Z" )
        ro = nodes[1]->getConf(2) - nodes[0]->getConf(2) ;

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


/**
 * @brief Destructor.
 */
ConstraintFixedAxis::~ConstraintFixedAxis()
{
}

/**
 * @brief Evaluates phi = rt - ro, the current minus the reference displacement along the fixed axis.
 */
void ConstraintFixedAxis::calcPhi()
{
    if( axisName == "x" || axisName == "X" )
        rt = nodes[1]->getConf(0) - nodes[0]->getConf(0) ;
    if( axisName == "y" || axisName == "Y" )
        rt = nodes[1]->getConf(1) - nodes[0]->getConf(1) ;
    if( axisName == "z" || axisName == "Z" )
        rt = nodes[1]->getConf(2) - nodes[0]->getConf(2) ;

    this->phi[0] = rt - ro ;
}

/**
 * @brief Computes the gradient of phi using the shape-function values of the support nodes.
 */
void ConstraintFixedAxis::calcPhiq()
{
    size_type i,j;
    for(i=0; i<nodes[0]->getSupportSize(0); ++i)
    {
        if( axisName == "x" || axisName == "X" )
            this->phi_q[0](dim*i+0) = -nodes[0]->getShapeFunValue(0,i);
        if( axisName == "y" || axisName == "Y" )
            this->phi_q[0](dim*i+1) = -nodes[0]->getShapeFunValue(0,i);
        if( axisName == "z" || axisName == "Z" )
            this->phi_q[0](dim*i+2) = -nodes[0]->getShapeFunValue(0,i);
    }
    // Now i=supportNodes[0].size()
    for(j=0; j<nodes[1]->getSupportSize(0); ++j)
    {
        if( axisName == "x" || axisName == "X" )
            this->phi_q[0](dim*(i+j)) = nodes[1]->getShapeFunValue(0,j);
        if( axisName == "y" || axisName == "Y" )
            this->phi_q[0](dim*(i+j)+1) = nodes[1]->getShapeFunValue(0,j);
        if( axisName == "z" || axisName == "Z" )
            this->phi_q[0](dim*(i+j)+2) = nodes[1]->getShapeFunValue(0,j);
    }
//  if( axisName == "x" || axisName == "X" ){
//    this->phi_q[0](0) = -1.0 ;
//    this->phi_q[0](1) = -0.0;
//    if (dim == 3)
//      this->phi_q[0](2) = -0.0 ;
//    this->phi_q[0](dim) = +1.0;
//    this->phi_q[0](dim+1) = +0.0;
//    if (dim == 3)
//      this->phi_q[0](dim+2) = +0.0;
//  }
//  else if( axisName == "y" || axisName == "Y" ){
//    this->phi_q[0](0) = -0.0 ;
//    this->phi_q[0](1) = -1.0;
//    if (dim == 3)
//      this->phi_q[0](2) = -0.0 ;
//    this->phi_q[0](dim) = +0.0;
//    this->phi_q[0](dim+1) = +1.0;
//    if (dim == 3)
//      this->phi_q[0](dim+2) = +0.0;
//  }
//  else if( axisName == "z" || axisName == "Z" ){
//    this->phi_q[0](0) = -0.0 ;
//    this->phi_q[0](1) = -0.0;
//    if (dim == 3)
//      this->phi_q[0](2) = -1.0 ;
//    this->phi_q[0](dim) = +0.0;
//    this->phi_q[0](dim+1) = +0.0;
//    if (dim == 3)
//      this->phi_q[0](dim+2) = +1.0;
//  }

}

/**
 * @brief Hessian of phi is zero (linear constraint); this is a no-op.
 */
void ConstraintFixedAxis::calcPhiqq()
{
}

}
