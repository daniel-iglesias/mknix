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
#include "bodyrigid2D.h"

#include <simulation/simulation.h>

namespace mknix
{

RigidBody2D::RigidBody2D()
    : RigidBody()
    , Ixx(0)
    , Iyy(0)
    , Ixy(0)
{
}

RigidBody2D::RigidBody2D( std::string title_in,
                          Node * nodeA_in,
                          Node * nodeB_in,
                          Node * nodeC_in
                        )
    : RigidBody( title_in )
    , Ixx(0)
    , Iyy(0)
    , Ixy(0)
{
    // Read frameNodes from theReader:
    // It is assumed that correspond to extremes of bar and
    // it is necessary to change them to CoG and director vector
    frameNodes.push_back( nodeA_in );
    frameNodes.push_back( nodeB_in );
    frameNodes.push_back( nodeC_in );
    // Create the boundary nodes and assign them the support frameNodes and constant function
    // Time to define the relationship between boundary nodes and formulation frameNodes

    this->localMassMatrix.resize( 3*dim, 3*dim );
    this->externalForces.resize( 3*dim );

}

RigidBody2D::~RigidBody2D()
{
}

void RigidBody2D::setInertia(double inertia_in, int axis)
{
    if( axis == 0 ) Ixx = inertia_in;
    else if( axis == 1 ) Iyy = inertia_in;
    else if( axis == 2 ) Ixy = inertia_in;
    else cerr << endl << "ERROR: Trying to set inertia out of bounds in RigidBody2D" << endl;
}
void RigidBody2D::setPosition(std::vector<double>& position)
{
    // TODO: check vector size. Should have 3 elements: CoG_x, CoG_y, rotation angle
    this->frameNodes[0]->setX( position[0] );
    this->frameNodes[0]->setY( position[1] );
    this->frameNodes[1]->setX( position[0] + std::cos(position[2]) );
    this->frameNodes[1]->setY( position[1] + std::sin(position[2]) );
    this->frameNodes[2]->setX( position[0] - std::sin(position[2]) );
    this->frameNodes[2]->setY( position[1] + std::cos(position[2]) );
}


void RigidBody2D::calcMassMatrix()
{
    // m00
    this->localMassMatrix(0,0) = mass + Ixx + Iyy;
    this->localMassMatrix(1,1) = mass + Ixx + Iyy;
    // m01
    this->localMassMatrix(0,Simulation::getDim()) = -Iyy;
    this->localMassMatrix(1,Simulation::getDim()+1) = -Iyy;
    // P00
    this->localMassMatrix(Simulation::getDim(),Simulation::getDim()) =  Iyy;
    this->localMassMatrix(Simulation::getDim()+1,Simulation::getDim()+1) =  Iyy ;
    // m10
    this->localMassMatrix(Simulation::getDim(),0) =  -Iyy;
    this->localMassMatrix(Simulation::getDim()+1,1) =  -Iyy;
    // m02
    this->localMassMatrix(0,2*Simulation::getDim()) = -Ixx;
    this->localMassMatrix(1,2*Simulation::getDim()+1) = -Ixx;
    // P11
    this->localMassMatrix(2*Simulation::getDim(),2*Simulation::getDim()) =  Ixx;
    this->localMassMatrix(2*Simulation::getDim()+1,2*Simulation::getDim()+1) =  Ixx ;
    // m20
    this->localMassMatrix(2*Simulation::getDim(),0) =  -Ixx;
    this->localMassMatrix(2*Simulation::getDim()+1,1) =  -Ixx;


    if (Simulation::getDim() == 3)   // It shouldn't happen, maybe throw an error
    {
        // m00
        this->localMassMatrix(2,2) = mass + Ixx + Iyy;
        // m01
        this->localMassMatrix(2,Simulation::getDim()+2) = -Iyy;
        // P00
        this->localMassMatrix(Simulation::getDim()+2,Simulation::getDim()+2) =  Iyy ;
        // m10
        this->localMassMatrix(Simulation::getDim()+2,2) =  -Iyy;
        // m02
        this->localMassMatrix(2,2*Simulation::getDim()+2) = -Ixx;
        // P11
        this->localMassMatrix(2*Simulation::getDim()+2,2*Simulation::getDim()+2) =  Ixx ;
        // m20
        this->localMassMatrix(2*Simulation::getDim()+2,2) =  -Ixx;
    }
}

void RigidBody2D::calcExternalForces()
{
    this->externalForces(0) = -mass * Simulation::getGravity(0);
    this->externalForces(1) = -mass * Simulation::getGravity(1);
    if (Simulation::getDim() == 3)
        this->externalForces(2) = -mass * Simulation::getGravity(2);
}

void RigidBody2D::addNode(Node* node_in)
{
    mknix::Body::addNode(node_in); // adds node_in to node vector
    this->nodes.back()->addSupportNode(this->frameNodes[0]);
    this->nodes.back()->addSupportNode(this->frameNodes[1]);
    this->nodes.back()->addSupportNode(this->frameNodes[2]);
    this->nodes.back()->setJacobian(0.5);
    this->nodes.back()->setShapeFunType("2D");
    this->nodes.back()->shapeFunSolve("2D", 1.);
}


}
