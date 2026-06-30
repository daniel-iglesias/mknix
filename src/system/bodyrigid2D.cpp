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

#include "bodyrigid2D.h"

#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor. Initializes the 2D rigid body with zero inertia terms.
 */
RigidBody2D::RigidBody2D()
    : RigidBody()
    , Ixx(0)
    , Iyy(0)
    , Ixy(0)
{
}

/**
 * @brief Constructor with title and three frame nodes (CoG, x-director, y-director).
 * @param title_in  Name identifier for this body.
 * @param nodeA_in  Pointer to the centre-of-gravity node.
 * @param nodeB_in  Pointer to the x-direction frame node.
 * @param nodeC_in  Pointer to the y-direction frame node.
 */
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

/**
 * @brief Destructor.
 */
RigidBody2D::~RigidBody2D()
{
}

/**
 * @brief Sets one component of the planar inertia tensor.
 * @param inertia_in Value to assign.
 * @param axis       Component index: 0=Ixx, 1=Iyy, 2=Ixy.
 */
void RigidBody2D::setInertia(double inertia_in, int axis)
{
    if( axis == 0 ) Ixx = inertia_in;
    else if( axis == 1 ) Iyy = inertia_in;
    else if( axis == 2 ) Ixy = inertia_in;
    else cerr << endl << "ERROR: Trying to set inertia out of bounds in RigidBody2D" << endl;
}
/**
 * @brief Sets the initial pose of the body from a centre-of-mass position and rotation angle.
 * @param position Vector [CoG_x, CoG_y, rotation_angle] defining the initial pose.
 */
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


/**
 * @brief Builds the local mass matrix for a planar rigid body using the body-inertia parameters.
 */
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

/**
 * @brief Computes the gravitational external force vector applied at the CoG node.
 */
void RigidBody2D::calcExternalForces()
{
    this->externalForces(0) = -mass * Simulation::getGravity(0);
    this->externalForces(1) = -mass * Simulation::getGravity(1);
    if (Simulation::getDim() == 3)
        this->externalForces(2) = -mass * Simulation::getGravity(2);
}

/**
 * @brief Adds a domain node, assigning the three frame nodes as support nodes
 *        and solving the 2D shape function.
 * @param node_in Pointer to the domain node to add.
 */
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
