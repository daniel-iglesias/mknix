/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * MkniX-develop
 * Copyright (C) Roberto Ortega 2008 <roberto.ortega@inbox.com>
 *
 * MkniX-develop is free software.
 *
 * You may redistribute it and/or modify it under the terms of the
 * GNU General Public License, as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * MkniX-develop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MkniX-develop.  If not, write to:
 * 	The Free Software Foundation, Inc.,
 * 	51 Franklin Street, Fifth Floor
 * 	Boston, MA  02110-1301, USA.
 */

#include "bodyrigid3D.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

RigidBody3D::RigidBody3D() : RigidBody()
{
}

RigidBody3D::RigidBody3D( std::string title_in,
                          Node* node0_in,
                          Node* node1_in,
                          Node* node2_in,
                          Node* node3_in
                        )
    : RigidBody( title_in )
    , Ixy(0), Iyz(0), Ixz(0)
{
    frameNodes.push_back( node0_in );
    frameNodes.push_back( node1_in );
    frameNodes.push_back( node2_in );
    frameNodes.push_back( node3_in );
    this->localMassMatrix.resize( 4*Simulation::getDim(), 4*Simulation::getDim() );
    this->externalForces.resize( 4*Simulation::getDim() );
}

RigidBody3D::~RigidBody3D()
{
}

void RigidBody3D::setInertia(double inertia_in, int axis)
{
    if( axis == 0 ) Ixx = inertia_in;
    else if( axis == 1 ) Iyy = inertia_in;
    else if( axis == 2 ) Izz = inertia_in;
    else if( axis == 3 ) Ixy = inertia_in;
    else if( axis == 4 ) Iyz = inertia_in;
    else if( axis == 5 ) Ixz = inertia_in;
    else cerr << endl << "ERROR: Trying to set inertia out of bounds in RigidBody3D" << endl;
}

void RigidBody3D::setPosition(std::vector<double>& position)
{
    // TODO: Read rotations and check vector size. Should have 3 elements: CoG_x, CoG_y, CoG_z
    this->frameNodes[0]->setX( position[0] );
    this->frameNodes[0]->setY( position[1] );
    this->frameNodes[0]->setZ( position[2] );

    this->frameNodes[1]->setX( position[0] + 1. );
    this->frameNodes[1]->setY( position[1] );
    this->frameNodes[1]->setZ( position[2] );

    this->frameNodes[2]->setX( position[0] );
    this->frameNodes[2]->setY( position[1] + 1. );
    this->frameNodes[2]->setZ( position[2] );

    this->frameNodes[3]->setX( position[0] );
    this->frameNodes[3]->setY( position[1] );
    this->frameNodes[3]->setZ( position[2] + 1. );
}

void RigidBody3D::calcMassMatrix()
{
    // TODO: Define matrix in generic inertia axes
    if( densityFactor != 1)
    {
        mass *= densityFactor;
        Ixx *= densityFactor;
        Iyy *= densityFactor;
        Izz *= densityFactor;
        Ixy *= densityFactor;
        Iyz *= densityFactor;
        Ixz *= densityFactor;
    }
    if( Ixy == 0 && Iyz == 0 && Ixz == 0 )   // Principal inertia axes defined
    {
        Pxx = .5*( -Ixx +Iyy +Izz );
        Pyy = .5*( +Ixx -Iyy +Izz );
        Pzz = .5*( +Ixx +Iyy -Izz );

        /* For Problem in 3D
        *  MassMatrix = | m00*I3 m01*I3 m02*I3 m03*I3 |
        *               | m10*I3 Pxx*I3   0*I3   0*I3 |
        *               | m20*I3   0*I3 Pyy*I3   0*I3 |
        *               | m30*I3   0*I3   0*I3 Pzz*I3 |
        */

        /* m00*I3 */
        this->localMassMatrix( 0, 0) = mass +Pxx +Pyy +Pzz ;
        this->localMassMatrix( 1, 1) = mass +Pxx +Pyy +Pzz ;
        this->localMassMatrix( 2, 2) = mass +Pxx +Pyy +Pzz ;

        /* m01*I3 */
        this->localMassMatrix( 0, 3) = -Pxx ;
        this->localMassMatrix( 1, 4) = -Pxx ;
        this->localMassMatrix( 2, 5) = -Pxx ;

        /* m02*I3 */
        this->localMassMatrix( 0, 6) = -Pyy ;
        this->localMassMatrix( 1, 7) = -Pyy ;
        this->localMassMatrix( 2, 8) = -Pyy ;

        /* m03*I3 */
        this->localMassMatrix( 0, 9) = -Pzz ;
        this->localMassMatrix( 1,10) = -Pzz ;
        this->localMassMatrix( 2,11) = -Pzz ;

        /* m10*I3 */
        this->localMassMatrix( 3, 0) = -Pxx ;
        this->localMassMatrix( 4, 1) = -Pxx ;
        this->localMassMatrix( 5, 2) = -Pxx ;

        /* m20*I3 */
        this->localMassMatrix( 6, 0) = -Pyy ;
        this->localMassMatrix( 7, 1) = -Pyy ;
        this->localMassMatrix( 8, 2) = -Pyy ;

        /* m30*I3 */
        this->localMassMatrix( 0, 0) = -Pzz ;
        this->localMassMatrix(10, 1) = -Pzz ;
        this->localMassMatrix(11, 2) = -Pzz ;

        /* Pxx*I3 */
        this->localMassMatrix( 3, 3) =  Pxx ;
        this->localMassMatrix( 4, 4) =  Pxx ;
        this->localMassMatrix( 5, 5) =  Pxx ;

        /* Pyy*I3 */
        this->localMassMatrix( 6, 6) =  Pyy ;
        this->localMassMatrix( 7, 7) =  Pyy ;
        this->localMassMatrix( 8, 8) =  Pyy ;

        /* Pzz*I3 */
        this->localMassMatrix( 9, 9) =  Pzz ;
        this->localMassMatrix(10,10) =  Pzz ;
        this->localMassMatrix(11,11) =  Pzz ;
    }
//    cout << localMassMatrix << endl;
}

void RigidBody3D::calcExternalForces()
{
    /*For Problem in 3D*/
    this->externalForces(0) = -mass*Simulation::getGravity(0);
    this->externalForces(1) = -mass*Simulation::getGravity(1);
    this->externalForces(2) = -mass*Simulation::getGravity(2);

//    cout << externalForces << endl;
}

void RigidBody3D::addNode(Node* node_in)
{
    mknix::Body::addNode(node_in); // adds node_in to node vector
    this->nodes.back()->addSupportNode(this->frameNodes[0]);
    this->nodes.back()->addSupportNode(this->frameNodes[1]);
    this->nodes.back()->addSupportNode(this->frameNodes[2]);
    this->nodes.back()->addSupportNode(this->frameNodes[3]);
    this->nodes.back()->setJacobian(1./6.);
    this->nodes.back()->setShapeFunType("3D");
    this->nodes.back()->shapeFunSolve("3D", 1.);
}


}
