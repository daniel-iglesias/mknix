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

#include "bodyrigid0D.h"

#include <simulation/simulation.h>

namespace mknix {

RigidBodyMassPoint::RigidBodyMassPoint() : RigidBody()
{
}

RigidBodyMassPoint::RigidBodyMassPoint( std::string title_in,
                                        Node * nodeA_in,
                                        double mass_in
                                      )
    : RigidBody( title_in )
{
    this->mass = mass_in;
    frameNodes.push_back( nodeA_in );
    this->localMassMatrix.resize( 1*dim, 1*dim );
    this->externalForces.resize( 1*dim );
}

RigidBodyMassPoint::~RigidBodyMassPoint()
{
}

void RigidBodyMassPoint::setInertia(double inertia_in, int axis)
{
    cerr << endl << "ERROR: Trying to set inertia in RigidBodyMassPoint" << endl;
}

void RigidBodyMassPoint::setPosition(std::vector<double>& position)
{
    // TODO: check vector size. Should have 2 elements: CoG_x, CoG_y
    this->frameNodes[0]->setX( position[0] );
    this->frameNodes[0]->setY( position[1] );
}

void RigidBodyMassPoint::calcMassMatrix()
{
    /* For Problem in 2D
     * MassMatrix = | m 0 |
     *              | 0 m |
     */
    this->localMassMatrix(0,0) = mass ;
    this->localMassMatrix(1,1) = mass ;

    /* For Problem in 3D
     * MassMatrix = | m 0 0 |
     *              | 0 m 0 |
     *              | 0 0 m |
     */
    if ( dim == 3 )
    {
        this->localMassMatrix(2,2) = mass ;
    }
}

void RigidBodyMassPoint::calcExternalForces()
{
    /*For Problem in 2D*/
    this->externalForces(0) = -mass * Simulation::getGravity(0);
    this->externalForces(1) = -mass * Simulation::getGravity(1);

    if ( dim == 3 )
    {
        this->externalForces(2) = -mass * Simulation::getGravity(2);
    }
}

// Node * RigidBodyMassPoint::getDomainNode( std::string name_in )
// {
//   if( name_in == "NODEA" )
//   {
// 	return this->frameNodes[0];
//   }
//   else
//   {
// 	cerr << "ERROR: NO NODE WITH THAT NAME IN BAR" << endl;
//   }
// }

}
