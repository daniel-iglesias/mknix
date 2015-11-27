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

#ifndef _RIGIDBODYMASSPOINT_H_
#define _RIGIDBODYMASSPOINT_H_

#include "bodyrigid.h"

namespace mknix {

class node;

class RigidBodyMassPoint: public RigidBody
{
public:

    RigidBodyMassPoint();

    RigidBodyMassPoint( std::string, Node*, double );

    ~RigidBodyMassPoint();

    void setInertia( double, int );
    void setPosition( std::vector<double>& );

    void calcMassMatrix();

    void calcExternalForces();

//     Node* getDomainNode( std::string );

    std::string getType() {
        return std::string("MASSPOINT");
    }

private:

};
}

#endif // _RIGIDBODYMASSPOINT_H_
