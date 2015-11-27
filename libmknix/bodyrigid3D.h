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

#ifndef _RIGIDBODY3D_H_
#define _RIGIDBODY3D_H_

#include "bodyrigid.h"

namespace mknix {
class Node;

class RigidBody3D: public RigidBody
{
public:

    RigidBody3D();

    RigidBody3D( std::string , Node*, Node*, Node*, Node*);

    ~RigidBody3D();

    std::string getType() {
        return std::string("GENERIC3D");
    }

    void setInertia( double, int );
    void setPosition( std::vector<double>& );

    void calcMassMatrix();

    void calcExternalForces();

    void addNode( Node* );

private:
    double Ixx, Iyy, Izz;
    double Ixy, Iyz, Ixz;
    double Pxx, Pyy, Pzz;
};

}
#endif // _RIGIDBODY3D_H_
