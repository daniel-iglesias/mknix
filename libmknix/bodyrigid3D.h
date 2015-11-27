/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of Nemesis.                                             *
 *                                                                            *
 *  Nemesis is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  Nemesis is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with Nemesis.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

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
