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
#ifndef MKNIXRIGIDBODY2D_H
#define MKNIXRIGIDBODY2D_H

#include "bodyrigid.h"

namespace mknix {
class Node;

/**
	@author AUTHORS <MAILS>
*/
class RigidBody2D : public RigidBody
{
public:
    RigidBody2D();

    RigidBody2D( std::string, Node*, Node*, Node* );

    ~RigidBody2D();

    std::string getType() {
        return std::string("GENERIC2D");
    }

    void setInertia( double, int );
    void setPosition( std::vector<double>& );

    void calcMassMatrix( );

    void calcExternalForces( );

    void addNode( Node* );

private:
    double Ixx, Iyy, Ixy;
    double rho;

};

}

#endif
