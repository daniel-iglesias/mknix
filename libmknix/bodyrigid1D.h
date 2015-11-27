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
#ifndef MKNIXRIGIDBAR_H
#define MKNIXRIGIDBAR_H

#include "bodyrigid.h"

namespace mknix {
class Node;

/**
	@author AUTHORS <MAILS>
*/
class RigidBar : public RigidBody
{
public:
    RigidBar();

    RigidBar( std::string, Node*, Node*, double );

    ~RigidBar();

    std::string getType() {
        return std::string("BAR");
    }

    void setInertia( double, int );
    void setPosition( std::vector<double>& );

    void calcMassMatrix( );

    void calcExternalForces( );

    void addNode( Node* );

//     Node* getDomainNode( std::string );

    void writeBoundaryNodes( std::vector<Node*>& );

    void writeBoundaryConnectivity(std::vector< std::vector<Node*> >&);


private:
    double lenght, rho, Jo;

};

}

#endif
