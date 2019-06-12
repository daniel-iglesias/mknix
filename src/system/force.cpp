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
#include "force.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

Force::Force()
    : Load()
{
}

Force::Force(Node * node_in, double fx_in, double fy_in, double fz_in )
    : Load()
{
    nodes.push_back( node_in );
    externalForces.resize( nodes.size()*Simulation::getDim() );
    externalForces(0) = fx_in;
    externalForces(1) = fy_in;
    if(Simulation::getDim() == 3)
        externalForces(2) = fz_in;
}


Force::~Force()
{
}

void Force::outputToFile(std::ofstream * outFile)
{
    // Nothing yet...
}

}
