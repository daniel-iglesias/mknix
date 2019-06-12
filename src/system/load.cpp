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
#include "load.h"

#include <simulation/simulation.h>
#include <core/node.h>

namespace mknix
{

Load::Load()
{
}


Load::~Load()
{
}


void Load::assembleExternalForces
(lmx::Vector< data_type > & globalExternalForces)
{
    int nodesSize = nodes.size();
    int i, m;
    for (i=0; i<nodesSize; ++i)
    {
        if (nodes[i]->getNumber() >= 0 )
        {
            for (m=0; m<Simulation::getDim(); ++m)
            {
                globalExternalForces( Simulation::getDim()*nodes[i]->getNumber() + m)
                += externalForces.readElement(Simulation::getDim()*i + m); // change of sign!!
            }
        }
    }
}


}
