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
