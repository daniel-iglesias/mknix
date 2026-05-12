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

#include "loadthermal.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

LoadThermal::LoadThermal()
{
}

LoadThermal::LoadThermal(Node * node_in, double fluence_in)
{
    nodes.push_back(node_in);
    externalHeat = fluence_in;
}

LoadThermal::~LoadThermal()
{
}

void LoadThermal::insertNodesXCoordinates(std::vector<double>& x_coordinates)
{
    auto nodesSize = nodes.size();
    for (auto i = 0u; i < nodesSize; ++i)
    {
        if (nodes[i]->getNumber() >= 0)
        {
            x_coordinates.push_back(nodes[i]->getX());
        }
    }

}

void LoadThermal::assembleExternalHeat
(lmx::Vector<data_type>& globalExternalHeat)
{
    auto nodesSize = nodes.size();
    for (auto i = 0u; i < nodesSize; ++i)
    {
        if (nodes[i]->getNumber() >= 0)
        {
//             for (m=0; m<Simulation::getDim(); ++m) {
            globalExternalHeat(nodes[i]->getNumber())
            += externalHeat; // change of sign!!
//             }
        }
    }
}

void LoadThermal::getMaxTemp(double& maxTemp_in)
{
    auto nodesSize = nodes.size();

    for (auto i = 0u; i < nodesSize; ++i)
    {
// can be improved using fmax()
        if (nodes[i]->getTemp() > maxTemp_in)
        {
            maxTemp_in = nodes[i]->getTemp();
        }
    }
}


}
