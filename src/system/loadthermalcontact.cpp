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

#include "loadthermalcontact.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Constructs a thermal contact load between two nodes.
 * @param nodeA           Pointer to the first node.
 * @param nodeB           Pointer to the second node.
 * @param filmCoefficient Film (conductance) coefficient used to scale the heat exchange.
 */
LoadThermalContact::LoadThermalContact(Node * nodeA, Node * nodeB, double filmCoefficient)
    : filmCoefficient(filmCoefficient)
{
    nodes.push_back(nodeA);
    nodes.push_back(nodeB);
}

/**
 * @brief Destructor.
 */
LoadThermalContact::~LoadThermalContact()
{
}

/**
 * @brief Assembles the heat exchanged between the two contact nodes into the global vector.
 * @param globalExternalHeat Reference to the global external heat vector.
 */
void LoadThermalContact::assembleExternalHeat
(lmx::Vector<data_type>& globalExternalHeat)
{
    double scale = 1.;
    if (m_hasTimeScale)
    {
        scale = interpolate1D(Simulation::getTime(), m_time);
    }

    const double Ta = nodes[0]->getTemp();
    const double Tb = nodes[1]->getTemp();
    const double heatA = -filmCoefficient * scale * (Ta - Tb);

    if (nodes[0]->getNumber() >= 0)
    {
        globalExternalHeat(nodes[0]->getNumber()) += heatA;
    }
    if (nodes[1]->getNumber() >= 0)
    {
        globalExternalHeat(nodes[1]->getNumber()) -= heatA;
    }
}

}
