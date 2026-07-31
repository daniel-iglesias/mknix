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
#include <fstream>

namespace mknix
{

/**
 * @brief Default constructor.
 */
LoadThermal::LoadThermal()
{
}

/**
 * @brief Constructor with node and fluence value.
 * @param node_in    Pointer to the node where the heat load is applied.
 * @param fluence_in Heat fluence (thermal load) value.
 */
LoadThermal::LoadThermal(Node * node_in, double fluence_in)
{
    nodes.push_back(node_in);
    externalHeat = fluence_in;
}

/**
 * @brief Destructor.
 */
LoadThermal::~LoadThermal()
{
}

/**
 * @brief Appends the X coordinates of the loaded nodes to the provided vector.
 * @param x_coordinates Vector to which node X coordinates are appended.
 */
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

/**
 * @brief Loads a time-dependent load scale factor from a file.
 *
 * Reads pairs of (time, load_factor) values from the specified file and
 * stores them in the internal time map for use during simulation.
 *
 * @param fileName Path to the file containing time and load factor pairs.
 */
void LoadThermal::loadTimeFile(const std::string& fileName)
{
    std::ifstream power;
    power.open(fileName);
    if (power.is_open())
    {
        double t, load;
        while (power >> t)
        {
            power >> load;
            m_time[t] = load;
        }
    }
    else
    {
        cerr << "ERROR: TIME FILE NOT FOUND!!!" << endl;
    }
    m_hasTimeScale = !m_time.empty();
}

/**
 * @brief Assembles the thermal load value into the global external heat vector.
 * @param globalExternalHeat Reference to the global external heat vector.
 */
void LoadThermal::assembleExternalHeat
(lmx::Vector<data_type>& globalExternalHeat)
{
    double load = externalHeat;
    if (m_hasTimeScale)
    {
        load *= interpolate1D(Simulation::getTime(), m_time);
    }
    auto nodesSize = nodes.size();
    for (auto i = 0u; i < nodesSize; ++i)
    {
        if (nodes[i]->getNumber() >= 0)
        {
            globalExternalHeat(nodes[i]->getNumber()) += load;
        }
    }
}

/**
 * @brief Updates maxTemp_in with the maximum temperature among all loaded nodes.
 * @param maxTemp_in Reference to the running maximum temperature value.
 */
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
