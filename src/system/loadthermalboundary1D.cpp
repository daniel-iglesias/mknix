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


#include "loadthermalboundary1D.h"

#include <core/point.h>
#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor.
 */
LoadThermalBoundary1D::LoadThermalBoundary1D()
{
}


/**
 * @brief Destructor.
 */
LoadThermalBoundary1D::~LoadThermalBoundary1D( /*double , double, double*/ )
{
}

/**
 * @brief Reads (X, load) pairs from a file and stores them in the spatial load map.
 * @param fileName Path to the spatial load data file.
 */
void LoadThermalBoundary1D::loadFile(std::string fileName)
{
    std::ifstream power;
    power.open(fileName);
    if (power.is_open())
    {
        double X, load;
        while (power >> X)
        {
            power >> load;
            loadmap[X] = load;
        }
    }
    else
    {
        cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
    }
}

/**
 * @brief Reads (time, load_factor) pairs from a file and stores them in the time-scale map.
 * @param fileName Path to the time-scale data file.
 */
void LoadThermalBoundary1D::loadTimeFile(std::string fileName)
{
    std::ifstream power;
    power.open(fileName);
    if (power.is_open())
    {
        double t, load;
        while (power >> t)
        {
            power >> load;
            timemap[t] = load;
        }
    }
    else
    {
        cerr << "ERROR: TIME FILE NOT FOUND!!!" << endl;
    }
}

/**
 * @brief Scales all spatial load values by a constant factor.
 * @param loadFactor_in Multiplicative factor to apply to every load entry.
 */
void LoadThermalBoundary1D::scaleLoad(double loadFactor_in)
{
    for (auto el : loadmap)
    {
        el.second *= loadFactor_in;
    }
    cout << "SCALE: " << loadFactor_in << " applied." << endl;
}


/**
 * @brief Returns the thermal boundary load at the given point's X coordinate,
 *        optionally scaled by the time-dependent factor.
 * @param thePoint Pointer to the boundary integration point.
 * @return Interpolated load value at the point's position and current simulation time.
 */
double LoadThermalBoundary1D::getLoadThermalBoundary1D(Point * thePoint)
{
//  cout << Simulation::getTime() << endl;

    if (loadmap.size() == 0) cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
    if (timemap.size() == 0)
    {
        return mknix::interpolate1D(thePoint->getX(), loadmap);
    }
    else
    {
        return mknix::interpolate1D(thePoint->getX(), loadmap)
               * mknix::interpolate1D(Simulation::getTime(), timemap);
    }
}

}
