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

#include "loadthermalbody.h"

#include <core/point.h>
#include <simulation/simulation.h>

namespace mknix
{

LoadThermalBody::LoadThermalBody()
    : constantSouce(0.)
{
}

LoadThermalBody::~LoadThermalBody( /*double , double, double*/ )
{
}

void LoadThermalBody::loadTimeFile(const std::string& fileName)
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
}

double LoadThermalBody::getLoadThermalBody( Point* point_in )
{
    double load = constantSouce;

    if (!m_source2D.empty())
    {
        load = interpolate2D(point_in->getX(), point_in->getY(), m_source2D);
    }
    else if (!m_source.empty())
    {
        load = interpolate1D(point_in->getX(), m_source);
    }

    if (m_time.empty())
    {
        return load;
    }

    return load * interpolate1D(Simulation::getTime(), m_time);
}

}
