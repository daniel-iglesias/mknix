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


#include "loadthermalbody.h"

#include <core/point.h>
#include <simulation/simulation.h>

#ifdef HAVE_VTK
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#endif

namespace mknix
{

/**
 * @brief Default constructor. Initializes the constant heat source value to zero.
 */
LoadThermalBody::LoadThermalBody()
    : constantSouce(0.)
{
}

/**
 * @brief Destructor.
 */
LoadThermalBody::~LoadThermalBody( /*double , double, double*/ )
{
}

/**
 * @brief Loads a time-dependent load scale factor from a file.
 *
 * Reads pairs of (time, load_factor) values from the specified file and
 * stores them in the internal time map for use during simulation.
 *
 * @param fileName Path to the file containing time and load factor pairs.
 */
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
    m_hasTimeScale = !m_time.empty();
}

/**
 * @brief Computes the thermal body load at a given point. *
 * Returns the heat source value at the location of the provided point.
 * If a 2D spatial distribution is available it is interpolated using the
 * point's X and Y coordinates; if only a 1D distribution is available it
 * is interpolated using the X coordinate; otherwise the constant source
 * value is used. The result is further scaled by the time-dependent factor
 * interpolated from the time map, if one has been loaded.
 *
 * @param point_in Pointer to the Point at which the load is evaluated.
 * @return The thermal body load value at the given point and current simulation time.
 */
double LoadThermalBody::getLoadThermalBody( Point* point_in )
{
    double load;
    switch (m_sourceType)
    {
#ifdef HAVE_VTK
        case SourceType::SOURCE_VTK:
            load = interpolateVTK(point_in->getX(), point_in->getY(), point_in->getZ(),
                                  m_vtkLocator, m_vtkCell, m_vtkScalars);
            break;
#endif
        case SourceType::SOURCE_3D:
            load = interpolate3D(point_in->getX(), point_in->getY(), point_in->getZ(), m_source3D);
            break;
        case SourceType::SOURCE_2D:
            load = interpolate2D(point_in->getX(), point_in->getY(), m_source2D);
            break;
        case SourceType::SOURCE_1D:
            load = interpolate1D(point_in->getX(), m_source);
            break;
        default:
            load = constantSouce;
            break;
    }

    if (!m_hasTimeScale)
    {
        return load;
    }

    return load * interpolate1D(Simulation::getTime(), m_time);
}

}
