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

#ifndef MKNIXLOADTHERMALBODYFLUENCE_H
#define MKNIXLOADTHERMALBODYFLUENCE_H

#include "common.h"
#include "LMX/lmx.h"
#include <map>

#include "loadthermalbody.h"

namespace mknix
{
    class LoadThermalBody;
    class Point;

/// Identifies which spatial source is active, in ascending priority order.
enum class SourceType { CONSTANT = 0, SOURCE_1D = 1, SOURCE_2D = 2, SOURCE_3D = 3, SOURCE_VTK = 4 };

/**
	@author AUTHORS <MAILS>
*/
class LoadThermalBody
{
public:
    LoadThermalBody();

//    LoadThermalBody( /*double, double, double*/ );

    virtual ~LoadThermalBody();

    virtual double getLoadThermalBody( Point* );

    virtual void setConstantValue( double value_in  )
    {
        constantSouce = value_in;
    }

    void loadTimeFile(const std::string& fileName);

    void addLoad( double distance_in, double load_in )
    {
        m_source[distance_in] = load_in;
        if (m_sourceType < SourceType::SOURCE_1D)
            m_sourceType = SourceType::SOURCE_1D;
    }

    void addLoad( double key1_in, double key2_in, double load_in )
    {
        m_source2D[key1_in][key2_in] = load_in;
        if (m_sourceType < SourceType::SOURCE_2D)
            m_sourceType = SourceType::SOURCE_2D;
    }

    void addLoad( double key1_in, double key2_in, double key3_in, double load_in )
    {
        m_source3D[key1_in][key2_in][key3_in] = load_in;
        m_sourceType = SourceType::SOURCE_3D;
    }

    void loadFile3D(const std::string& fileName)
    {
        readFile3D(fileName, m_source3D);
        m_sourceType = SourceType::SOURCE_3D;
    }

#ifdef HAVE_VTK
    void loadFileVTK(const std::string& fileName)
    {
        if (readFileVTK(fileName, m_vtkGrid, m_vtkLocator, m_vtkCell, m_vtkScalars))
        {
            m_sourceType = SourceType::SOURCE_VTK;
        }
    }
#endif

protected:
    double constantSouce;
    SourceType m_sourceType = SourceType::CONSTANT;
    bool m_hasTimeScale = false;
    std::map<double, double> m_source;
    std::map<double, std::map<double, double>> m_source2D;
    std::map<double, std::map<double, std::map<double, double>>> m_source3D;
    std::map<double, double> m_time;

#ifdef HAVE_VTK
    vtkSmartPointer<vtkUnstructuredGrid> m_vtkGrid;
    vtkSmartPointer<vtkCellLocator>      m_vtkLocator;
    vtkSmartPointer<vtkGenericCell>      m_vtkCell;
    vtkDataArray*                        m_vtkScalars = nullptr;
#endif
};

}

#endif
