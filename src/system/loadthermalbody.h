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
    }

    void addLoad( double key1_in, double key2_in, double load_in )
    {
        m_source2D[key1_in][key2_in] = load_in;
    }


protected:
    double constantSouce;
    std::map<double, double> m_source;
    std::map<double, std::map<double, double>> m_source2D;
    std::map<double, double> m_time;
};

}

#endif
