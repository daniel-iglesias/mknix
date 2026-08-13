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

#ifndef MKNIXLOADTHERMALBOUNDARY1D_H
#define MKNIXLOADTHERMALBOUNDARY1D_H

#include "common.h"
#include "LMX/lmx.h"
#include <map>


namespace mknix
{
class Point;

/**
	@author AUTHORS <MAILS>
*/
class LoadThermalBoundary1D
{
public:
    LoadThermalBoundary1D();

//    LoadThermalBoundary1D( /*double, double, double*/ );

    /*virtual */~LoadThermalBoundary1D();

    void loadFile(const std::string& fileName);

    void loadFile(const std::string& fileName, double /*key1*/, double /*key2*/);

    void loadFile2D(const std::string& fileName);

    void loadTimeFile(const std::string& fileName);

    void scaleLoad(double);

    double getLoadThermalBoundary1D( Point* );

protected:
    std::map<double, double> loadmap;
    std::map<double, std::map<double, double>> loadmap2D;
    std::map<double, double> timemap;
    double scaleFactor;
};

}

#endif
