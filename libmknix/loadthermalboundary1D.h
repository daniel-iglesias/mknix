/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of Nemesis.                                             *
 *                                                                            *
 *  Nemesis is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  Nemesis is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with Nemesis.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#ifndef MKNIXLOADTHERMALBOUNDARY1D_H
#define MKNIXLOADTHERMALBOUNDARY1D_H

#include "common.h"
#include "LMX/lmx.h"
#include <map>


namespace mknix {
class Point;

/**
	@author AUTHORS <MAILS>
*/
class LoadThermalBoundary1D {
public:
    LoadThermalBoundary1D();

//    LoadThermalBoundary1D( /*double, double, double*/ );

    /*virtual */~LoadThermalBoundary1D();

    void loadFile(std::string);
    
    void loadTimeFile(std::string);
    
    void scaleLoad(double);

    double getLoadThermalBoundary1D( Point* );
    
protected:
    std::map<double, double> loadmap;
    std::map<double, double> timemap;
};

}

#endif
