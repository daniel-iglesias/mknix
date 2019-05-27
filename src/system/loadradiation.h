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

#ifndef MKNIXRADIATION_H
#define MKNIXRADIATION_H

#include <map>

#include "load.h"

namespace mknix
{

/**
	@author AUTHORS <MAILS>
*/
class Radiation : public Load
{
public:
    Radiation();

    ~Radiation();

    void addVoxel( double, double, double, double );

    void outputToFile( std::ofstream* );

private:
    std::map< double, std::map< double, std::map< double, double > > > radMap; // x,y,z,value
};

}

#endif
