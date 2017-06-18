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

#ifndef MKNIXLOADTHERMAL_H
#define MKNIXLOADTHERMAL_H

#include "common.h"
#include "LMX/lmx.h"

namespace mknix {
class Node;

/**
	@author AUTHORS <MAILS>
*/
class LoadThermal {
public:
    LoadThermal();

    LoadThermal( Node*, double );

    virtual ~LoadThermal();

    virtual void insertNodesXCoordinates( std::vector<double>& );

    virtual void updateLoad(double load_in)
    { externalHeat = load_in; }

    virtual void assembleExternalHeat( VectorX<data_type> & );

    virtual void outputToFile( std::ofstream* )
    {}

    void getMaxTemp(double&);

protected:
    std::vector<Node*> nodes;
    data_type externalHeat;
};

}

#endif
