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

#ifndef MKNIXLOAD_H
#define MKNIXLOAD_H

#include "common.h"
#include "LMX/lmx.h"

namespace mknix {
class Node;

/**
	@author AUTHORS <MAILS>
*/
class Load {
public:
    Load();

    virtual ~Load();

    virtual void assembleExternalForces( lmx::Vector<data_type> & );
    //virtual void assembleExternalForces( VectorX<data_type> & );

    virtual void outputToFile( std::ofstream* )=0;

protected:
    std::vector<Node*> nodes;
    lmx::Vector<data_type> externalForces;
    lmx::Vector<data_type> externalHeat;
};

}

#endif
