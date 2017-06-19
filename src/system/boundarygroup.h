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

#ifndef MKNIXBOUNDARYGROUP_H
#define MKNIXBOUNDARYGROUP_H

#include "common.h"
#include "LMX/lmx.h"
#include <gpu/cpu_run_type.h>


namespace mknix {

class Point;

class Node;

class CellBoundary;

class LoadThermalBoundary1D;

/**
	@author AUTHORS <MAILS>
*/
class BoundaryGroup
{

public:
    BoundaryGroup();

//     BoundaryGroup( std::string );

    virtual ~BoundaryGroup();

    virtual void initialize();

    virtual void calcExternalHeat();

    virtual void assembleExternalHeat(lmx::Vector<data_type>&);
    virtual void assembleExternalHeat(VectorX<data_type>&);

//     virtual void assembleExternalForces( lmx::Vector<data_type> & ) = 0;

    virtual void addNode(Node * node_in)
    {
        this->nodes.push_back(node_in);
    }

    virtual void addCell(CellBoundary * cell_in)
    {
        int num = cells.size();
        this->cells[num] = cell_in;
    }

    // Temporary, should be a pointer to a load class
    virtual void setLoadThermal(LoadThermalBoundary1D * theLoad)
    {
        loadThermalBoundaryGroup = theLoad;
    }

    // Temporary, should be a pointer to a load class
    virtual LoadThermalBoundary1D* getLoadThermal()
    {
      std::cout<< "  virtual LoadThermalBoundary1D* getLoadThermal()" << std::endl;
        return loadThermalBoundaryGroup;
    }


protected:
    std::vector<Node *> nodes;
    std::map<int, CellBoundary *> cells;
    /**< Map of integration cells. */
    LoadThermalBoundary1D * loadThermalBoundaryGroup;

};

}

#endif
