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
#ifndef MKNIXBOUNDARYGROUP_H
#define MKNIXBOUNDARYGROUP_H

#include "common.h"
#include "LMX/lmx.h"


namespace mknix {

class Point;
class Node;
class CellBoundary;
class LoadThermalBoundary1D;

/**
	@author AUTHORS <MAILS>
*/
class BoundaryGroup{

public:
    BoundaryGroup();

//     BoundaryGroup( std::string );

    ~BoundaryGroup();

    virtual void initialize( );

    virtual void calcExternalHeat( );

    virtual void assembleExternalHeat( lmx::Vector<data_type> & );

//     virtual void assembleExternalForces( lmx::Vector<data_type> & ) = 0;

    virtual void addNode( Node* node_in )
    {
        this->nodes.push_back( node_in );
    }

    virtual void addCell( CellBoundary* cell_in )
    {
      int num = cells.size();
      this->cells[num] = cell_in;
    }

    // Temporary, should be a pointer to a load class
    virtual void setLoadThermal( LoadThermalBoundary1D* theLoad )
    {
        loadThermalBoundaryGroup = theLoad;
    }



protected:
    std::vector<Node*> nodes;
    std::map<int,CellBoundary*> cells; /**< Map of integration cells. */
    LoadThermalBoundary1D* loadThermalBoundaryGroup;

};

}

#endif
