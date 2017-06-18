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

#ifndef MKNIXTHERMALBODY_H
#define MKNIXTHERMALBODY_H

#include "LMX/lmx.h"
#include "common.h"

namespace mknix {

class Node;
class Cell;
class LoadThermalBody;

/**
	@author AUTHORS <MAILS>
*/
class ThermalBody {

public:

    ThermalBody();

    ThermalBody( std::string );

    virtual ~ThermalBody();

//     virtual std::string getType() = 0;
//
//     virtual void setType( std::string type_in ) = 0;
//
//     virtual void setFormulation( std::string formulation_in )
//     {
//       formulation = formulation_in;
//       cout << formulation << endl;
//     }

    virtual void initialize( );

    virtual void calcCapacityMatrix( );

    virtual void calcConductivityMatrix( );

    virtual void calcExternalHeat( );

    virtual void assembleCapacityMatrix( SparseMatrix<data_type> & );

    virtual void assembleConductivityMatrix( SparseMatrix<data_type> & );

    virtual void assembleExternalHeat( VectorX<data_type> & );

    virtual void setOutput( std::string );

    virtual void outputStep
    ( const VectorX<data_type>&, const VectorX<data_type>& );

    virtual void outputStep( const VectorX<data_type>& );

    virtual void outputToFile( std::ofstream* );

    virtual void addNode( Node* node_in )
    {
        nodes.push_back( node_in );
    }

    virtual Node* getNode( int node_number )
    {
        return nodes[node_number];
    }

    virtual Node* getLastNode( )
    {
        return nodes.back();
    }

    virtual void addCell( int num, Cell* cell_in )
    {
        cells[num] = cell_in;
    }

    // Temporary, should be a pointer to a load class
    virtual void setLoadThermal( LoadThermalBody* theLoad )
    {
        loadThermalBody = theLoad;
    }

protected:
    std::string title;
//     std::string formulation;
    bool computeEnergy;
    std::vector<Node*> nodes;
    std::map<int,Cell*> cells; /**< Map of integration cells. */
    std::vector< VectorX<data_type>* > temperature;
    LoadThermalBody* loadThermalBody;

};

}

#endif
