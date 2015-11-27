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

    virtual void assembleExternalHeat( lmx::Vector<data_type> & );

    virtual void outputToFile( std::ofstream* )
    {}
    
    void getMaxTemp(double&);

protected:
    std::vector<Node*> nodes;
    data_type externalHeat;
};

}

#endif
