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
#ifndef MKNIXREADERFLEX_H
#define MKNIXREADERFLEX_H

#include <iostream>
#include <fstream>

namespace mknix {

class Simulation;
class System;

/**
	@author AUTHORS <MAILS>
*/
class ReaderFlex
{
public:
    ReaderFlex();
    ReaderFlex( Simulation*, std::ofstream &, std::ifstream & );
    ~ReaderFlex();
    void readFlexBodies( System* );

private:
    void readFlexbody( System*, std::string& );
    void readMultiFlexBody( System* );

private:
    Simulation* theSimulation;
    std::ofstream * output;
    std::ifstream * input; // file to read points from
};

}

#endif
