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

#ifndef MKNIXREADERCONSTRAINTS_H
#define MKNIXREADERCONSTRAINTS_H

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

namespace mknix {

class Node;
class Simulation;
class System;

/**
	@author AUTHORS <MAILS>
*/
class ReaderConstraints
{
public:
    ReaderConstraints();
    ReaderConstraints( Simulation*, std::ofstream &, std::ifstream & );
    ~ReaderConstraints();
    void readConstraints( System* );

private:
    void readNodeName( std::string & , std::string & );
    void assignConstraintNodes( System* ,
                                std::string & ,
                                std::string & ,
                                std::string & ,
                                std::string &
                              );
    void outputConstraintNode( System* ,
                               std::string & ,
                               const char * ,
                               std::string & ,
                               std::string & ,
                               int
                             );
    void outputConstraintThermalNode( System* ,
                               std::string & ,
                               const char * ,
                               std::string & ,
                               std::string & ,
                               int
                             );

private:
    Simulation* theSimulation;
    std::ofstream * output;
    std::ifstream * input; // file to read points from
    /* temporary pointers to constraint nodes */
    Node * p_nodeA;
    Node * p_nodeB;
};

}

#endif
