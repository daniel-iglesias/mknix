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
