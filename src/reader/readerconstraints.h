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
#if defined(_WIN32) || defined(WIN32)
#  include <io.h>
#else
#  include <unistd.h>
#endif

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

    ReaderConstraints(Simulation*, std::ofstream&, std::ifstream&);

    ~ReaderConstraints();

    void readConstraints(System*);

private:
    void readNodeName(std::string&, std::string&);

    void assignConstraintNodes(System* system_in,
                               const std::string& consName,
                               const std::string& bodyTitleA,
                               const std::string& nodeA,
                               const std::string& bodyTitleB,
                               const std::string& nodeB);

    void outputConstraintNode(System* system_in,
                              const std::string& consTitle,
                              const std::string& nodeName,
                              const std::string& bodyTitle,
                              const std::string& node,
                              std::size_t i);

    void outputConstraintThermalNode(System* system_in,
                                     const std::string& consTitle,
                                     const std::string& nodeName,
                                     const std::string& bodyTitle,
                                     const std::string& node,
                                     std::size_t i);

private:
    Simulation* theSimulation;
    std::ofstream* output;
    std::ifstream* input; // file to read points from
    /* temporary pointers to constraint nodes */
    Node* p_nodeA;
    Node* p_nodeB;
};

}

#endif
