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

#ifndef MKNIXREADER_H
#define MKNIXREADER_H

#include <iostream>
#include <fstream>
#include <string>
#include "sectionreader.h"
#if defined(_WIN32) || defined(WIN32)
#  include <io.h>
#  define chdir _chdir
#  define getcwd _getcwd
#else
#  include <unistd.h>
#endif

namespace mknix {

class ReaderRigid;

class ReaderFlex;

class ReaderConstraints;

class Simulation;

class System;

/**
	@author AUTHORS <MAILS>
*/
class Reader
{
public:
    Reader();

    Reader(Simulation *);

    ~Reader();

    void inputFromFile(const std::string& fileIn);

private:
    void readSystem(System *);

    void readBodyPoints(System *);

    void readLoads(System *);

    void readEnvironment(System *);

    void readMotion(System *);

    void readAnalysis();

    void readSignals(System *);

//    void readNodeName( std::string & , std::string & );
//    void createConstraintNodes( System* ,
//                                std::string & ,
//                                std::string & ,
//                                std::string & ,
//                                std::string &
//                              );
//    void outputConstraintNode( System* ,
//                              std::string & ,
//                              char * ,
//                              std::string & ,
//                              std::string &
//                             );

private:
    Simulation * theSimulation;
    std::ofstream output;
    std::ifstream input; // file to read points from
    ReaderRigid * theReaderRigid; //used to read RIGIDBODIES
    ReaderFlex * theReaderFlex; //used to read FLEXBODIES
    ReaderConstraints * theReaderConstraints; //used to read CONSTRAINTS
    SectionReader sectionReader;
};

}

#endif
