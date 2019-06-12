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

#ifndef MKNIXREADERRIGID_H
#define MKNIXREADERRIGID_H

#include <iostream>
#include <fstream>
#include <string>
#if defined(_WIN32) || defined(WIN32) || _WIN64
#  include <io.h>
#else
#  include <unistd.h>
#endif

namespace mknix
{

class Simulation;
class System;

/**
	@author AUTHORS <MAILS>
*/
class ReaderRigid
{
public:
    ReaderRigid();
    ReaderRigid( Simulation*, std::ofstream &, std::ifstream & );
    ~ReaderRigid();
    void readRigidBodies( System* );

private:
    void readRigidBody0DMassPoint( System* );
    void readRigidBody1DBar( System* );
    void readRigidBody1DChain( System* );
    void readRigidBody2DMesh( System* );
    void readRigidBody3DGeneric( System* );

    void readNode( double &, double &, double &, std::string );

private:
    Simulation* theSimulation;
    std::ofstream * output;
    std::ifstream * input; // file to read points from
};

}

#endif
