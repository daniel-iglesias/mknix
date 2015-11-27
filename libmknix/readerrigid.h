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
#ifndef MKNIXREADERRIGID_H
#define MKNIXREADERRIGID_H

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

namespace mknix {

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

    void readNode( double & , double & , double & , std::string );

private:
    Simulation* theSimulation;
    std::ofstream * output;
    std::ifstream * input; // file to read points from
};

}

#endif
