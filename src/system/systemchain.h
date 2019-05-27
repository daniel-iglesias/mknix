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

#ifndef MKNIXSYSTEMCHAIN_H
#define MKNIXSYSTEMCHAIN_H

#include "system.h"
#include <string>

namespace mknix
{

class Simulation;

class ConstraintDistance;

/**
  @author AUTHORS <MAILS>
*/
class SystemChain : public System
{

//   friend class Reader;
//   friend class ReaderFlex;
//   friend class ReaderRigid;
//   friend class ReaderConstraints;
//   friend class Contact;

public:
    SystemChain();

    SystemChain(const char*);

    ~SystemChain();

    void setInterfaceNodeA(double x, double y, double z)
    {
        x0 = x;
        y0 = y;
        z0 = z;
    }

    void setInterfaceNodeB(double x, double y, double z)
    {
        x1 = x;
        y1 = y;
        z1 = z;
    }

    void setProperties(int segments_in, double length_in)
    {
        segments = segments_in;
        length = length_in;
    }

    void setMass(double mass_in)
    {
        mass = mass_in;
    }

    // Note: timeLength[0] is updated in populate() function
    void setTimeLengths(std::map<double, double>& timeLenghts_in)
    {
        timeLenghts = timeLenghts_in;
    }

    Node* getNode(size_t) override;

    void addTimeLenght(double time_in, double length_in)
    {
        timeLenghts[time_in] = length_in;
    }

    void populate(Simulation*, std::string&);

    void update(double) override;

private:
    double length, currentLength, mass;
    int segments;
    double x0, y0, z0;
    double x1, y1, z1;
    std::map<double, double> timeLenghts;
    std::map<std::string, ConstraintDistance*> localConstraintDistance;

};

}

#endif
