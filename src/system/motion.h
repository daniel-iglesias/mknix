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

#ifndef MKNIXMOTION_H
#define MKNIXMOTION_H

#include "system.h"

namespace mknix {

class Node;

/**
  @author AUTHORS <MAILS>
*/
class Motion : public System
{

public:
    Motion();

    Motion(Node *);

    ~Motion();

    void setNode(Node * node_in)
    {
        theNode = node_in;
    }

    // Note: timeLength[0] is updated in populate() function
    void setTimeUx(std::map<double, double>& timeUX_in)
    {
        timeUx = timeUX_in;
    }

    void setTimeUy(std::map<double, double>& timeUY_in)
    {
        timeUy = timeUY_in;
    }

    void setTimeUz(std::map<double, double>& timeUZ_in)
    {
        timeUz = timeUZ_in;
    }

    void update(double);

private:
    Node * theNode;
    std::map<double, double> timeUx;
    std::map<double, double> timeUy;
    std::map<double, double> timeUz;

};

}

#endif
