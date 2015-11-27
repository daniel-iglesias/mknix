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
