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
#include "motion.h"
#include "system.h"
#include "node.h"

namespace mknix {

Motion::Motion()
{
}

Motion::Motion( Node* node_in )
    : theNode(node_in)
{   // Initializing motions to zero
    timeUx[0.]=0.;
    timeUy[0.]=0.;
    timeUz[0.]=0.;
}

Motion::~Motion()
{
}

void Motion::update(double theTime)
{
//   cout << "Node: "<< theNode->getNumber()
//        << ": U=(" << interpolate1D(theTime, timeUx)
//        << ", " << interpolate1D(theTime, timeUy)
//        << ", " << interpolate1D(theTime, timeUz)
//        << ") " << endl;
    theNode->setUx( interpolate1D(theTime, timeUx) );
    theNode->setUy( interpolate1D(theTime, timeUy) );
    theNode->setUz( interpolate1D(theTime, timeUz) );
}


}
