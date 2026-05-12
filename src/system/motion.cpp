/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                             *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#include "motion.h"

#include <core/node.h>

namespace mknix
{

Motion::Motion()
{
}

Motion::Motion(Node * node_in)
    : theNode(node_in)
{
    // Initializing motions to zero
    timeUx[0.] = 0.;
    timeUy[0.] = 0.;
    timeUz[0.] = 0.;
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
    theNode->setUx(interpolate1D(theTime, timeUx));
    theNode->setUy(interpolate1D(theTime, timeUy));
    theNode->setUz(interpolate1D(theTime, timeUz));
}


}
