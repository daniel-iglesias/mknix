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

#include "loadthermalbody.h"

#include <core/point.h>
#include <simulation/simulation.h>

namespace mknix {

LoadThermalBody::LoadThermalBody()
{
    std::ifstream power;
    power.open("POWER.txt");
    if (power.is_open()) {

        double keyword, keyword_2;

        while(power >> keyword) {
            power >> keyword_2;
            srim[keyword] = keyword_2;
        }
    }
    else {
        cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
    }

}


LoadThermalBody::~LoadThermalBody( /*double , double, double*/ )
{
}

double LoadThermalBody::getLoadThermalBody( Point* thePoint )
{
  std::cout << "--------------------------THIS GETLOADTHERMALBODY!!!------------------------------" << std::endl;
//  cout << Simulation::getTime() << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////
//     For thermal slits:
//   if (srim.size() == 0) cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
//   if (Simulation::getTime() <= 1.E-4){
// //    if ( thePoint->getX() < 10.E-3 )
//     typedef std::map<double, double>::const_iterator i_t;
//
//     i_t i=srim.upper_bound(thePoint->getX()/4.);
//     if(i==srim.end())
//     {
//       return (--i)->second;
//     }
//     if (i==srim.begin())
//     {
//       return i->second;
//     }
//     i_t l=i; --l;
//
//     const double delta=(thePoint->getX()/4.- l->first)/(i->first - l->first);
//     return (delta*i->second +(1-delta)*l->second)/4.;//    else return 0.;
//   }
//   else return 0.;
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
//     For thermal pendulum:
    /////////////////////////////////////////////////////////////////////////////////////////////
    if (srim.size() == 0) cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
    if (Simulation::getTime() <= 0.1) { // permanent
        if ( thePoint->getX() < 5. ) {
            return mknix::interpolate1D(thePoint->getX(), srim); // else return 0.;
        }
    }
    else return 0.;

  return 0;
}

}
