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

#include "loadthermalboundary1D.h"
#include "simulation.h"
#include "point.h"

namespace mknix {

LoadThermalBoundary1D::LoadThermalBoundary1D()
{
}


LoadThermalBoundary1D::~LoadThermalBoundary1D( /*double , double, double*/ )
{
}

void LoadThermalBoundary1D::loadFile(std::string fileName)
{
    std::ifstream power;
    power.open(fileName);
    if (power.is_open()) {
        double X, load;
        while(power >> X) {
            power >> load;
            loadmap[X] = load;
        }
    }
    else {
        cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
    }
}

void LoadThermalBoundary1D::loadTimeFile(std::string fileName)
{
    std::ifstream power;
    power.open(fileName);
    if (power.is_open()) {
        double t, load;
        while(power >> t) {
            power >> load;
            timemap[t] = load;
        }
    }
    else {
        cerr << "ERROR: TIME FILE NOT FOUND!!!" << endl;
    }
}

void LoadThermalBoundary1D::scaleLoad(double loadFactor_in)
{
    std::map<double, double>::iterator it_loadmap;
    for(it_loadmap = loadmap.begin();
	it_loadmap!= loadmap.end();
	++it_loadmap){
      it_loadmap->second *= loadFactor_in;
    }
cout << "SCALE: " << loadFactor_in << " applied." << endl;
}



double LoadThermalBoundary1D::getLoadThermalBoundary1D( Point* thePoint )
{
//  cout << Simulation::getTime() << endl;

    if (loadmap.size() == 0) cerr << "ERROR: LOAD FILE NOT FOUND!!!" << endl;
    if (timemap.size() == 0) { 
      return mknix::interpolate1D(thePoint->getX(), loadmap);
    }
    else{
      return mknix::interpolate1D(thePoint->getX(), loadmap) 
	   * mknix::interpolate1D(Simulation::getTime(), timemap) ;
    }
}

}
