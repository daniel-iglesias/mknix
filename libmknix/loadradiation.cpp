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
#include "loadradiation.h"
#include "simulation.h"

namespace mknix {

Radiation::Radiation()
    : Load()
{
}

Radiation::~Radiation()
{
}


void Radiation::addVoxel( double x_in, double y_in, double z_in, double value_in )
{
    radMap[z_in][y_in][x_in] = value_in;
//  cout << "voxel read: (" << x_in << ", " << y_in << ", " << z_in << ") = " << value_in << endl;
}

void Radiation::outputToFile(std::ofstream * outFile)
{
    std::map<double, std::map<double, std::map<double, double> > >::iterator it_z;
    std::map<double, std::map<double, double> >::iterator it_y;
    std::map<double, double>::iterator it_x;


    *outFile << "RADIATION" << "\t";
    // Computing size of regular grid
    *outFile
            << radMap.begin()->second.begin()->second.size() << "\t"
            << radMap.begin()->second.size() << "\t"
            << radMap.size() << endl;

    for( it_z = radMap.begin();
            it_z!= radMap.end();
            ++it_z)
    {
        for(it_y = it_z->second.begin();
                it_y!= it_z->second.end();
                ++it_y)
        {
            for(it_x = it_y->second.begin();
                    it_x!= it_y->second.end();
                    ++it_x)
            {
                *outFile << it_x->first << "\t"
                         << it_y->first << "\t"
                         << it_z->first << "\t"
                         << it_x->second << endl;
            }
        }
    }
}

}
