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

namespace mknix {

Radiation::Radiation()
        : Load()
{
}

Radiation::~Radiation()
{
}


void Radiation::addVoxel(double x_in, double y_in, double z_in, double value_in)
{
    radMap[z_in][y_in][x_in] = value_in;
//  cout << "voxel read: (" << x_in << ", " << y_in << ", " << z_in << ") = " << value_in << endl;
}

void Radiation::outputToFile(std::ofstream * outFile)
{
    *outFile << "RADIATION" << "\t";
    // Computing size of regular grid
    *outFile
            << radMap.begin()->second.begin()->second.size() << "\t"
            << radMap.begin()->second.size() << "\t"
            << radMap.size() << endl;

    for (auto& z : radMap) {
        for (auto& y : z.second) {
            for (auto& x : y.second) {
                *outFile << x.first << "\t"
                        << y.first << "\t"
                        << z.first << "\t"
                        << x.second << endl;
            }
        }
    }
}

}
