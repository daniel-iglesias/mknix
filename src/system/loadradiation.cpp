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

#include "loadradiation.h"

namespace mknix
{

/**
 * @brief Default constructor.
 */
Radiation::Radiation()
    : Load()
{
}

/**
 * @brief Destructor.
 */
Radiation::~Radiation()
{
}


/**
 * @brief Stores a radiation voxel value at the given spatial coordinates.
 * @param x_in     X coordinate of the voxel.
 * @param y_in     Y coordinate of the voxel.
 * @param z_in     Z coordinate of the voxel.
 * @param value_in Radiation intensity value at this voxel.
 */
void Radiation::addVoxel(double x_in, double y_in, double z_in, double value_in)
{
    radMap[z_in][y_in][x_in] = value_in;
//  cout << "voxel read: (" << x_in << ", " << y_in << ", " << z_in << ") = " << value_in << endl;
}

/**
 * @brief Writes the radiation map (regular voxel grid) to the output file.
 * @param outFile Pointer to the output file stream.
 */
void Radiation::outputToFile(std::ofstream * outFile)
{
    *outFile << "RADIATION" << "\t";
    // Computing size of regular grid
    *outFile
            << radMap.begin()->second.begin()->second.size() << "\t"
            << radMap.begin()->second.size() << "\t"
            << radMap.size() << endl;

    for (auto& z : radMap)
    {
        for (auto& y : z.second)
        {
            for (auto& x : y.second)
            {
                *outFile << x.first << "\t"
                         << y.first << "\t"
                         << z.first << "\t"
                         << x.second << endl;
            }
        }
    }
}

}
