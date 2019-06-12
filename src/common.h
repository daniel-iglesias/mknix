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

#ifndef MKNIXCOMMON_H
#define MKNIXCOMMON_H

#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>

using namespace std;

namespace mknix
{

typedef double data_type;

double interpolate1D( double, const std::map<double,double>& );

/*!
 * Stand-in for std::make_unique included in C++14
 */
template<class T, class... Args>
std::unique_ptr<T> make_unique(Args&& ... args)
{
    return std::unique_ptr<T>(new T{ std::forward<Args>(args)... });
};

class boxFIR
{
    std::size_t numCoeffs;
    vector<double> b; //Filter coefficients
    vector<double> m; //Filter memories

public:
    boxFIR(int);
    void filter(vector<double> &);
};

std::vector<double> doubles_in_vector( const std::string& );

std::vector< std::vector<double> > read_lines( std::istream& );

}

#endif
