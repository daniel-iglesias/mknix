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

#ifdef HAVE_VTK
#include <vtkSmartPointer.h>
class vtkUnstructuredGrid;
class vtkCellLocator;
class vtkDataArray;
class vtkGenericCell;
#endif

namespace mknix
{

typedef double data_type;

double interpolate1D( double, const std::map<double,double>& );
double interpolate2D(double key1,
                     double key2,
                     const std::map<double, std::map<double, double>>& the2DMap);
double interpolate3D(double key1,
                     double key2,
                     double key3,
                     const std::map<double, std::map<double, std::map<double, double>>>& the3DMap);

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

bool isNumber(const std::string& str);

void readFile2D(const std::string& fileName,
                std::map<double, std::map<double, double>>& dest);

void readFile3D(const std::string& fileName,
                std::map<double, std::map<double, std::map<double, double>>>& dest);

#ifdef HAVE_VTK
bool readFileVTK(const std::string& fileName,
                 vtkSmartPointer<vtkUnstructuredGrid>& grid,
                 vtkSmartPointer<vtkCellLocator>& locator,
                 vtkSmartPointer<vtkGenericCell>& cell,
                 vtkDataArray*& scalars);

double interpolateVTK(double x, double y, double z,
                      vtkCellLocator* locator,
                      vtkGenericCell* cell,
                      vtkDataArray* scalars);
#endif
}

#endif
