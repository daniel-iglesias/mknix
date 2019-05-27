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

#ifndef CELLBOUNDARY_H
#define CELLBOUNDARY_H

#include <vector>
#include <string>
#include "common.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file cell.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx
{
template <typename T> class Vector;
template <typename T> class Matrix;
}

namespace mknix
{

class LoadThermalBoundary1D;
class GaussPointBoundary;
class Node;
class Point;

/**
@author Daniel Iglesias
*/
class CellBoundary
{

public:
    CellBoundary();

    CellBoundary( std::string, double, int );

    virtual ~CellBoundary();

    virtual void initialize( std::vector<Node*> & );

    virtual void computeShapeFunctions(  );

    void computeQextGaussPoints( LoadThermalBoundary1D* );

    void assembleQextGaussPoints( lmx::Vector<data_type> & );

//     void outputConnectivityToFile(std::ofstream*);

//     virtual void gnuplotOut( std::ofstream&, std::ofstream& ) = 0;

//     void gnuplotOutStress( std::ofstream& );

protected:
    std::string formulation;
    double alpha;
    int nGPoints; /**< number of Gauss Points */
    std::vector<GaussPointBoundary*> gPoints;
    double jacobian;
    std::vector< Point* > bodyPoints;
    double dc;

};

} //Namespace mknix

#endif
