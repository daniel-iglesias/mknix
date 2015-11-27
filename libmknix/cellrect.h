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

#ifndef CELLRECT_H
#define CELLRECT_H

#include "LMX/lmx_mat_dense_matrix.h"
#include "cell.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file cellrect.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {
template<typename T>
class DenseMatrix;
}

namespace mknix {

/**
@author Daniel Iglesias
*/
class CellRect : public Cell
{
private:
    double Ax, Ay, Az;
    double minX, minY, minZ, maxX, maxY, maxZ;
    lmx::DenseMatrix<double> points; /**< position of vertex points */

public:
    CellRect();

    CellRect(Material&,
             std::string,
             double, int,
             double, double,
             double, double,
             double, double,
             double, double,
             double, double,
             double, double,
             double, double);

    CellRect(Material&,
             std::string,
             double, int,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double, double, double
    );

    ~CellRect();

//    void initialize( std::map<int,Node*> & );

    void initialize(std::vector<Node *>&);

    void gnuplotOut(std::ofstream&, std::ofstream&);

private:
    void createGaussPoints(double, double, double = 0);


};

} //Namespace mknix

#endif
