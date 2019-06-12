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

#ifndef CELLTETRAHEDRON_H
#define CELLTETRAHEDRON_H

#include "LMX/cofe_TensorRank2.h"
#include "cell.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file celltetrahedron.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace mknix
{

/**
@author Daniel Iglesias
*/
class CellTetrahedron : public Cell
{
protected:
    cofe::TensorRank2<3, double> points; /**< position of vertex points */

public:
    CellTetrahedron();

    CellTetrahedron(Material&,
                    std::string,
                    double, int,
                    Point *,
                    Point *,
                    Point *,
                    Point *
                   );

    virtual ~CellTetrahedron();

    void gnuplotOut(std::ofstream&, std::ofstream&);

protected:
    void createGaussPoints();


};

} //Namespace mknix

#endif
