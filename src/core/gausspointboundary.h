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

#ifndef MKNIXGAUSSPOINTBOUNDARY_H
#define MKNIXGAUSSPOINTBOUNDARY_H

#include <LMX/lmx.h>
#include <core/point.h>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file gausspoint.h

  \brief Point for numerical integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix
{

class LoadThermalBoundary1D;

/**
@author Daniel Iglesias
*/
class GaussPointBoundary : public Point
{

public:
    GaussPointBoundary();

    GaussPointBoundary( int dim_in,
                        double alpha_in,
                        double weight_in,
                        double jacobian_in,
                        int num_in,
                        double coor_x,
                        double dc_in );

    GaussPointBoundary( int dim_in,
                        double alpha_in,
                        double weight_in,
                        double jacobian_in,
                        int num_in,
                        double coor_x,
                        double coor_y,
                        double dc_in );

    ~GaussPointBoundary();

    virtual void shapeFunSolve( std::string, double ) override;

    void computeQext( LoadThermalBoundary1D* );
    void assembleQext( lmx::Vector<data_type> & );

//     void gnuplotOutStress( std::ofstream& );

protected:
    void initializeMatVecs();

protected:
    int num;
    double weight;

    lmx::Vector<data_type> Qext;
};

} //Namespace mknix

#endif
