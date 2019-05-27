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

#ifndef MKNIXGAUSSPOINT2D_H
#define MKNIXGAUSSPOINT2D_H

#include "LMX/lmx.h"
#include "gausspoint.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file gausspoint.2Dh

  \brief Point for numerical integration in 2D cells.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix
{

class Material;
class LoadThermal;

/**
@author Daniel Iglesias
*/
class GaussPoint2D : public GaussPoint
{

public:
    GaussPoint2D();

    GaussPoint2D( double alpha_in, double weight_in, double jacobian_in,
                  Material* mat_in, int num, double coor_x, double coor_y,
                  double dc_in, bool stressPoint_in );

    GaussPoint2D( double alpha_in, double weight_in, double jacobian_in,
                  Material* mat_in, int num, double coor_x, double coor_y,
                  double coor_z, double dc_in, bool stressPoint_in );

    ~GaussPoint2D();

    void shapeFunSolve( std::string, double ) override;

    void fillFEmatrices( ) override;

    void computeMij( ) override;

    void computeKij( ) override;

    void computeStress( ) override;

    void computeNLStress( ) override;

    void computeFint( ) override;

    void computeFext( ) override;

    void computeNLFint( ) override;

    void computeNLKij( ) override;

    void assembleMij( lmx::Matrix<data_type> & ) override;

    void assembleKij( lmx::Matrix<data_type> & ) override;

    void assembleRi( lmx::Vector<data_type> &, int ) override;

    void assembleFint( lmx::Vector<data_type> & ) override;

    void assembleFext( lmx::Vector<data_type> & ) override;

    double calcPotentialE( const lmx::Vector<data_type> & ) override;

    double calcKineticE( const lmx::Vector<data_type> & ) override;

    double calcElasticE( ) override;

private:
    void initializeMatVecs();

private:
    cofe::TensorRank2<2,double> F2;
    cofe::TensorRank2Sym<2,double> S2;
    cofe::TensorRank2<2,double> P2;
    cofe::TensorRank2Sym<2,double> sigma2;
};

} //Namespace mknix

#endif
