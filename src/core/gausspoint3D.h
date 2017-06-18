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

#ifndef MKNIXGAUSSPOINT3D_H
#define MKNIXGAUSSPOINT3D_H

#include "LMX/lmx.h"
#include "gausspoint.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file gausspoint.3Dh

  \brief Point for numerical integration in 3D cells.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

class Material;
class LoadThermal;

/**
@author Daniel Iglesias
*/
class GaussPoint3D : public GaussPoint {

public:
    GaussPoint3D();

    GaussPoint3D( double alpha_in, double weight_in, double jacobian_in,
                  Material* mat_in, int num, double coor_x, double coor_y,
                  double coor_z, double dc_in, bool stressPoint_in );

    ~GaussPoint3D();

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

    void assembleMij( SparseMatrix<data_type> & ) override;

    void assembleKij( SparseMatrix<data_type> & ) override;

    void assembleRi( VectorX<data_type> &, int ) override;

    void assembleFint( VectorX<data_type> & ) override;

    void assembleFext( VectorX<data_type> & ) override;

    double calcPotentialE( const VectorX<data_type> & ) override;

    double calcKineticE( const VectorX<data_type> & ) override;

    double calcElasticE( ) override;

private:
    void initializeMatVecs();

private:
    cofe::TensorRank2<3,double> F3;
    cofe::TensorRank2Sym<3,double> S3;
    cofe::TensorRank2<3,double> P3;
    cofe::TensorRank2Sym<3,double> sigma3;
};

} //Namespace mknix

#endif
