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

#ifndef MKNIXMATERIAL_H
#define MKNIXMATERIAL_H

#include "LMX/lmx.h"
// #include <gmm/gmm_solver_idgmres.h>
#include "common.h"

namespace mknix {

/**
  @author Daniel Iglesias
 */
class Material {
private:
    int dim; /**< Dimension of space */
    double capacity; /**< Thermal specific capacity */
    double kappa; /**< Thermal conductivity */
    double young; /**< Young Modulus. */
    double poisson; /**< Poisson Modulus. */
    double beta; /**< Thermal expansion*/
    double lambda, mu; /**< Lame's coefficients. */
    double density; /**< Density. */
    std::map<double, double> m_capacity;
    std::map<double, double> m_kapppa;
    lmx::DenseMatrix<double> D; /**< Constitutive Linear */
    lmx::DenseMatrix<double> C; /**< Constitutive Saint-Venant Kirchoff*/
    cofe::TensorRank2Sym<2,double> E;
    cofe::TensorRank2Sym<3,double> E3;

public:
    Material();

    ~Material();

    double getE() {
        return young;
    }
    double getMu() {
        return poisson;
    }
    double getDensity() {
        return density;
    }
    double getCapacity(double temp_in=0) {
      if (m_capacity.empty()) return capacity;
      else return interpolate1D(temp_in, m_capacity);
    }
    double getKappa(double temp_in=0) {
      if (m_kapppa.empty()) return kappa;
      else return interpolate1D(temp_in, m_kapppa);
    }
    lmx::DenseMatrix<double>& getD() {
        return D;    // be careful, returns a writable reference!!!
    }
    lmx::DenseMatrix<double>& getC() {
        return C;    // be careful, returns a writable reference!!!
    }

    void setThermalProps( double capacity_in, double kappa_in, double beta_in, double density_in );
    void setMechanicalProps( int dim_in, double young_in, double poisson_in, double density_in );

    void addThermalCapacity( double temp_in, double capacity_in)
    { m_capacity[temp_in] = capacity_in; }
    void addThermalConductivity( double temp_in, double conductivity_in)
    { m_kapppa[temp_in] = conductivity_in; }
    
    void computeD();

    void computeC();

    double getCsym( int& i, int& j, int& k, int& l );

    void computeS( cofe::TensorRank2Sym<2,double>& S, const cofe::TensorRank2<2,double>& F, double );
    void computeS( cofe::TensorRank2Sym<3,double>& S, const cofe::TensorRank2<3,double>& F );

    double computeEnergy( const cofe::TensorRank2<2,double>& S );
    double computeEnergy( const cofe::TensorRank2<3,double>& S );

private:
    double Cijkl( int& i, int& j, int& k, int& l );

};

}

#endif
