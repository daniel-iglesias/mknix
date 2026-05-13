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

#ifndef MKNIXMATERIAL_H
#define MKNIXMATERIAL_H

#include "LMX/lmx.h"
// #include <gmm/gmm_solver_idgmres.h>
#include "common.h"

namespace mknix
{
    class Point;

/**
  @author Daniel Iglesias
 */
class Material
{
private:
    int dim; /**< Dimension of space */
    bool b_porous; /**< Is material porous? */
    double capacity; /**< Thermal specific capacity */
    double kappa; /**< Thermal conductivity */
    double young; /**< Young Modulus. */
    double poisson; /**< Poisson Modulus. */
    double beta; /**< Thermal expansion*/
    double lambda, mu; /**< Lame's coefficients. */
    double density; /**< Density. */
    double resistance; /**< Volumetric thermal resistance for porosity formulation. */
    double fluid_temperature; /**< Temperature of the fluid for porosity formulation. */
    double porosityCapacity = 0.0; /**< Optional porosity capacity for porous materials. */
    std::map<double, double> m_capacity;
    std::map<double, double> m_kapppa;
    std::map<double, double> m_beta;
    std::map<double, double> m_density;
    std::map<double, double> m_resistance; /**< Volumetric thermal resistance for porosity formulation, variable with distance. */
    std::map<double, std::map<double, double>> m_resistance2D; /**< Variable porosity resistance in 2D coordinate tables. */
    std::string m_resistanceCoords; /**< Coordinate mode used for variable resistance interpolation: x, y, z, xy, xz, yz. */
    std::vector<double> fluidTempHistory; /**< History of fluid temperatures at each converged iteration. */
    lmx::DenseMatrix<double> D; /**< Constitutive Linear */
    lmx::DenseMatrix<double> C; /**< Constitutive Saint-Venant Kirchoff*/
    cofe::TensorRank2Sym<2,double> E;
    cofe::TensorRank2Sym<3,double> E3;

public:
    Material();

    ~Material();

    double getE()
    {
        return young;
    }
    double getMu()
    {
        return poisson;
    }
    double getDensity(double temp_in=0)
    {
        if (m_density.empty()) return density;
        else return interpolate1D(temp_in, m_density);
    }
    double getCapacity(double temp_in=0)
    {
        if (m_capacity.empty()) return capacity;
        else return interpolate1D(temp_in, m_capacity);
    }
    double getKappa(double temp_in=0)
    {
        if (m_kapppa.empty()) return kappa;
        else return interpolate1D(temp_in, m_kapppa);
    }
    double getBeta(double temp_in=0)
    {
        if (m_beta.empty()) return beta;
        else return interpolate1D(temp_in, m_beta);
    }
    lmx::DenseMatrix<double>& getD()
    {
        return D;    // be careful, returns a writable reference!!!
    }
    lmx::DenseMatrix<double>& getC()
    {
        return C;    // be careful, returns a writable reference!!!
    }
    double getPorosityResistance(Point*);
    double computePorosityLoad(Point*);

    double getPorosityFluidTemp()
    {
        return fluid_temperature;
    }

    void update(int convergence)
    {
        if (b_porous)
        {
            if (convergence == 1)
            {
                fluidTempHistory.push_back(fluid_temperature);
            }
            else if (convergence == 0)
            {
                    fluid_temperature = fluidTempHistory.back();
            }
        }
    }

    const std::vector<double>& getFluidTempHistory() const
    {
        return fluidTempHistory;
    }

    void setPorosityFluidTemp(double temp_in)
    {
        fluid_temperature = temp_in;
    }

    inline bool isPorous()
    {
        return b_porous;
    }

    void setThermalProps( double capacity_in, double kappa_in, double beta_in, double density_in );
    void setPorosityProps(double resistance_in, double fluid_temperature_in, double porosityCapacity_in = 0.0);
    void setMechanicalProps( int dim_in, double young_in, double poisson_in, double density_in );

    void addThermalCapacity( double temp_in, double capacity_in)
    {
        m_capacity[temp_in] = capacity_in;
    }
    void addThermalConductivity( double temp_in, double conductivity_in)
    {
        m_kapppa[temp_in] = conductivity_in;
    }
    void addThermalExpansion( double temp_in, double beta_in)
    {
        m_beta[temp_in] = beta_in;
    }
    void addThermalDensity( double temp_in, double density_in)
    {
        m_density[temp_in] = density_in;
    }
    void addVariableResistance( double distance_in, double resistance_in)
    {
        m_resistance[distance_in] = resistance_in;
    }
    void addVariableResistance( double key1_in, double key2_in, double resistance_in)
    {
        m_resistance2D[key1_in][key2_in] = resistance_in;
    }
    void setVariableResistanceCoords(const std::string& coords_in)
    {
        m_resistanceCoords = coords_in;
    }

    void computeD();

    void computeC();

    double getCsym( int& i, int& j, int& k, int& l );

    void computeS( cofe::TensorRank2Sym<2,double>& S, const cofe::TensorRank2<2,double>& F, double );
    void computeS( cofe::TensorRank2Sym<3,double>& S, const cofe::TensorRank2<3,double>& F );

    double computeEnergy( const cofe::TensorRank2<2,double>& S );
    double computeEnergy( const cofe::TensorRank2<3,double>& S );

    void outputToFile(std::ofstream * outFile);

private:
    double Cijkl( int& i, int& j, int& k, int& l );

};

}

#endif
