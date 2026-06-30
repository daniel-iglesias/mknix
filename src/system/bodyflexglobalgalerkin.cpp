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

#include "bodyflexglobalgalerkin.h"

#include <core/cell.h>
#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor.
 */
FlexGlobalGalerkin::FlexGlobalGalerkin()
    : FlexBody()
{
}


/**
 * @brief Constructor with title.
 * @param title_in Name identifier for this body.
 */
FlexGlobalGalerkin::FlexGlobalGalerkin(std::string title_in)
    : FlexBody(title_in)
{
}


/**
 * @brief Destructor.
 */
FlexGlobalGalerkin::~FlexGlobalGalerkin()
{
}

/**
 * @brief Computes the local mass and external force matrices for all cells.
 *        External forces are computed here because they depend on mass and gravity.
 */
void FlexGlobalGalerkin::calcMassMatrix()
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->computeMGaussPoints();
        this->cells[i]->computeFextGaussPoints(); // Depends on M and gravity
    }

}

/**
 * @brief Computes the local internal force vector for all cells (linear or nonlinear formulation).
 */
void FlexGlobalGalerkin::calcInternalForces()
{
    auto end_int = this->cells.size();
    if (formulation == "NONLINEAR")
    {
//         #pragma omp parallel for
        for (auto i = 0u; i < end_int; ++i)
        {
            this->cells[i]->computeNLFintGaussPoints();
        }
    }
    else if (formulation == "LINEAR")
    {
//         #pragma omp parallel for
        for (auto i = 0u; i < end_int; ++i)
        {
            this->cells[i]->computeFintGaussPoints();
        }
    }

}

/**
 * @brief No-op in GlobalGalerkin formulation; external forces are computed in calcMassMatrix().
 */
void FlexGlobalGalerkin::calcExternalForces()
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i)   // Moved to calcMassMatrix(), as M and Gravity are constant
    {
//         this->cells[i]->computeFextGaussPoints( );
    }

}

/**
 * @brief Computes the local tangent stiffness matrix for all cells (linear or nonlinear formulation).
 */
void FlexGlobalGalerkin::calcTangentMatrix()
{
    auto end_int = this->cells.size();
    if (formulation == "NONLINEAR")
    {
        #pragma omp parallel for
        for (auto i = 0u; i < end_int; ++i)
        {
            this->cells[i]->computeNLKGaussPoints();
        }
    }
    else if (formulation == "LINEAR")
    {
        #pragma omp parallel for
        for (auto i = 0u; i < end_int; ++i)
        {
            this->cells[i]->computeKGaussPoints();
        }
    }

}

/**
 * @brief Assembles the local mass matrix into the global mass matrix.
 * @param globalMass Reference to the global mass matrix.
 */
void FlexGlobalGalerkin::assembleMassMatrix
(lmx::Matrix<data_type>& globalMass)
{
    auto end_int = this->cells.size();
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->assembleMGaussPoints(globalMass);
    }

    // Chapuza...
    if (Simulation::getSmoothingType() == "GLOBAL")
    {
        if (smoothingMassMatrix.rows() == 0)
        {
            // Compute smoothing matrix (local mass matrix in 1D)
            smoothingMassMatrix.resize(nodes.size(), nodes.size());
            for (auto i = 0u; i < nodes.size(); ++i)
            {
                for (auto j = 0u; j < nodes.size(); ++j)
                {
                    smoothingMassMatrix.writeElement(globalMass.readElement
                                                     (Simulation::getDim() * (nodes[i]->getNumber()),
                                                      Simulation::getDim() * (nodes[j]->getNumber())),
                                                     i, j);
                }
            }
        }
    }
}

/**
 * @brief Assembles the local internal forces into the global internal force vector.
 * @param globalInternalForces Reference to the global internal force vector.
 */
void FlexGlobalGalerkin::assembleInternalForces
(lmx::Vector<data_type>& globalInternalForces)
{
    auto end_int = this->cells.size();
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->assembleFintGaussPoints(globalInternalForces);
    }

}

/**
 * @brief Assembles the local external forces into the global external force vector.
 * @param globalExternalForces Reference to the global external force vector.
 */
void FlexGlobalGalerkin::assembleExternalForces
(lmx::Vector<data_type>& globalExternalForces)
{
    auto end_int = this->cells.size();
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->assembleFextGaussPoints(globalExternalForces);
    }

}

/**
 * @brief Assembles the local tangent stiffness matrix into the global tangent matrix.
 * @param globalTangent Reference to the global tangent matrix.
 */
void FlexGlobalGalerkin::assembleTangentMatrix(lmx::Matrix<data_type>& globalTangent)
{
    auto end_int = this->cells.size();
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->assembleKGaussPoints(globalTangent);
    }

}

/**
 * @brief Postprocess and store step results for dynamic analysis
 *
 * @param q Global configuration vector
 * @param qdot Global configuration first derivative vector
 * @return void
 **/
void FlexGlobalGalerkin::outputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    Body::outputStep();
    int stressVectorSize = (Simulation::getDim() == 2) ? 3 : 6;

    if (computeStress)
    {
        if (formulation == "LINEAR")
        {

            stresses.push_back(lmx::Vector<data_type>(stressVectorSize * nodes.size()));

            auto end_int = this->cells.size();
            for (auto i = 0u; i < end_int; ++i)
            {
                this->cells[i]->assembleRGaussPoints(stresses.back(), nodes[0]->getNumber());
            }
        }
        else if (formulation == "NONLINEAR")
        {

            stresses.push_back(lmx::Vector<data_type>(stressVectorSize * nodes.size()));

            auto end_int = this->cells.size();
            for (auto i = 0u; i < end_int; ++i)
            {
                this->cells[i]->assembleNLRGaussPoints(stresses.back(), nodes[0]->getNumber());
            }
        }
        recoverStressField(stressVectorSize);
    }

    if (computeEnergy)
    {
        energies.push_back(new lmx::Vector<data_type>(4)); //potential, kinetic, elastic, total

        energies.back()->fillIdentity(0.);

        auto end_int = this->cells.size();

        for (auto i = 0u; i < end_int; ++i)
        {
            energies.back()->operator()(0) += this->cells[i]->calcPotentialEGaussPoints(q); //potential
            energies.back()->operator()(1) += this->cells[i]->calcKineticEGaussPoints(qdot); //kinetic
            energies.back()->operator()(2) += this->cells[i]->calcElasticEGaussPoints(); //elastic
//total
        }
        energies.back()->operator()(3) +=
            energies.back()->readElement(0) + energies.back()->readElement(1) + energies.back()->readElement(2);
    }
}

/**
 * @brief Postprocess and store step results for static analysis
 *
 * @param q Global configuration vector
 * @return void
 **/
void FlexGlobalGalerkin::outputStep(const lmx::Vector<data_type>& q)
{
    Body::outputStep();
    int stressVectorSize = (Simulation::getDim() == 2) ? 3 : 6;

    if (computeStress)
    {
        if (formulation == "LINEAR")
        {

            stresses.push_back(lmx::Vector<data_type>(stressVectorSize * nodes.size()));

            auto end_int = this->cells.size();
//             #pragma omp parallel for
            for (auto i = 0u; i < end_int; ++i)
            {
                this->cells[i]->assembleRGaussPoints(stresses.back(), nodes[0]->getNumber());
            }
        }
        else if (formulation == "NONLINEAR")
        {

            stresses.push_back(lmx::Vector<data_type>(stressVectorSize * nodes.size()));

            auto end_int = this->cells.size();
//             #pragma omp parallel for
            for (auto i = 0u; i < end_int; ++i)
            {
                this->cells[i]->assembleNLRGaussPoints(stresses.back(), nodes[0]->getNumber());
            }
        }
        recoverStressField(stressVectorSize);
    }

    if (computeEnergy)
    {
        energies.push_back(new lmx::Vector<data_type>(4)); //potential, kinetic, elastic, total

        energies.back()->fillIdentity(0.);

        auto end_int = this->cells.size();
        #pragma omp parallel for
        for (auto i = 0u; i < end_int; ++i)
        {
            energies.back()->operator()(0) += this->cells[i]->calcPotentialEGaussPoints(q); //potential
            // no kinetic
            energies.back()->operator()(2) += this->cells[i]->calcElasticEGaussPoints(); //elastic
        }
        energies.back()->operator()(3) +=
            energies.back()->readElement(0) + energies.back()->readElement(1) + energies.back()->readElement(2); //total
    }
}

/**
 * @brief Performs stress field smoothing using global L2-projection if GLOBAL smoothing mode is active.
 * @param stressVectorSize Number of stress components per node (3 for 2D, 6 for 3D).
 */
void FlexGlobalGalerkin::recoverStressField(int stressVectorSize)
{
    if (Simulation::getSmoothingType() == "GLOBAL")
    {
        int j = 0;
        // Solve linear systems for stress smoothing:
        lmx::Vector<data_type> rhs(nodes.size());
        lmx::Vector<data_type> smooth(nodes.size());
//         for(j=0; j<stressVectorSize; ++j){
        for (auto i = 0u; i < nodes.size(); ++i)
        {
            rhs.writeElement(stresses.back().readElement(stressVectorSize * i + j), i);
        }
        lmx::LinearSystem<data_type> theLSolver(smoothingMassMatrix, smooth, rhs);
//           cout << smoothingMassMatrix << rhs << endl;
        theLSolver.solveYourself();
        for (auto i = 0u; i < nodes.size(); ++i)
        {
            stresses.back().writeElement(smooth.readElement(i), stressVectorSize * i + j);
        }
//           cout << smooth << endl;
    }
//     }
}

}

