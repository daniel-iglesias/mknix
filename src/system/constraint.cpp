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

#include "constraint.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor. Reads alpha and method from the global Simulation settings.
 */
Constraint::Constraint()
    : dim(Simulation::getDim())
    , iter_augmented(0)
    , alpha(Simulation::getAlpha())
    , augmentedTolerance(Simulation::getAugmentedTolerance())
    , method()
    , title()
{
}

/**
 * @brief Constructor with explicit alpha and method parameters.
 * @param alpha_in  Augmented-Lagrangian penalty parameter.
 * @param method_in Constraint enforcement method (e.g. "LAGRANGE", "AUGMENTED").
 */
Constraint::Constraint(double& alpha_in, std::string& method_in)
    : dim(Simulation::getDim())
    , iter_augmented(0)
    , alpha(alpha_in)
    , augmentedTolerance(Simulation::getAugmentedTolerance())
    , method(method_in)
    , title()
{
}

/**
 * @brief Constructor with explicit parameters including problem dimension.
 * @param alpha_in  Augmented-Lagrangian penalty parameter.
 * @param method_in Constraint enforcement method.
 * @param dim_in    Problem dimension (2 or 3).
 */
Constraint::Constraint(double& alpha_in, std::string& method_in, int dim_in)
    : dim(dim_in)
    , iter_augmented(0)
    , alpha(alpha_in)
    , augmentedTolerance(Simulation::getAugmentedTolerance())
    , method(method_in)
    , title()
{
}

/**
 * @brief Destructor.
 */
Constraint::~Constraint()
{
}


/**
 * @brief Writes the joint type, title and node numbers to the output file.
 * @param outfile Pointer to the output file stream.
 */
void Constraint::writeJointInfo(std::ofstream * outfile)
{
    int i, nodesSize(nodes.size());

    if (!title.empty())
    {
        *outfile << "\t" << method << " " << title << " ";
        for (i = 0; i < nodesSize; ++i)
        {
            *outfile << nodes[i]->getNumber() << " ";
        }
        *outfile << endl;
    }

}

/**
 * @brief Computes the local internal force vector from the current constraint residual
 *        and Lagrange multiplier estimates.
 */
void Constraint::calcInternalForces()
{
    internalForces.reset();
    calcPhi();
    calcPhiq();

    unsigned int i;
    for (i = 0; i < phi.size(); ++i)
    {
        internalForces += phi_q[i] * (alpha * phi[i] + lambda[i]);
//     cout << "i: " << i << ", phi_q[i]= "<< phi_q[i];
    }
}

/**
 * @brief Computes the local tangent stiffness matrix from the second derivatives of the constraint.
 */
void Constraint::calcTangentMatrix()
{
    stiffnessMatrix.reset();
    calcPhiqq();

    unsigned int i, j, k;
    for (i = 0; i < stiffnessMatrix.rows(); ++i)
    {
        for (j = 0; j < stiffnessMatrix.cols(); ++j)
        {
            for (k = 0; k < phi.size(); ++k)
            {
                stiffnessMatrix(i, j) += phi_qq[k](i, j) * (alpha * phi[k] + lambda[k]);
                stiffnessMatrix(i, j) += phi_q[k](i) * (alpha * phi_q[k](j));
            }
        }
    }
//   cout << "phi" << phi[0];
//   cout << "phi_q" << phi_q[0];
//   cout << "phi_qq" << phi_qq[0];
//    cout << "stiffnessMatrix: " << stiffnessMatrix << endl;
}


/**
 * @brief Assembles the local internal forces into the global force vector using shape-function mapping.
 * @param globalInternalForces Reference to the global internal force vector.
 */
void Constraint::assembleInternalForces
(lmx::Vector<data_type>& globalInternalForces)
{
    int nodesSize = nodes.size();
    int i, m;
    size_t k, counter(0);
    for (i = 0; i < nodesSize; ++i)
    {
        for (k = 0; k < nodes[i]->getSupportSize(0); ++k)
        {
            if (nodes[i]->getNumber() >= 0)
            {
                for (m = 0; m < Simulation::getDim(); ++m)
                {
//           cout << endl <<
//                   "Simulation::getDim()*nodes[i]->getSupportNodeNumber(0,k) + m = " <<
//                   Simulation::getDim()*nodes[i]->getSupportNodeNumber(0,k) + m << endl <<
//                   "Simulation::getDim()*counter + m = " <<
//                   Simulation::getDim()*counter + m << endl;
                    globalInternalForces(Simulation::getDim() * nodes[i]->getSupportNodeNumber(0, k) + m)
                    += internalForces.readElement(Simulation::getDim() * counter + m);
                }
            }
            ++counter;
        }
//        for (m=0; m<Simulation::getDim(); ++m){
//          globalInternalForces( Simulation::getDim()*nodes[i]->getNumber() + m)
//            += internalForces.readElement(Simulation::getDim()*i + m);
    }
}

/**
 * @brief Assembles the local tangent matrix into the global tangent matrix using shape-function mapping.
 * @param globalTangent Reference to the global tangent (stiffness) matrix.
 */
void Constraint::assembleTangentMatrix
(lmx::Matrix<data_type>& globalTangent)
{
    int nodesSize = nodes.size();
    size_t supportNodesSize_i, supportNodesSize_j;
    int i, j, m, n;
    size_t k, l, counter_i(0), counter_j(0);
    for (i = 0; i < nodesSize; ++i)
    {
        supportNodesSize_i = nodes[i]->getSupportSize(0);
        for (k = 0; k < supportNodesSize_i; ++k)
        {
            for (j = 0; j < nodesSize; ++j)
            {
                supportNodesSize_j = nodes[j]->getSupportSize(0);
                for (l = 0; l < supportNodesSize_j; ++l)
                {
                    if (nodes[i]->getNumber() >= 0)
                    {
                        if (nodes[j]->getNumber() >= 0)
                        {
                            for (m = 0; m < Simulation::getDim(); ++m)
                            {
                                for (n = 0; n < Simulation::getDim(); ++n)
                                {
                                    globalTangent(Simulation::getDim() * nodes[i]->getSupportNodeNumber(0, k) + m,
                                                  Simulation::getDim() * nodes[j]->getSupportNodeNumber(0, l) + n)
                                    += stiffnessMatrix.readElement(Simulation::getDim() * counter_i + m,
                                                                   Simulation::getDim() * counter_j + n);
                                }
                            }
                        }
                    }
                    ++counter_j;
                }
//              for ( m=0; m<Simulation::getDim(); ++m ){
//                for ( n=0; n<Simulation::getDim(); ++n ){
//                  globalTangent(Simulation::getDim()*nodes[i]->getNumber() + m,
//                                Simulation::getDim()*nodes[j]->getNumber() + n)
//                  += stiffnessMatrix.readElement(Simulation::getDim()*i + m,
//                                                 Simulation::getDim()*j + n);
//                }
//              }
            }
            counter_j = 0;
            ++counter_i;
        }
    }
}

/**
 * @brief Checks whether the AUGMENTED method has converged. Updates Lagrange multipliers if not.
 * @return True if converged or method is not AUGMENTED; false otherwise.
 */
bool Constraint::checkAugmented()
{
    if (method == "AUGMENTED")
    {
//         double energy(0);
        double delta(0);
        for (auto i = 0u; i < phi.size(); ++i)
        {
            delta += std::abs(phi[i]);
//             energy += phi[i] * phi[i];
//             cout << endl << "Phi: " << phi[0] << "\t\t\t\t";
        }
//         energy *= 0.5 * alpha;
//         if (energy <= 2E5) {
        if (delta <= augmentedTolerance)
        {
//             cout << endl << "Energy: " << energy << "\t\t\t\t";
            return 1;
        }
        else
        {
//             cout << endl << "energy: " << energy << "\t\t\t\t";
            for (auto i = 0u; i < phi.size(); ++i)
            {
                lambda[i] += alpha * phi[i];
            }
            return 0;
        }
    }
    else
    {
        return 1;
    }
}

}

/**
 * @brief Resets the Lagrange multipliers to zero for the next augmented-Lagrangian cycle.
 */
void mknix::Constraint::clearAugmented()
{
//   if( nodes[0]->getNumber() < 0 || nodes[0]->getNumber() < 0 ){
    /*  cout << "Convergence in constraint... NODEA: " << endl
          << nodes[0]->getNumber() << ","
          << nodes[0]->getx() << ","
          << nodes[0]->gety() << ","
          << nodes[0]->getz() << ","
          << "NODEB: " << endl
          << nodes[1]->getNumber() << ","
          << nodes[1]->getx() << ","
          << nodes[1]->gety() << ","
          << nodes[1]->getz() << ","
          << endl;
    //   }*/

    if (method == "AUGMENTED")
    {
        for (auto i = 0u; i < phi.size(); ++i)
        {
            lambda[i] = 0.;
        }
    }
}

/**
 * @brief Stores the current internal force vector for post-processing (dynamic step).
 * @param q    Global configuration vector (unused here).
 * @param qdot Global velocity vector (unused here).
 */
void mknix::Constraint::outputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    if (!title.empty())   // do not store internal constraints of rigid bodies
    {
        internalForcesOutput.push_back(lmx::Vector<data_type>(internalForces.size()));
        internalForcesOutput.back() = internalForces;
    }
}

/**
 * @brief Stores the current internal force vector for post-processing (static step).
 * @param q Global configuration vector (unused here).
 */
void mknix::Constraint::outputStep(const lmx::Vector<data_type>& q)
{
    if (!title.empty())   // do not store internal constraints of rigid bodies
    {
        internalForcesOutput.push_back(lmx::Vector<data_type>(internalForces.size()));
        internalForcesOutput.back() = internalForces;
    }
}

/**
 * @brief Writes stored constraint force history to the output file.
 * @param outFile Pointer to the output file stream.
 */
void mknix::Constraint::outputToFile(std::ofstream * outFile)
{
    auto vectorSize = internalForces.size();

    if (internalForcesOutput.size() > 0)
    {

        *outFile << "FORCES " << title << " " << vectorSize << endl;
        for (auto& force : internalForcesOutput)
        {
            for (auto i = 0u; i < vectorSize; ++i)
            {
                *outFile << force.readElement(i) << " ";
            }
            *outFile << endl;
        }
    }

}