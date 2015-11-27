/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   http://code.google.com/p/mknix                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "bodyflexglobalgalerkin.h"
#include "simulation.h"
#include "cell.h"
#include "node.h"

namespace mknix {

FlexGlobalGalerkin::FlexGlobalGalerkin()
    : FlexBody()
{
}


FlexGlobalGalerkin::FlexGlobalGalerkin( std::string title_in )
    : FlexBody( title_in )
{
}


FlexGlobalGalerkin::~FlexGlobalGalerkin()
{
}

void FlexGlobalGalerkin::calcMassMatrix()
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->computeMGaussPoints( );
        this->cells[i]->computeFextGaussPoints( ); // Depends on M and gravity
    }

}

void FlexGlobalGalerkin::calcInternalForces()
{
    int end_int = this->cells.size();
    if( formulation == "NONLINEAR" ) {
//         #pragma omp parallel for
        for (int i=0;
                i < end_int;
                ++i)
        {
            this->cells[i]->computeNLFintGaussPoints( );
        }
    }
    else if( formulation == "LINEAR" ) {
//         #pragma omp parallel for
        for (int i=0;
                i < end_int;
                ++i)
        {
            this->cells[i]->computeFintGaussPoints( );
        }
    }

}

void FlexGlobalGalerkin::calcExternalForces()
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i=0;
            i < end_int;
            ++i)
    { // Moved to calcMassMatrix(), as M and Gravity are constant
//         this->cells[i]->computeFextGaussPoints( );
    }

}

void FlexGlobalGalerkin::calcTangentMatrix()
{
    int end_int = this->cells.size();
    if( formulation == "NONLINEAR" ) {
        #pragma omp parallel for
        for (int i=0;
                i < end_int;
                ++i)
        {
            this->cells[i]->computeNLKGaussPoints( );
        }
    }
    else if( formulation == "LINEAR" ) {
        #pragma omp parallel for
        for (int i=0;
                i < end_int;
                ++i)
        {
            this->cells[i]->computeKGaussPoints( );
        }
    }

}

void FlexGlobalGalerkin::assembleMassMatrix
(lmx::Matrix< data_type > & globalMass)
{
    int end_int = this->cells.size();
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->assembleMGaussPoints( globalMass );
    }

    // Chapuza...
    if( Simulation::getSmoothingType() == "GLOBAL") {
        int i,j;
        if( smoothingMassMatrix.rows() == 0 ) {
            // Compute smoothing matrix (local mass matrix in 1D)
            smoothingMassMatrix.resize( nodes.size(), nodes.size() );
            for(i=0; i<nodes.size(); ++i) {
                for(j=0; j<nodes.size(); ++j) {
                    smoothingMassMatrix.writeElement( globalMass.readElement
                                                      ( Simulation::getDim() * (nodes[i]->getNumber()),
                                                        Simulation::getDim() * (nodes[j]->getNumber()) ),
                                                      i, j );
                }
            }
        }
    }
}

void FlexGlobalGalerkin::assembleInternalForces
(lmx::Vector< data_type > & globalInternalForces)
{
    int end_int = this->cells.size();
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->assembleFintGaussPoints( globalInternalForces );
    }

}

void FlexGlobalGalerkin::assembleExternalForces
(lmx::Vector< data_type > & globalExternalForces)
{
    int end_int = this->cells.size();
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->assembleFextGaussPoints( globalExternalForces );
    }

}

void FlexGlobalGalerkin::assembleTangentMatrix(lmx::Matrix< data_type > & globalTangent)
{
    int end_int = this->cells.size();
    for (int i=0;
            i < end_int;
            ++i)
    {
        this->cells[i]->assembleKGaussPoints( globalTangent );
    }

}

/**
 * @brief Postprocess and store step results for dynamic analysis
 *
 * @param q Global configuration vector
 * @param qdot Global configuration first derivative vector
 * @return void
 **/
void FlexGlobalGalerkin::outputStep( const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot )
{
    Body::outputStep();
    int stressVectorSize = (Simulation::getDim()==2) ? 3 : 6;

    if( computeStress ) {
        if ( formulation=="LINEAR" ) {

            stress.push_back( lmx::Vector<data_type>( stressVectorSize*nodes.size() ) );

            int end_int = this->cells.size();
            for (int i=0;
                    i < end_int;
                    ++i)
            {
                this->cells[i]->assembleRGaussPoints( stress.back(), nodes[0]->getNumber() );
            }
        }
        else if( formulation=="NONLINEAR" ) {

            stress.push_back( lmx::Vector<data_type>( stressVectorSize*nodes.size() ) );

            int end_int = this->cells.size();
            for (int i=0;
                    i < end_int;
                    ++i)
            {
                this->cells[i]->assembleNLRGaussPoints( stress.back(), nodes[0]->getNumber() );
            }
        }
        recoverStressField(stressVectorSize);
    }

    if( computeEnergy ) {
        energy.push_back( new lmx::Vector<data_type>( 4 ) ); //potential, kinetic, elastic, total

        energy.back()->fillIdentity( 0. );

        int end_int = this->cells.size();

        for (int i=0;
                i < end_int;
                ++i)
        {
            energy.back()->operator()(0) += this->cells[i]->calcPotentialEGaussPoints( q ); //potential
            energy.back()->operator()(1) += this->cells[i]->calcKineticEGaussPoints( qdot ); //kinetic
            energy.back()->operator()(2) += this->cells[i]->calcElasticEGaussPoints( ); //elastic
//total
        }
        energy.back()->operator()(3) += energy.back()->readElement(0) + energy.back()->readElement(1) + energy.back()->readElement(2);
    }
}

/**
 * @brief Postprocess and store step results for static analysis
 *
 * @param q Global configuration vector
 * @return void
 **/
void FlexGlobalGalerkin::outputStep( const lmx::Vector<data_type>& q )
{
    Body::outputStep();
    int stressVectorSize = (Simulation::getDim()==2) ? 3 : 6;

    if( computeStress ) {
        if ( formulation=="LINEAR" ) {

            stress.push_back( lmx::Vector<data_type>( stressVectorSize*nodes.size() ) );

            int end_int = this->cells.size();
//             #pragma omp parallel for
            for (int i=0;
                    i < end_int;
                    ++i)
            {
                this->cells[i]->assembleRGaussPoints( stress.back(), nodes[0]->getNumber() );
            }
        }
        else if( formulation=="NONLINEAR" ) {

            stress.push_back( lmx::Vector<data_type>( stressVectorSize*nodes.size() ) );

            int end_int = this->cells.size();
//             #pragma omp parallel for
            for (int i=0;
                    i < end_int;
                    ++i)
            {
                this->cells[i]->assembleNLRGaussPoints( stress.back(), nodes[0]->getNumber() );
            }
        }
        recoverStressField(stressVectorSize);
    }

    if( computeEnergy ) {
        energy.push_back( new lmx::Vector<data_type>( 4 ) ); //potential, kinetic, elastic, total

        energy.back()->fillIdentity( 0. );

        int end_int = this->cells.size();
        #pragma omp parallel for
        for (int i=0;
                i < end_int;
                ++i)
        {
            energy.back()->operator()(0) += this->cells[i]->calcPotentialEGaussPoints( q ); //potential
            // no kinetic
            energy.back()->operator()(2) += this->cells[i]->calcElasticEGaussPoints( ); //elastic
        }
        energy.back()->operator()(3) += energy.back()->readElement(0) + energy.back()->readElement(1) + energy.back()->readElement(2); //total
    }
}

void FlexGlobalGalerkin::recoverStressField(int stressVectorSize)
{
    if( Simulation::getSmoothingType() == "GLOBAL") {
        int i,j(0);
        // Solve linear systems for stress smoothing:
        lmx::Vector<data_type> rhs(nodes.size());
        lmx::Vector<data_type> smooth(nodes.size());
//         for(j=0; j<stressVectorSize; ++j){
        for(i=0; i<nodes.size(); ++i)
            rhs.writeElement( stress.back().readElement(stressVectorSize*i+j), i );
        lmx::LinearSystem<data_type> theLSolver( smoothingMassMatrix, smooth, rhs );
//           cout << smoothingMassMatrix << rhs << endl;
        theLSolver.solveYourself();
        for(i=0; i<nodes.size(); ++i)
            stress.back().writeElement( smooth.readElement(i), stressVectorSize*i+j );
//           cout << smooth << endl;
    }
//     }
}

}

