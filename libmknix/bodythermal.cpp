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
#include "bodythermal.h"
#include "simulation.h"
#include "node.h"
#include "cell.h"


namespace mknix {

ThermalBody::ThermalBody()
        : computeEnergy(0)
//   , formulation( "NONLINEAR" )
{
}


/**
 * @brief Constructor with 1 parameter
 *
 * @param title_in Name of body in the system. Will be the same as the associated material body
 **/
ThermalBody::ThermalBody(std::string title_in)
        : title(title_in)
        , computeEnergy(0)
//   , formulation( "NONLINEAR" )
{
}


ThermalBody::~ThermalBody()
{
}

/**
 * @brief Cascade initialization funtion. Calls the initialize methods for each of the Cells
 *        and tells them to compute their shapefunction values. Both loops are parallelized.
 *
 * @return void
 **/
void ThermalBody::initialize()
{
    int end_int = this->cells.size();

//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->initialize(this->nodes);
    }

//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->computeShapeFunctions();
    }

// Checking the output of a shapefunction:
//   int mid_int = this->cells.size()/2;
//   std::ofstream cell_data("cell_data.dat");
//   std::ofstream gpoint_data("cell_gpoint_data.dat");
//
//   this->cells[mid_int]->gnuplotOut(cell_data, gpoint_data); // Bus error

    // Initialize output files
}

/**
 * @brief Computes the local Capacity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void ThermalBody::calcCapacityMatrix()
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->computeCapacityGaussPoints();
    }
}

/**
 * @brief Computes the local Conductivity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void ThermalBody::calcConductivityMatrix()
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->computeConductivityGaussPoints();
    }
}

/**
 * @brief Computes the local volumetric heat vector of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void ThermalBody::calcExternalHeat()
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->computeQextGaussPoints(this->loadThermalBody);
    }

}

/**
 * @brief Assembles the local conductivity into the global matrix by calling each cell's cascade function.
 *
 * @param globalCapacity Reference to the global matrix of the thermal simulation.
 * @return void
 **/
void ThermalBody::assembleCapacityMatrix(lmx::Matrix<data_type>& globalCapacity)
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->assembleCapacityGaussPoints(globalCapacity);
    }
}

/**
 * @brief Assembles the local conductivity into the global matrix by calling each cell's cascade function.
 *
 * @param globalConductivity Reference to the global matrix of the thermal simulation.
 * @return void
 **/
void ThermalBody::assembleConductivityMatrix(lmx::Matrix<data_type>& globalConductivity)
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->assembleConductivityGaussPoints(globalConductivity);
    }
}

/**
 * @brief Assembles the local volumetric heat into the global heat load vector by calling each cell's cascade function.
 *
 * @return void
 **/
void ThermalBody::assembleExternalHeat(lmx::Vector<data_type>& globalExternalHeat)
{
    int end_int = this->cells.size();
//     #pragma omp parallel for
    for (int i = 0;
         i < end_int;
         ++i) {
        this->cells[i]->assembleQextGaussPoints(globalExternalHeat);
    }
}


/**
 * @brief Activates a flag for output data at the end of the analysis.
 *
 * @see outputToFile()
 * @see outputStep()
 * @param outputType_in Keyword of the flag. Options are: [ENERGY]
 * @return void
 **/
void ThermalBody::setOutput(std::string outputType_in)
{
    if (outputType_in == "ENERGY") {
        computeEnergy = 1;
    }
}


/**
 * @brief Streams the data stored during the analysis to a file.
 * The idea is that all thermalBodies will be linked to a solid body with the same name (title).
 * The postprocessor will treat this values as if it was the same body, so it must take into account
 * that it can need to read several ENERGY(.*) sections for the same body.
 *
 * @param outFile Output files
 * @return void
 **/
void ThermalBody::outputToFile(std::ofstream * outFile)
{
//   if( computeEnergy ){
//     std::vector< lmx::Vector<data_type>* >::iterator itEnergy;
//     int i, vectorSize;
//
//     *outFile << "ENERGY.THERMAL " << title << endl;
//     for( itEnergy = energy.begin();
//          itEnergy!= energy.end();
//          ++itEnergy
//        )
//     {
//       vectorSize = (*itEnergy)->size();
//       for( i=0; i<vectorSize; ++i){
//         *outFile << (*itEnergy)->readElement(i) << " ";
//       }
//       *outFile << endl;
//     }
//   }

    if (temperature.size() != 0) {
        *outFile << "TEMPERATURE " << title << endl;
        for (auto& temp : temperature) {
            auto vectorSize = temp->size();
            for (auto i = 0u; i < vectorSize; ++i) {
                *outFile << temp->readElement(i) << " ";
            }
            *outFile << endl;
        }
    }
}

/**
 * @brief Postprocess and store step results for dynamic analysis
 *
 * @param q Global configuration vector
 * @param qdot Global configuration first derivative vector
 * @return void
 **/
void ThermalBody::outputStep(const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    if (computeEnergy) { // TODO: store thermal energy
//     energy.push_back( new lmx::Vector<data_type>( 4 ) ); //potential, kinetic, elastic, total
//
//     energy.back()->fillIdentity( 0. );
//
//     int end_int = this->cells.size();
// #pragma omp parallel for
//     for (int i=0;
//          i < end_int;
//          ++i)
//     {
//       energy.back()->operator()(0) += this->cells[i]->calcPotentialEGaussPoints( q ); //potential
//       energy.back()->operator()(1) += this->cells[i]->calcKineticEGaussPoints( qdot ); //kinetic
//       energy.back()->operator()(2) += this->cells[i]->calcElasticEGaussPoints( ); //elastic
// //total
//     }
//     energy.back()->operator()(3) += energy.back()->readElement(0) + energy.back()->readElement(1) + energy.back()->readElement(2);
    }
}

/**
 * @brief Postprocess and store step results for static analysis
 *
 * @param q Global configuration vector
 * @return void
 **/
void ThermalBody::outputStep(const lmx::Vector<data_type>& q)
{

    if (computeEnergy) { // TODO: see above
    }
}

}
