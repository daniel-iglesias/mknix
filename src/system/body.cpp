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
#include "body.h"

#include <core/cell.h>

namespace mknix {

Body::Body()
        : computeEnergy(0)
        , isThermal(1)
{
}

/**
 * @brief Constructor with 1 parameter
 *
 * @param title_in Name of body in the system. Will be the same as the associated material body
 **/
Body::Body(std::string title_in)
        : title(title_in)
        , lastNode(0)
        , computeEnergy(0)
        , isThermal(1)
{
}


Body::~Body()
{
    for (auto& temp : temperature) {
        delete temp;
    }
    for (auto& cell : cells) {
        delete cell.second;
    }
    /*
    for (auto& node : nodes) {
        delete node;
    }
    for (auto& node : bondedBodyNodes) {
        delete node;
    }
    */
    for (auto& group : boundaryGroups) {
        delete group.second;
    }
}

/**
 * @brief Cascade initialization funtion. Calls the initialize methods for each of the Cells
 *        and tells them to compute their shapefunction values. Both loops are parallelized.
 *
 * @return void
 **/
void Body::initialize()
{
    lastNode = nodes.back();
    auto end_int = this->cells.size();

    nodes.insert(nodes.end(), bondedBodyNodes.begin(), bondedBodyNodes.end());

//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->initialize(this->nodes);
    }

//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->computeShapeFunctions();
    }

//  // Checking the output of a shapefunction:
//   int mid_int = this->cells.size()/2;
//   // Initialize individual output files
//   std::ofstream cell_data(std::string("cell_data_"+title+".dat").c_str());
//   std::ofstream gpoint_data(std::string("cell_gpoint_data_"+title+".dat").c_str());
//   this->cells[mid_int]->gnuplotOut(cell_data, gpoint_data); // Bus error

//The iteration on nodes MUST be done AFTER the cells.
    end_int = this->nodes.size();

//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        if (this->nodes[i]->getShapeFunType() == "RBF" ||
            this->nodes[i]->getShapeFunType() == "MLS") {
                this->nodes[i]->findSupportNodes(this->nodes);
        }
    }

//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        if (this->nodes[i]->getShapeFunType() == "RBF") {
            this->nodes[i]->shapeFunSolve("RBF", 1.03);
        }
        if (this->nodes[i]->getShapeFunType() == "MLS") {
            this->nodes[i]->shapeFunSolve("MLS", 1.03);
        }
    }
    std::map<std::string, BoundaryGroup *>::iterator it_boundaryGroups;
    for (auto& group : boundaryGroups) {
        group.second->initialize();
    }
}

/**
 * @brief Computes the local Capacity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::calcCapacityMatrix()
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->computeCapacityGaussPoints();
    }
}

/**
 * @brief Computes the local Conductivity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::calcConductivityMatrix()
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->computeConductivityGaussPoints();
    }
}

/**
 * @brief Computes the local volumetric heat vector of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::calcExternalHeat()
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    if (loadThermalBody) {
        for (auto i = 0u; i < end_int; ++i) {
            this->cells[i]->computeQextGaussPoints(this->loadThermalBody);
        }
    }
    for (auto group : boundaryGroups) {
        group.second->calcExternalHeat();
    }
}

/**
 * @brief Assembles the local conductivity into the global matrix by calling each cell's cascade function.
 *
 * @param globalCapacity Reference to the global matrix of the thermal simulation.
 * @return void
 **/
void Body::assembleCapacityMatrix(lmx::Matrix<data_type>& globalCapacity)
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleCapacityGaussPoints(globalCapacity);
    }
}

/**
 * @brief Assembles the local conductivity into the global matrix by calling each cell's cascade function.
 *
 * @param globalConductivity Reference to the global matrix of the thermal simulation.
 * @return void
 **/
void Body::assembleConductivityMatrix(lmx::Matrix<data_type>& globalConductivity)
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleConductivityGaussPoints(globalConductivity);
    }
}

/**
 * @brief Assembles the local volumetric heat into the global heat load vector by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::assembleExternalHeat(lmx::Vector<data_type>& globalExternalHeat)
{
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleQextGaussPoints(globalExternalHeat);
    }
    for (auto group : boundaryGroups) {
        group.second->assembleExternalHeat(globalExternalHeat);
    }
//     cout << globalExternalHeat << endl;
}

void Body::setTemperature(double temp_in)
{
    for (auto& node : nodes) {
        node->setqt(temp_in);
    }
}


/**
 * @brief Postprocess and store thermal step results for any analysis
 *
 * @return void
 **/
void Body::outputStep()
{
    if (isThermal && nodes.size() != 0) {
        temperature.push_back(new lmx::Vector<data_type>(nodes.size())); //temperature
        for (auto i = 0u; i < nodes.size(); ++i) {
            temperature.back()->writeElement(nodes[i]->getqt(), i);
        }
    }

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
 * @brief Streams the data stored during the analysis to a file.
 *
 * @param outFile Output files
 * @return void
 **/
void Body::outputToFile(std::ofstream * outFile)
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

    if (boundaryConnectivity.size() > 0) {
        *outFile << "BOUNDARY " << title << " " << boundaryConnectivity.size() << endl;
        for (auto& boundary : boundaryConnectivity) {
            for (auto& segment : boundary) {
                *outFile << segment << " ";
            }
            *outFile << endl;
        }
    }

    if (temperature.size() != 0) {
        *outFile << "TEMPERATURE " << title << endl;
        for (auto& temp : temperature) {
            auto  vectorSize = temp->size();
            for (auto i = 0u; i < vectorSize; ++i) {
                *outFile << temp->readElement(i) << " ";
            }
            *outFile << endl;
        }
    }
}

void Body::addBoundaryConnectivity(std::vector<int> connectivity_in)
{
    this->boundaryConnectivity.push_back(std::vector<int>(connectivity_in));
}


void Body::translate(double x_in, double y_in, double z_in)
{
    for (auto node : nodes) {
        node->setX(node->getX() + x_in);
        node->setY(node->getY() + y_in);
        node->setZ(node->getZ() + z_in);
    }
}

}
