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
#include <gpu/chTimer.h>
#include <gpu/assembly_cpu.h>
#ifdef HAVE_CUDA
#include <gpu/assembly_kernels.h>
#include <gpu/cuda_helper.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif



namespace mknix {

Body::Body()
        : computeEnergy(0)
        , isThermal(1)
{
  _use_gpu = false;
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
  _use_gpu = false;
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
    ///////////////////////////////////
    free(_h_presence_matrix);
    free(_h_globalCapacity);
    free(_h_globalConductivity);
    ////////////////////////////////////
  double cpuAvg = 0.0;
  int cpuMeasures = microCPU1.size();
  for(int i = 0; i < cpuMeasures; i++)
        cpuAvg += microCPU1.at(i);
  cpuAvg /= cpuMeasures;

  std::cout << "AssembleCapacityMatrix time measurements = "<< cpuMeasures << std::endl;
  std::cout << "CPU Avg:  " << cpuAvg << std::endl;
  //std::cout << "CPU Avg:  " << cpuAvg << " vs GPU Avg: " << gpuAvg
  cpuMeasures = microCPU2.size();
  cpuAvg = 0.0;
  for(int i = 0; i < cpuMeasures; i++) cpuAvg += microCPU2.at(i);
  cpuAvg /= cpuMeasures;

  std::cout << "AssembleConductivityMatrix time measurements = "<< cpuMeasures << std::endl;


#ifdef HAVE_CUDA
    if(_use_gpu){
      cudaFree(_d_globalCapacityf);
      cudaFree(_d_globalConductivityf);
      cudaFree(_d_capacity_map);
    }
    int gpuMeasures = microGPU1.size();

    if(cpuMeasures != gpuMeasures)
      std::cout << "WARNING we have " << cpuMeasures << " CPU measures and " << gpuMeasures << " GPU ones, the should be the same!"<< std::endl;

   double gpuAvg = 0.0;
   for(int i = 0; i < gpuMeasures; i++)
         gpuAvg += microGPU1.at(i);
   gpuAvg /= gpuMeasures;
   std::cout << "CPU Avg:  " << cpuAvg << " vs GPU Avg: " << gpuAvg << std::endl;
   std::cout << "Represent and average of  " << cpuAvg/gpuAvg << "x speedup" << std::endl;

   gpuMeasures = microGPU2.size();
   if(cpuMeasures != gpuMeasures)
     std::cout << "WARNING we have " << cpuMeasures << " CPU measures and " << gpuMeasures << " GPU ones, the should be the same!"<< std::endl;

   gpuAvg = 0.0;
   for(int i = 0; i < gpuMeasures; i++)
         gpuAvg += microGPU2.at(i);
   gpuAvg /= gpuMeasures;

    std::cout << "Represent and average of  " << cpuAvg/gpuAvg << "x speedup" << std::endl;
#endif

////////////////////////////////////////////////////////////////////////




  /////////////////////////////////////////////////////////////////////////////////////
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
    _number_points = end_int;
    //suppossition all nodes have same supportNodesSize
    _support_node_size = 4;
    std::cout << "last node " << lastNode << std::endl;
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
    _number_nodes = end_int;
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

    _h_presence_matrix = (int*) malloc(_number_nodes * _number_nodes * sizeof(int));
    _locaThermalNumbers.resize(_number_points * _support_node_size );

    mapConectivityCapacityMatrix();
    map_global_matrix(_full_map,
                      _vec_ind,
                      _cvec_ptr,
                      _h_presence_matrix,
                      _number_nodes,
                      _number_nodes);

    //map_vector_thermal_numbers(_locaThermalNumbers);
    //GPU part
  #ifdef HAVE_CUDA
    if(_use_gpu){
       _sparse_matrix_size = _cvec_ptr[_number_nodes];//convention
       std::cout << "sparse matrix has " << _sparse_matrix_size << " elements" << std::endl;
       CudaSafeCall(cudaMalloc((void**)&_d_globalCapacityf, _sparse_matrix_size * sizeof(float)));
       CudaSafeCall(cudaMalloc((void**)&_d_globalConductivityf, _sparse_matrix_size * sizeof(float)));
       CudaSafeCall(cudaMalloc((void**)&_d_localCapacityf, _number_points * _support_node_size * _support_node_size * sizeof(float)));
       CudaSafeCall(cudaMalloc((void**)&_d_localConductivityf, _number_points * _support_node_size * _support_node_size * sizeof(float)));
       CudaSafeCall(cudaMalloc((void**)&_d_capacity_map, _number_nodes * _number_nodes*sizeof(int)));
       CudaSafeCall(cudaMemcpy(_d_capacity_map, _full_map.data(),_number_nodes * _number_nodes*sizeof(int), cudaMemcpyHostToDevice));
       CudaCheckError();
    }
 #endif
 _h_globalCapacity = (data_type*)malloc(_sparse_matrix_size * sizeof(data_type));
 _h_globalConductivity = (data_type*)malloc(_sparse_matrix_size * sizeof(data_type));
}

/**
 * @brief Computes the local Capacity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::calcCapacityMatrix()
{
    auto end_int = this->cells.size();
//#pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->computeCapacityGaussPoints();
    }
}

/**
 * @brief Computes the local Capacity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::mapConectivityCapacityMatrix()
{

    auto end_int = this->cells.size();
    std::cout << "Mapping capacity matrix with "<< end_int << " cells" <<std::endl;
//#pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->presenceCapacityGaussPoints(_h_presence_matrix, _number_nodes);
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
//#pragma omp parallel for
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
  cpuClock cck1;
  cpuTick(&cck1);
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleCapacityGaussPoints(globalCapacity);
    }
  cpuTock(&cck1, "CPU assembleCapacityMatrix");
  microCPU1.push_back(cck1.elapsedMicroseconds);
  //new CPU assembly function
  cpuClock cck1b;
  cpuTick(&cck1b);
    init_host_array_to_value(_h_globalCapacity, 0.0, _sparse_matrix_size);
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleCapacityGaussPointsWithMap(_h_globalCapacity,
                                                           _full_map.data(),
                                                           _support_node_size);
    }
  cpuTock(&cck1b, "CPU assembleCapacityMatrixWithMap");
  microCPU1b.push_back(cck1b.elapsedMicroseconds);
  //std::cout << "miiiii " << cck1.elapsedMicroseconds <<std::endl;
#ifdef HAVE_CUDA
cudaClock gck1;
 if(_use_gpu){
   init_array_to_value(_d_globalCapacityf, 0.0f, _sparse_matrix_size,128);
   cudaTick(&gck1);
     gpu_assemble_global_matrix(_d_globalCapacityf,
                                _d_capacity_map,
                                _d_localCapacityf,
                                this->cells.size(),
                                _support_node_size,
                                _number_points,
                                128);
   cudaTock(&gck1, "GPU assembleCapacityMatrix");
   microGPU1.push_back(gck1.elapsedMicroseconds);
 }
 #endif

}

/**
 * @brief Assembles the local conductivity into the global matrix by calling each cell's cascade function.
 *
 * @param globalConductivity Reference to the global matrix of the thermal simulation.
 * @return void
 **/
void Body::assembleConductivityMatrix(lmx::Matrix<data_type>& globalConductivity)
{ cpuClock cck2;
    cpuTick(&cck2);
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleConductivityGaussPoints(globalConductivity);
    }
    cpuTock(&cck2, "CPU assembleConductivityMatrix");
    microCPU2.push_back(cck2.elapsedMicroseconds);
    //new CPU assembly function
    cpuClock cck2b;
    cpuTick(&cck2b);
    init_host_array_to_value(_h_globalConductivity, 0.0, _sparse_matrix_size);
      for (auto i = 0u; i < end_int; ++i) {
          this->cells[i]->assembleConductivityGaussPointsWithMap(_h_globalConductivity,
                                                                 _full_map.data(),
                                                                 _support_node_size);
      }
    cpuTock(&cck2b, "CPU assembleConductivityGaussPointsWithMap");
    microCPU2b.push_back(cck2b.elapsedMicroseconds);
  #ifdef HAVE_CUDA
    cudaClock gck2;
    if(_use_gpu){
      init_array_to_value(_d_globalConductivityf, 0.0f, _sparse_matrix_size,128);
      cudaTick(&gck2);
      gpu_assemble_global_matrix(_d_globalConductivityf,
                                _d_capacity_map,
                                _d_localConductivityf,
                                this->cells.size(),
                                _support_node_size,
                                _number_points,
                                128);
    cudaTock(&gck2, "GPU assembleConductivityMatrix");
    microGPU2.push_back(gck2.elapsedMicroseconds);
    }
  #endif

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
  #ifdef HAVE_CUDA
    if(_use_gpu){
      init_array_to_value(_d_globalExternalHeat, 0.0f, _sparse_matrix_size,128);
      gpu_assemble_global_vector(_d_globalExternalHeat,
                                _d_locaThermalNumbers,
                                _d_localHeatf,
                                _support_node_size,
                                _number_points,
                                128);
      gpu_assemble_global_vector(_d_globalExternalHeat,
                                _d_locaThermalNumbers,
                                _d_localHeatf,
                                _support_node_size,
                                _number_points,
                                128);
    }
  #endif
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
    std::vector<std::vector<int> >::iterator itBoundary;
    std::vector<int>::iterator itOneBoundarySegment;
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
