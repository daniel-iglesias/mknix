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
#include <cpu/chTimer.h>
#include <cpu/cpu_run_type.h>
#include <cpu/assembly_cpu.h>

#include <omp.h>

//#ifdef HAVE_CUDA
#include "gpu/assembly_kernels.h"
#include "gpu/cuda_helper.h"
//#include <gpu/device_helper.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
//#endif

#define DEBUG_CELL 1121

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
    ///////////////////////////////////
    //freeMaterialTableMemory(&_h_materials);
    free(_h_materials);
  //  free_shape_functions_table(&_h_shapeFunctionTable);
    free(_h_shapeFunctionTable);
  //  freeThermalBoundaryTableMemory(&_h_thermal_boundaries);
    free(_h_thermal_boundaries);

    free(_h_materials_cap_ids);
    free(_h_materials_cond_ids);
    free(_h_presence_matrix_cap);
    free(_h_presence_matrix_cond);
    free(_h_localCapacityf);
    //free(_h_globalConductivity);
    free(_h_localConductivityf);
    free(_h_local_temperatures_cap_array);
    free(_h_local_jacobian_cap_array);
    free(_h_local_weight_cap_array);
    free(_h_local_temperatures_cond_array);
    free(_h_local_jacobian_cond_array);
    free(_h_local_weight_cond_array);
    if(MULTICPU){
    free(_param_array_capacity);
    free(_param_array_conductivity);
    free(_param_calc_capacity);
    free(_param_calc_conductivity);
  }
    ////////////////////////////////////
  if(OLD_CODE){
    clockFullStats(microCPU_old_conductivity, "Existing CPU calcConductivityMatrix");
    clockFullStats(microCPU_old_capacity,     "Existing CPU calcCapacityMatrix");
    clockFullStats(microCPU1,                 "Existing AssembleCapacityMatrix");
    clockFullStats(microCPU2,                 "ExistingAssembleConductivityMatrix");
  }
  if(NEWCPU){
    clockFullStats(microCPU_single_conductivity, "New Single CPU calcConductivityMatrix");
    clockFullStats(microCPU_single_capacity,     "New Single CPU calcCapacityMatrix");
    clockFullStats(microCPU1b,                   "New AssembleCapacityMatrixWithMap");
    clockFullStats(microCPU2b,                   "New AssembleConductivityMatrixWithMap");
  }
 if(MULTICPU){
   clockFullStats(microCPU_multi_conductivity,  "Multi-CPU calcConductivityMatrix");
   clockFullStats(microCPU_multi_capacity,      "Multi-CPU calcCapacityMatrix");
   clockFullStats(microCPU1c,                   "Multi-CPU AssembleCapacityMatrixWithMap");
   clockFullStats(microCPU2c,                   "Multi-CPU AssembleConductivityMatrixWithMap");
 }

#ifdef HAVE_CUDA
    if(GPU){
      /*freeGPU(_d_globalCapacityf);
      freeGPU(_d_globalConductivityf);
      freeGPU(_d_capacity_map_cap);
      freeGPU(_d_capacity_map_cond);*/
      //cudaFree(_d_materials);
      cudaFree(_d_shapeFunctionTable);
      cudaFree(_d_thermal_boundaries);

      cudaFree(_d_materials_cap_ids);
      cudaFree(_d_materials_cond_ids);
      cudaFree(_d_localCapacityf);
      //free(_h_globalConductivity);
      cudaFree(_d_localConductivityf);
      cudaFree(_d_local_temperatures_cap_array);
      cudaFree(_d_local_jacobian_cap_array);
      cudaFree(_d_local_weight_cap_array);
      cudaFree(_d_local_temperatures_cond_array);
      cudaFree(_d_local_jacobian_cond_array);
      cudaFree(_d_local_weight_cond_array);
      clockFullStats(microGPU_conductivity, "GPU SOA calcConductivityMatrix");
      clockFullStats(microGPU_capacity,     "GPU SOA calcCapacityMatrix");
      clockFullStats(microGPU1,             "GPU AssembleCapacityMatrix");
      clockFullStats(microGPU2,             "GPU AssembleConductivityMatrix");
    }

#endif

}

/**
 * @brief Cascade initialization funtion. Calls the initialize methods for each of the Cells
 *        and tells them to compute their shapefunction values. Both loops are parallelized.
 *
 * @return void
 **/
void Body::initialize()
{
    //std::cout << "1. Body::initialize() :Initializing body tables "<< std::endl;
   _h_materials = (MaterialTable*) malloc(1*sizeof(MaterialTable));
   _h_shapeFunctionTable = (ShapeFunctionTable*) malloc(1*sizeof(ShapeFunctionTable));
   //_h_thermal_boundaries = (ThermalBoundaryTable*)malloc(1*sizeof(ThermalBoundaryTable));
    lastNode = nodes.back();
    auto end_int = this->cells.size();
    _number_cells = end_int;//
    //std::cout<< std::endl << "2. Body::initialize() Number of Cells " << _number_cells << std::endl ;
    _MC_points_per_cell = 3;//TODO MAGIC NUMBER HERE!!!!!!!
    _points_per_cell = 1;//TODO MAGIC NUMBER HERE!!!!!!!
    _number_points_MC = _number_cells * _MC_points_per_cell;
    _number_points = _number_cells * _points_per_cell;

    std::cout << "last node " << lastNode << std::endl;
    nodes.insert(nodes.end(), bondedBodyNodes.begin(), bondedBodyNodes.end());

    //std::cout << "3. Body::initialize() Initializing nodes " << std::endl;
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->initialize(this->nodes);
    }
   //std::cout << std::endl<< "SUPPORT SIZE MC = " << this->cells[0]->getSupportSizeMC() << std::endl;
   //std::cout << std::endl<< "SUPPORT SIZE = " << this->cells[0]->getSupportSize() << std::endl;
   _support_node_size = this->cells[0]->getSupportSize();//TODO MAGIC NUMBER HERE!!!!!!!
   _dim = 2; //TODO MAGIC NUMBER HERE!!!!!!!

  //std::cout << "4. Body::initialize() Counting MC Gausspoints" << std::endl;
        int nump_mc = 0;
        for (auto i = 0u; i < _number_cells; ++i) {
          nump_mc +=  this->cells[i]->getNumPoints_MC();
        }

    //std::cout << "5. Body::initialize() Counting Gausspoints" << std::endl;
        int nump = 0;
        for (auto i = 0u; i < _number_cells; ++i) {
            nump +=  this->cells[i]->getNumPoints();
        }
  //std::cout<< std::endl << "assumed MC points " << _number_points_MC << " vs counted " << nump_mc <<std::endl;
  //std::cout<< std::endl << "assumed points " << _number_points << " vs counted " << nump <<std::endl;

//std::cout << "6. Body::initialize() Computing shape shape functions" << std::endl;
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->computeShapeFunctions();
    }
  //setupShapeTables();
   //std::cout << "7. Body::initialize() Setting Material Ids" << std::endl;
   _h_materials_cap_ids = (int*) malloc( _number_points_MC * sizeof(int));
    for (int i = 0; i < _number_cells; ++i) {
        int cellMaterial = this->cells[i]->getMaterialId();
        for(int p = 0; p < _MC_points_per_cell; p++)
          _h_materials_cap_ids[i * _MC_points_per_cell + p] = cellMaterial;
    }

    _h_materials_cond_ids = (int*) malloc( _number_points * sizeof(int));
     for (int i = 0; i < _number_cells; ++i) {
         int cellMaterial = this->cells[i]->getMaterialId();
         for(int p = 0; p < _points_per_cell; p++)
           _h_materials_cond_ids[i * _points_per_cell + p] = cellMaterial;
     }

  if(GPU){
    CudaSafeCall(cudaMalloc((void**)&_d_materials_cap_ids,  _number_points_MC * sizeof(int)));
    CudaSafeCall(cudaMalloc((void**)&_d_materials_cond_ids, _number_points * sizeof(int)));

    CudaSafeCall(cudaMemcpy(_d_materials_cap_ids,  _h_materials_cap_ids,  _number_points_MC * sizeof(int), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(_d_materials_cond_ids, _h_materials_cond_ids, _number_points * sizeof(int),    cudaMemcpyHostToDevice));
  }

//std::cout << "8. Body::initialize() Iteration in nodes" << std::endl;
//The iteration on nodes MUST be done AFTER the cells.
    end_int = this->nodes.size();
    _number_nodes = end_int;
//     #pragma omp parallel for
    for (auto i = 0u; i < _number_nodes; ++i) {
        if (this->nodes[i]->getShapeFunType() == "RBF" ||
            this->nodes[i]->getShapeFunType() == "MLS") {
                this->nodes[i]->findSupportNodes(this->nodes);
        }
    }
//std::cout << "9. Body::initialize() Solving shape funcs" << std::endl;
//     #pragma omp parallel for
    for (auto i = 0u; i < _number_nodes; ++i) {
        if (this->nodes[i]->getShapeFunType() == "RBF") {
            this->nodes[i]->shapeFunSolve("RBF", 1.03);
        }
        if (this->nodes[i]->getShapeFunType() == "MLS") {
            this->nodes[i]->shapeFunSolve("MLS", 1.03);
        }
    }
//std::cout << "9b. Body::initialize() LoadThermalBoundary1D" << std::endl;
    std::map<std::string, BoundaryGroup *>::iterator it_boundaryGroups;
    std::map<int,LoadThermalBoundary1D*> thermal_boundaries_map;
    int ib = 0;
    for (auto& group : boundaryGroups) {
        group.second->initialize();
        //std::cout << "load map size " << group.second->getLoadThermal()<< std::endl;
        //thermal_boundaries_map.insert( {ib, group.second->getLoadThermal()});
        ib++;
    }

    setupShapeTables();

    _h_thermal_map_MC.resize(_number_points_MC * _support_node_size);
    _h_thermal_map.resize(_number_points * _support_node_size);
    _h_node_map_MC.resize(_number_points_MC * _support_node_size * _support_node_size);
    _h_node_map.resize(_number_points * _support_node_size * _support_node_size);

//std::cout<< " 11. Body::initialize() Mapping all nodes"<< std::endl;
    for (int i = 0; i < _number_cells; ++i) {
      //1D MAP
      this->cells[i]->mapThermalNodesMC(_h_thermal_map_MC.data(),
                                        _support_node_size,
                                        i);
      this->cells[i]->mapThermalNodes(_h_thermal_map.data(),
                                      _support_node_size,
                                      i);
      //2D MAP
      this->cells[i]->mapNodesMC(_h_node_map_MC.data(),
                                 _support_node_size,
                                 i,
                                 _number_nodes);

      this->cells[i]->mapNodes(_h_node_map.data(),
                               _support_node_size,
                               i,
                               _number_nodes);
    }
//std::cout<< " 12. Body::initialize() Allocating presence matrices"<< std::endl;
    _h_presence_matrix_cap = (int*) malloc(_number_nodes * _number_nodes * sizeof(int));
    _h_presence_matrix_cond = (int*) malloc(_number_nodes * _number_nodes * sizeof(int));
//std::cout<< " 13. Body::initialize() Mapping Sparse Matrices for direct assembly"<< std::endl;
    mapConectivityCapacityMatrix();
    map_global_matrix(_full_map_cap,
                      _vec_ind_cap,
                      _cvec_ptr_cap,
                      _h_presence_matrix_cap,
                      _number_nodes,
                      _number_nodes,
                      USECSC);

    mapConectivityConductivityMatrix();
    map_global_matrix(_full_map_cond,
                      _vec_ind_cond,
                      _cvec_ptr_cond,
                      _h_presence_matrix_cond,
                      _number_nodes,
                      _number_nodes,
                      USECSC);
    _sparse_matrix_size = _cvec_ptr_cap[_number_nodes];

//std::cout<< " 14. Body::initialize() Allocating local matrices arrays"<< std::endl;
  _h_globalCapacity.resize(_sparse_matrix_size);
  _h_localCapacityf      = (data_type*)malloc(_number_points_MC * _support_node_size * _support_node_size * sizeof(data_type));
  _h_globalConductivity.resize(_sparse_matrix_size);
  _h_localConductivityf  = (data_type*)malloc(_number_points * _support_node_size * _support_node_size * sizeof(data_type));

//std::cout<< " 14. Body::initialize() Allocating local nodes arrays"<< std::endl;
  _h_local_temperatures_cap_array   = (data_type*)malloc(_number_points_MC * _support_node_size * sizeof(data_type));
  _h_local_jacobian_cap_array       = (data_type*)malloc(_number_points_MC * sizeof(data_type));
  _h_local_weight_cap_array         = (data_type*)malloc(_number_points_MC * sizeof(data_type));
  _h_local_temperatures_cond_array   = (data_type*)malloc(_number_points * _support_node_size * sizeof(data_type));
  _h_local_jacobian_cond_array       = (data_type*)malloc(_number_points * sizeof(data_type));
  _h_local_weight_cond_array         = (data_type*)malloc(_number_points * sizeof(data_type));

for(int i = 0; i < _number_points * _support_node_size; i++)
      _h_local_temperatures_cond_array[i] = 1.0;
for(int i = 0; i < _number_points_MC * _support_node_size; i++)
      _h_local_temperatures_cap_array[i] = 1.0;

//std::cout<< " 14. Body::initialize() populating local nodes arrays"<< std::endl;
  for(int eachcell = 0; eachcell < this->cells.size(); eachcell++){
    for (int eachpoint = 0; eachpoint < _MC_points_per_cell; eachpoint++){
      int pindex = eachcell * _MC_points_per_cell + eachpoint;
      _h_local_jacobian_cap_array[pindex] = this->cells[eachcell]->getJacobianMC(eachpoint);
      _h_local_weight_cap_array[pindex] = this->cells[eachcell]->getWeightMC(eachpoint);
    }
  }

  for(int eachcell = 0; eachcell < this->cells.size(); eachcell++){
    for (int eachpoint = 0; eachpoint < _points_per_cell; eachpoint++){
      int pindex = eachcell * _points_per_cell + eachpoint;
      _h_local_jacobian_cond_array[pindex] = this->cells[eachcell]->getJacobianP(eachpoint);
      _h_local_weight_cond_array[pindex] = this->cells[eachcell]->getWeight(eachpoint);
    }
  }
    //map_vector_thermal_numbers(_locaThermalNumbers);
    //GPU part
  #ifdef HAVE_CUDA
    if(GPU){
       //_sparse_matrix_size = _cvec_ptr[_number_nodes];//convention
       std::cout << "sparse matrix has " << _sparse_matrix_size << " elements" << std::endl;
       CudaSafeCall(cudaMalloc((void**)&_d_globalCapacityf,     _sparse_matrix_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_globalConductivityf, _sparse_matrix_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_localCapacityf,      _number_points_MC * _support_node_size * _support_node_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_localConductivityf,  _number_points * _support_node_size * _support_node_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_full_map_cap,        _number_nodes * _number_nodes * sizeof(uint)));
       CudaSafeCall(cudaMalloc((void**)&_d_full_map_cond,       _number_nodes * _number_nodes * sizeof(uint)));
       CudaSafeCall(cudaMalloc((void**)&_d_vec_ind_cap,         _vec_ind_cap.size() * sizeof(uint)));
       CudaSafeCall(cudaMalloc((void**)&_d_vec_ind_cond,        _vec_ind_cond.size() * sizeof(uint)));
       CudaSafeCall(cudaMalloc((void**)&_d_cvec_ptr_cap,        _cvec_ptr_cap.size() * sizeof(uint)));
       CudaSafeCall(cudaMalloc((void**)&_d_cvec_ptr_cond,       _cvec_ptr_cond.size() * sizeof(uint)));
       CudaSafeCall(cudaMemcpy(_d_full_map_cap,   _full_map_cap.data(),   _number_nodes * _number_nodes * sizeof(uint), cudaMemcpyHostToDevice));
       CudaSafeCall(cudaMemcpy(_d_full_map_cond,  _full_map_cond.data(),  _number_nodes * _number_nodes * sizeof(uint), cudaMemcpyHostToDevice));
       CudaSafeCall(cudaMemcpy(_d_vec_ind_cap,    _vec_ind_cap.data(),    _vec_ind_cap.size() * sizeof(uint), cudaMemcpyHostToDevice));
       CudaSafeCall(cudaMemcpy(_d_vec_ind_cond,   _vec_ind_cond.data(),   _vec_ind_cond.size() * sizeof(uint), cudaMemcpyHostToDevice));
       CudaSafeCall(cudaMemcpy(_d_cvec_ptr_cap,   _cvec_ptr_cap.data(),   _cvec_ptr_cap.size() * sizeof(uint), cudaMemcpyHostToDevice));
       CudaSafeCall(cudaMemcpy(_d_cvec_ptr_cond,  _cvec_ptr_cond.data(),  _cvec_ptr_cond.size() * sizeof(uint), cudaMemcpyHostToDevice));

       CudaSafeCall(cudaMalloc((void**)&_d_local_temperatures_cap_array, _number_points_MC * _support_node_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_local_temperatures_cond_array, _number_points * _support_node_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_local_jacobian_cap_array, _number_points_MC * _support_node_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_local_jacobian_cond_array, _number_points * _support_node_size * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_local_weight_cap_array, _number_points_MC * sizeof(data_type)));
       CudaSafeCall(cudaMalloc((void**)&_d_local_weight_cond_array, _number_points * sizeof(data_type)));

       CudaCheckError();
    }
 #endif
//std::cout<< " 15. Body::initialize() Setting Up SOA boundary tables"<< std::endl;
 //setupThermalBoundaryTable(&_h_thermal_boundaries,
//                           thermal_boundaries_map);
}

void Body::setupPthreadsParameters()
{
/*  std::cout<< " MCPU. Body::initialize() Creating pthreads structures"<< std::endl;
  _param_array_capacity     = (p_struct*)malloc(MAX_THREADS * sizeof(p_struct));
  _param_array_conductivity = (p_struct*)malloc(MAX_THREADS * sizeof(p_struct));
  for(int i = 0; i < MAX_THREADS; i++){
      _param_array_capacity[i].globalMatrix = (std::atomic<double>*)_h_globalCapacity.data();
      _param_array_capacity[i].fullMap = &_full_map_cap;
      _param_array_capacity[i].local_matrices_array = _h_localCapacityf;
      _param_array_capacity[i].numCells = _number_points_MC;//this->cells.size();
      _param_array_capacity[i].supportNodeSize = _support_node_size;
      _param_array_capacity[i].thread_id = i;
      _param_array_capacity[i].max_threads = MAX_THREADS;
  }
  for(int i = 0; i < MAX_THREADS; i++){
      _param_array_conductivity[i].globalMatrix = (std::atomic<double>*)_h_globalConductivity.data();
      _param_array_conductivity[i].fullMap = &_full_map_cond;
      _param_array_conductivity[i].local_matrices_array = _h_localConductivityf;
      _param_array_conductivity[i].numCells = _number_points;//this->cells.size();
      _param_array_conductivity[i].supportNodeSize = _support_node_size;
      _param_array_conductivity[i].thread_id = i;
      _param_array_conductivity[i].max_threads = MAX_THREADS;
  }
  //struct p_calc_SOA_cap_struct
  _param_calc_capacity = (p_calc_SOA_cap_struct*)malloc(MAX_THREADS * sizeof(p_calc_SOA_cap_struct));
  for(int i = 0; i < MAX_THREADS; i++){
      _param_calc_capacity[i].local_capacity_array = _h_localCapacityf;
      _param_calc_capacity[i].local_temperatures_cap_array =_h_local_temperatures_cap_array ;
      _param_calc_capacity[i].local_weight_cap_array = _h_local_weight_cap_array;
      _param_calc_capacity[i].local_jacobian_cap_array = _h_local_jacobian_cap_array;
      _param_calc_capacity[i].local_shapes_phis_array = _h_local_shapes_cap_0.data();
      _param_calc_capacity[i].material_ids = _h_materials_cap_ids;
      _param_calc_capacity[i].materials = _h_materials;
      _param_calc_capacity[i].number_points = _number_points_MC;
      _param_calc_capacity[i].supportNodeSize = _support_node_size;
      _param_calc_capacity[i].thread_id = i;
      _param_calc_capacity[i].max_threads = MAX_THREADS;
    }
    //struct p_calc_SOA_cond_struct{
    _param_calc_conductivity = (p_calc_SOA_cond_struct*)malloc(MAX_THREADS * sizeof(p_calc_SOA_cond_struct));
    for(int i = 0; i < MAX_THREADS; i++){
        _param_calc_conductivity[i].local_conductivity_array = _h_localConductivityf;
        _param_calc_conductivity[i].local_temperatures_cond_array = _h_local_temperatures_cond_array;
        _param_calc_conductivity[i].local_weight_cond_array = _h_local_weight_cond_array;
        _param_calc_conductivity[i].local_jacobian_cond_array = _h_local_jacobian_cond_array;
        _param_calc_conductivity[i].local_shapes_phis_array = _h_local_shapes_cond_0.data();
        _param_calc_conductivity[i].local_shapes_phis_dim_array = _h_local_shapes_cond_dim.data();
        _param_calc_conductivity[i].material_ids = _h_materials_cond_ids;
        _param_calc_conductivity[i].materials = _h_materials;
        //std::cout << "Materials pointer = " << _h_materials << std::endl;
        //std::cout << "Materials pointer = " << _param_calc_conductivity[i].materials << std::endl;
        //std::cout << "Materials counter = " << _h_materials->number_materials << std::endl;
        //std::cout << "Materials counter = " << _param_calc_conductivity[i].materials->number_materials << std::endl;
        _param_calc_conductivity[i].number_points = _number_points;
        _param_calc_conductivity[i].supportNodeSize = _support_node_size;
        _param_calc_conductivity[i].thread_id = i;
        _param_calc_conductivity[i].max_threads = MAX_THREADS;
    }*/
}


void Body::setupShapeTables()
{
  /////////////////////////////////////////////////////////////////////////////
  //std::cout << "Resizing ShapeTables" << std::endl;
  _h_local_shapes_cap_0.resize(_number_points_MC * _support_node_size);
  ///////////CAPACITY PART ////////////////
//  std::cout << "Copying ShapeTable 0: " << std::endl;
      for(int eachcell = 0; eachcell < _number_cells; eachcell++){
          for (int lp = 0; lp < _MC_points_per_cell; lp++){
             for(int lNode = 0; lNode < _support_node_size ; lNode++){
               int shapeIndex = (eachcell * _MC_points_per_cell * _support_node_size) + (lp * _support_node_size) + lNode;
               _h_local_shapes_cap_0[shapeIndex] = this->cells[eachcell]->getCapPhi(lp,0,lNode);
          }
        }
      }

  _h_local_shapes_cond_0.resize(_number_points * _support_node_size);
  _h_local_shapes_cond_dim.resize(_number_points * _support_node_size * _support_node_size);
    ///////////CONDUCTIVITY PART ////////////////
      //std::cout << "Copying ShapeTable 0: " << std::endl;
          for(int eachcell = 0; eachcell < _number_cells; eachcell++){
              for (int eachpoint = 0; eachpoint < _points_per_cell; eachpoint++){
                 for(int lNode = 0; lNode < _support_node_size ; lNode++){
                   int shapeIndex = eachcell * _points_per_cell * _support_node_size + eachpoint * _support_node_size + lNode;
                   _h_local_shapes_cond_0[shapeIndex] = this->cells[eachcell]->getCondPhi(eachpoint,0,lNode);
              }
            }
          }
  //std::cout << "Copying ShapeTable Dim" << std::endl;
     int spns2 = _support_node_size * _support_node_size;
      for(int eachcell = 0; eachcell < _number_cells; eachcell++){
        for (int eachpoint = 0; eachpoint < _points_per_cell; eachpoint++){
            for(int rowNode = 0; rowNode < _support_node_size ; rowNode++){
              for(int colNode = 0; colNode < _support_node_size ; colNode++){
                  int shapeIndex = (eachcell * _points_per_cell * spns2) + (eachpoint * spns2) + rowNode * _support_node_size + colNode;
                  for(int idim = 1; idim < 3 ; idim++){
                    _h_local_shapes_cond_dim[shapeIndex] += this->cells[eachcell]->getCondPhi(eachpoint,idim,rowNode) * this->cells[eachcell]->getCondPhi(eachpoint,idim,colNode);
                  }
              }
            }
      }
      }

  if(GPU){
    CudaSafeCall(cudaMalloc((void**) &_d_local_shapes_cap_0,    _number_points_MC * _support_node_size * sizeof(double)));
    CudaSafeCall(cudaMalloc((void**) &_d_local_shapes_cond_0,   _number_points * _support_node_size * sizeof(double)));
    CudaSafeCall(cudaMalloc((void**) &_d_local_shapes_cond_dim, _number_points * _support_node_size * _support_node_size * sizeof(double)));

    CudaSafeCall(cudaMemcpy(_d_local_shapes_cap_0,    _h_local_shapes_cap_0.data(),    _number_points_MC * _support_node_size * sizeof(double), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(_d_local_shapes_cond_0,   _h_local_shapes_cond_0.data(),   _number_points * _support_node_size * sizeof(double), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(_d_local_shapes_cond_dim, _h_local_shapes_cond_dim.data(), _number_points * _support_node_size * _support_node_size * sizeof(double), cudaMemcpyHostToDevice));
  }
}

void Body::setThermalBoundaryTable(ThermalBoundaryTable *tb_ptr)
{
  _h_thermal_boundaries = tb_ptr;
}

void Body::setMaterialTable(MaterialTable* mt_ptr)
{
  _h_materials = mt_ptr;
  _d_materials = (MaterialTable*)malloc(sizeof(MaterialTable));
  if(GPU){
    allocateTransferMaterialTable( &_d_materials, _h_materials);
  }
 //debug_printMaterialTable(_h_materials);
}

void Body::setTemperatureVector(lmx::Vector<data_type>& q)
{
  for(int ip = 0 ; ip < _number_points_MC ; ip++){
    for(int lnode = 0 ; lnode < _support_node_size; lnode++){
      int myIndex = ip *_support_node_size + lnode;
      int myNode = _h_thermal_map_MC[myIndex];
      _h_local_temperatures_cap_array[myIndex] = q.readElement(myNode);
    }
  }
  for(int ip = 0 ; ip < _number_points ; ip++){
    for(int lnode = 0 ; lnode < _support_node_size; lnode++){
      int myIndex = ip *_support_node_size + lnode;
      int myNode = _h_thermal_map[myIndex];
      _h_local_temperatures_cond_array[myIndex] = q.readElement(myNode);
    }
  }
}

void Body::setQVector(const lmx::Vector<data_type>& q)
{
  for(int ip = 0 ; ip < _number_points_MC ; ip++){
    for(int lnode = 0 ; lnode < _support_node_size; lnode++){
      int myIndex = ip *_support_node_size + lnode;
      int myNode = _h_thermal_map_MC[myIndex];
      _h_local_temperatures_cap_array[myIndex] = q.readElement(myNode);
    }
  }
  for(int ip = 0 ; ip < _number_points ; ip++){
    for(int lnode = 0 ; lnode < _support_node_size; lnode++){
      int myIndex = ip *_support_node_size + lnode;
      int myNode = _h_thermal_map[myIndex];
      _h_local_temperatures_cond_array[myIndex] = q.readElement(myNode);
    }
  }
}

/**
 * @brief Computes the local Capacity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::calcCapacityMatrix()
{
  if(OLD_CODE)
  {
    cpuClock cck1;
    cpuTick(&cck1);
      auto end_int = this->cells.size();
  //#pragma omp parallel for
      for (auto i = 0u; i < end_int; ++i) {
          this->cells[i]->computeCapacityGaussPoints();
      }
      cpuTock(&cck1, "Existing CPU calcCapacityMatrix ");
      microCPU_old_capacity.push_back(cck1.elapsedMicroseconds);

  }
  else if(NEWCPU)
  {

    cpuClock cck;
    cpuTick(&cck);
    computeSOACapacityMatrix(_h_localCapacityf,
                             _h_local_temperatures_cap_array,
                             _h_local_weight_cap_array,
                             _h_local_jacobian_cap_array,
                             _h_local_shapes_cap_0.data(),
                             _h_materials_cap_ids,
                             _h_materials,
                             _number_points_MC,
                             _support_node_size,
                             0,
                             1);
    cpuTock(&cck, "New Single CPU calcCapacityMatrix ");
    microCPU_single_capacity.push_back(cck.elapsedMicroseconds);

  }
  //
  else if(MULTICPU)
  {
    //pthread_t _threads[MAX_THREADS];
      omp_set_num_threads(MAX_THREADS);
    //std::cout << "About to run New Multi CPU calcConductivityMatrix " << std::endl;
    cpuClock cck;
    cpuTick(&cck);

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      computeSOACapacityMatrix(_h_localCapacityf,
                               _h_local_temperatures_cap_array,
                               _h_local_weight_cap_array,
                               _h_local_jacobian_cap_array,
                               _h_local_shapes_cap_0.data(),
                               _h_materials_cap_ids,
                               _h_materials,
                               _number_points_MC,
                               _support_node_size,
                               tid,
                               MAX_THREADS);
    }

    /*for(int i = 0; i < MAX_THREADS; i++)
        pthread_create(&_threads[i],NULL,computeCapacityThreadWrapper,(void*)&(_param_calc_capacity[i]));
    for(int i = 0; i < MAX_THREADS; i++)
        pthread_join(_threads[i],NULL);*/

    cpuTock(&cck, "New Multi CPU calcCapacityMatrix ");
    microCPU_multi_capacity.push_back(cck.elapsedMicroseconds);

  }
  //
  if(GPU)
  {
    cudaClock gck1;
    cudaTick(&gck1);
    gpu_computeSOACapacityMatrix(_d_localCapacityf,
                                 _d_local_temperatures_cap_array,
                                 _d_local_weight_cap_array,
                                 _d_local_jacobian_cap_array,
                                 _d_local_shapes_cap_0,
                                 _d_materials_cap_ids,
                                 _d_materials,
                                 _number_points_MC,
                                 _support_node_size,
                                 128,
                                 0);
    cudaTock(&gck1, "GPU calcCapacityMatrix");
    microGPU_capacity.push_back(gck1.elapsedMicroseconds);
    CudaCheckError();
  }

}

/**
 * @brief Computes the local Conductivity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::calcConductivityMatrix()
{
   if(OLD_CODE)
   {
    auto end_int = this->cells.size();
     cpuClock cck1;
     cpuTick(&cck1);
      for (auto i = 0u; i < end_int; ++i) {
          this->cells[i]->computeConductivityGaussPoints();
      }
      cpuTock(&cck1, "\nExisting CPU calcConductivityMatrix ");
      microCPU_old_conductivity.push_back(cck1.elapsedMicroseconds);

  }
  else if(NEWCPU)
  {
    //std::cout << " Body::calcConductivityMatrix() " << std::endl;
    int number_points = this->cells.size();
    cpuClock cck;
    cpuTick(&cck);

    computeSOAConductivityMatrix(_h_localConductivityf,
                                 _h_local_temperatures_cond_array,
                                 _h_local_weight_cond_array,
                                 _h_local_jacobian_cond_array,
                                 _h_local_shapes_cond_0.data(),
                                 _h_local_shapes_cond_dim.data(),
                                 _h_materials_cond_ids,
                                 _h_materials,
                                 _number_points,
                                 _support_node_size,
                                 0,
                                 1);

    cpuTock(&cck, "\nNew Single CPU calcConductivityMatrix ");
    microCPU_single_conductivity.push_back(cck.elapsedMicroseconds);

  }
  //
  else if(MULTICPU)
  {
    //pthread_t _threads[MAX_THREADS];
    //std::cout << "About to run New Multi CPU calcConductivityMatrix " << std::endl;
      omp_set_num_threads(MAX_THREADS);
    cpuClock cck;
    cpuTick(&cck);
    /*for(int i = 0; i < MAX_THREADS; i++)
        pthread_create(&_threads[i],NULL,computeConductivityThreadWrapper,(void*)&(_param_calc_conductivity[i]));
   for(int i = 0; i < MAX_THREADS; i++)
       pthread_join(_threads[i],NULL);*/
  #pragma omp parallel
  {
       int tid = omp_get_thread_num();
       //std::cout << "tid = " << tid << std::endl;
       computeSOAConductivityMatrix(_h_localConductivityf,
                                    _h_local_temperatures_cond_array,
                                    _h_local_weight_cond_array,
                                    _h_local_jacobian_cond_array,
                                    _h_local_shapes_cond_0.data(),
                                    _h_local_shapes_cond_dim.data(),
                                    _h_materials_cond_ids,
                                    _h_materials,
                                    _number_points,
                                    _support_node_size,
                                    tid,
                                    MAX_THREADS);
    }
    cpuTock(&cck, "\nNew Multi CPU calcConductivityMatrix ");
    microCPU_multi_conductivity.push_back(cck.elapsedMicroseconds);
  }
  //
  if(GPU)
  {
    CudaCheckError();
    cudaClock gck1;
    cudaTick(&gck1);
    gpu_computeSOAConductivityMatrix(_d_localConductivityf,
                                     _d_local_temperatures_cond_array,
                                     _d_local_weight_cond_array,
                                     _d_local_jacobian_cond_array,
                                     _d_local_shapes_cond_0,
                                     _d_local_shapes_cond_dim,
                                     _d_materials_cond_ids,
                                     _d_materials,
                                     _number_points,
                                     _support_node_size,
                                     128,
                                     0);
    cudaTock(&gck1, "GPU SOA calcConductivityMatrix");
    microGPU_conductivity.push_back(gck1.elapsedMicroseconds);
    CudaCheckError();
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
        this->cells[i]->presenceCapacityGaussPoints(_h_presence_matrix_cap, _number_nodes);
    }
}

/**
 * @brief Computes the local Capacity of the material body by calling each cell's cascade function.
 *
 * @return void
 **/
void Body::mapConectivityConductivityMatrix()
{
    auto end_int = this->cells.size();
    std::cout << "Mapping conductivity matrix with "<< end_int << " cells" <<std::endl;
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->presenceConductivityGaussPoints(_h_presence_matrix_cond, _number_nodes);
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
   //std::cout << " Inside  Body::calcExternalHeat()"<< std::endl;
    if (loadThermalBody) {
       //std::cout << "loadThermalBody!!!"<< std::endl;//not used
        for (auto i = 0u; i < end_int; ++i) {
            this->cells[i]->computeQextGaussPoints(this->loadThermalBody);
        }
    }
    for (auto group : boundaryGroups) {
      //std::cout << "AHA!!! for (auto group : boundaryGroups)"<< std::endl;//this is the one
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
  if(OLD_CODE){
    cpuClock cck1;
    cpuTick(&cck1);
    auto end_int = this->cells.size();
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleCapacityGaussPoints(globalCapacity);
    }
  cpuTock(&cck1, "CPU assembleCapacityMatrix");
  microCPU1.push_back(cck1.elapsedMicroseconds);

}else if(NEWCPU){

    cpuClock cck1b;
    cpuTick(&cck1b);

    init_host_array_to_value(_h_globalCapacity.data(), 0.0, _sparse_matrix_size);//now reset in simulat

    AssembleGlobalMatrix(_h_globalCapacity,
                         _full_map_cap,
                         _h_node_map_MC,
                         _h_localCapacityf,
                         _number_points_MC,
                         _support_node_size,
                         USECSC);
    cpuTock(&cck1b, "CPU assembleCapacityMatrixWithMap");
    microCPU1b.push_back(cck1b.elapsedMicroseconds);

    cast_into_lmx_type(globalCapacity,
                       _h_globalCapacity,
                       _vec_ind_cap,
                       _cvec_ptr_cap,
                       _number_nodes,
                       _number_nodes,
                       USECSC);

}
else if(MULTICPU){
  ///////multi-cpu part
  //pthread_t _threads[MAX_THREADS];
    omp_set_num_threads(MAX_THREADS);
    cpuClock cck1c;
    cpuTick(&cck1c);
   init_host_array_to_value(_h_globalCapacity, 0.0, _sparse_matrix_size);
  /* for(int i = 0; i < MAX_THREADS; i++)
       pthread_create(&_threads[i],NULL,threadWrapper,(void*)&(_param_array_capacity[i]));
   for(int i = 0; i < MAX_THREADS; i++)
       pthread_join(_threads[i],NULL);*/
     #pragma omp parallel
      {
        int tid = omp_get_thread_num();
        atomicAssembleGlobalMatrix((std::atomic<double>*)_h_globalCapacity.data(),
                             _full_map_cap,
                             _h_node_map_MC,
                             _h_localCapacityf,
                             _number_points_MC,
                             _support_node_size,
                             tid,
                             MAX_THREADS,
                             USECSC);

      }

   cpuTock(&cck1c, "MULTI CPU assembleCapacityGaussPointsWithMap");
   microCPU1c.push_back(cck1c.elapsedMicroseconds);

   cast_into_lmx_type(globalCapacity,
                      _h_globalCapacity,
                      _vec_ind_cap,
                      _cvec_ptr_cap,
                      _number_nodes,
                      _number_nodes,
                      USECSC);

  }
  if(GPU){
  #ifdef HAVE_CUDA
  cudaClock gck1;

    cudaTick(&gck1);
     init_array_to_value(_d_globalCapacityf, 0.0, _sparse_matrix_size,128);
     gpu_assemble_global_matrix(_d_globalCapacityf,
                                  _d_capacity_map_cap,
                                  _d_localCapacityf,
                                  this->cells.size(),
                                  _support_node_size,
                                  _number_points,
                                  128);
     cudaTock(&gck1, "GPU assembleCapacityMatrix");
     microGPU1.push_back(gck1.elapsedMicroseconds);
     CudaCheckError();
   #endif
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
if(OLD_CODE) {
    cpuClock cck2;
    cpuTick(&cck2);
    auto end_int = this->cells.size();
//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i) {
        this->cells[i]->assembleConductivityGaussPoints(globalConductivity);
    }
    cpuTock(&cck2, "CPU assembleConductivityMatrix");
    microCPU2.push_back(cck2.elapsedMicroseconds);
    /*std::cout << " SumSum of globalConductivity = " << globalConductivity.sumSum() << std::endl;
    std::cout << " Trace of globalConductivity = " << globalConductivity.trace() << std::endl;
    std::cout << " StdDev of globalConductivity = " << globalConductivity.stddev() << std::endl;
    std::cout << " Max of globalConductivity = " << globalConductivity.max() << std::endl;
    std::cout << " Min of globalConductivity = " << globalConductivity.min() << std::endl;*/
    //std::cout << " NonZeroes of globalConductivity = " << globalConductivity.numNonZeros() << std::endl;
    //new CPU assembly function
} else if(NEWCPU){
  //std::cout << "inside assembleConductivityMatrix" << std::endl;
    cpuClock cck2b;
    cpuTick(&cck2b);
    //std::cout << "before init_host_array_to_value" << std::endl;
    //auto end_int = this->cells.size();
    init_host_array_to_value(_h_globalConductivity, 0.0, _sparse_matrix_size);
    //std::cout << "before assembleConductivityGaussPointsWithMap" << std::endl;
    AssembleGlobalMatrix(_h_globalConductivity,
                         _full_map_cond,
                         _h_node_map,
                         _h_localConductivityf,
                         _number_points,
                         _support_node_size,
                         USECSC);
    cpuTock(&cck2b, "CPU assembleConductivityGaussPointsWithMap");
    /*double debugio =0.0;
    for (auto& el : _h_globalConductivity){
        debugio += el;}
  std::cout << "_h_globalConductivity sumsum = " << debugio << std::endl;*/

    microCPU2b.push_back(cck2b.elapsedMicroseconds);
    cpuClock cck2b1;
    cpuTick(&cck2b1);
    cast_into_lmx_type(globalConductivity,
                       _h_globalConductivity,
                       _vec_ind_cond,
                       _cvec_ptr_cond,
                       _number_nodes,
                       _number_nodes,
                       USECSC);

    cpuTock(&cck2b1, "CPU cast_into_lmx_type");
    /*std::cout << " after the cast, leaving assembleConductivityMatrix " << std::endl;
    std::cout << " SumSum of globalConductivity = " << globalConductivity.sumSum() << std::endl;
    std::cout << " Trace of globalConductivity = " << globalConductivity.trace() << std::endl;
    std::cout << " StdDev of globalConductivity = " << globalConductivity.stddev() << std::endl;
    std::cout << " Max of globalConductivity = " << globalConductivity.max() << std::endl;
    std::cout << " Min of globalConductivity = " << globalConductivity.min() << std::endl;*/
    //std::cout << " NonZeroes of globalConductivity = " << globalConductivity.numNonZeros() << std::endl;
  } else if(MULTICPU){
    //pthread_t _threads[MAX_THREADS];
      omp_set_num_threads(MAX_THREADS);
      cpuClock cck2c;
      cpuTick(&cck2c);
      init_host_array_to_value(_h_globalConductivity, 0.0, _sparse_matrix_size);
      #pragma omp parallel
       {
         int tid = omp_get_thread_num();
         atomicAssembleGlobalMatrix((std::atomic<double>*)_h_globalConductivity.data(),
                                    _full_map_cond,
                                    _h_node_map,
                                    _h_localConductivityf,
                                    _number_points,
                                    _support_node_size,
                                    tid,
                                    MAX_THREADS,
                                    USECSC);
       }


     /*for(int i = 0; i < MAX_THREADS; i++)
         pthread_create(&_threads[i],NULL,threadWrapper,(void*)&(_param_array_conductivity[i]));
     for(int i = 0; i < MAX_THREADS; i++)
         pthread_join(_threads[i],NULL);*/
     cpuTock(&cck2c, "MULTI CPU assembleConductivityGaussPointsWithMap");
     microCPU2c.push_back(cck2c.elapsedMicroseconds);

     cast_into_lmx_type(globalConductivity,
                        _h_globalConductivity,
                        _vec_ind_cond,
                        _cvec_ptr_cond,
                        _number_nodes,
                        _number_nodes,
                        USECSC);
}
////////
  #ifdef HAVE_CUDA
    cudaClock gck2;
    if(GPU){
      init_array_to_value(_d_globalConductivityf, 0.0, _sparse_matrix_size,128);
      cudaTick(&gck2);
      gpu_assemble_global_matrix(_d_globalConductivityf,
                                _d_capacity_map_cond,
                                _d_localConductivityf,
                                this->cells.size(),
                                _support_node_size,
                                _number_points,
                                128);
    cudaTock(&gck2, "GPU assembleConductivityMatrix");
    microGPU2.push_back(gck2.elapsedMicroseconds);
    CudaCheckError();
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
  if(OLD_CODE){

      auto end_int = this->cells.size();
     for (auto i = 0u; i < end_int; ++i) {
        //std::cout <<"Inside Body::assembleExternalHeat  for (auto i = 0u; i < end_int; ++i)" <<std::endl;
          this->cells[i]->assembleQextGaussPoints(globalExternalHeat);
      }
      //std::cout << globalExternalHeat << std::endl;
      for (auto group : boundaryGroups) {
        //std::cout <<"Inside Body::assembleExternalHeat for (auto group : boundaryGroups)" <<std::endl;
          group.second->assembleExternalHeat(globalExternalHeat);
      }
      //std::cout << globalExternalHeat << std::endl;
  } else if(NEWCPU) {
    //std::cout <<"Inside Body::assembleExternalHeat" <<std::endl;
      auto end_int = this->cells.size();
      //std::cout <<"Body::assembleExternalHeat Initializing global vector" <<std::endl;
      //globalExternalHeat.fillIdentity(0.0f);
      //std::cout <<"Body::assembleExternalHeat from Cells" <<std::endl;
      for (auto i = 0u; i < end_int; ++i) {
          this->cells[i]->assembleQextGaussPoints(globalExternalHeat);

      }
    //  std::cout <<"Body::assembleExternalHeat from Boundary Groups" <<std::endl;
      for (auto group : boundaryGroups) {
          group.second->assembleExternalHeat(globalExternalHeat);
      }
  } else if(MULTICPU){
    //std::cout <<"Inside Body::assembleExternalHeat" <<std::endl;
      auto end_int = this->cells.size();
      //std::cout <<"Body::assembleExternalHeat Initializing global vector" <<std::endl;
      //globalExternalHeat.fillIdentity(0.0f);
      //std::cout <<"Body::assembleExternalHeat from Cells" <<std::endl;
      for (auto i = 0u; i < end_int; ++i) {
          this->cells[i]->assembleQextGaussPoints(globalExternalHeat);

      }
    //  std::cout <<"Body::assembleExternalHeat from Boundary Groups" <<std::endl;
      for (auto group : boundaryGroups) {
          group.second->assembleExternalHeat(globalExternalHeat);
      }
  }else if(GPU){
  #ifdef HAVE_CUDA
    if(GPU){
      //init_array_to_value(_d_globalExternalHeat, 0.0, _number_nodes, 128);
      /*gpu_assemble_global_vector(_d_globalExternalHeat,
                                _d_locaThermalNumbers,
                                _d_localHeatf,
                                _support_node_size,
                                _number_points,
                                128);*/

    }
  #endif
 }
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
      //std::cout << "\n\nOUTPUT STEP \n\n"<< std::endl;
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
