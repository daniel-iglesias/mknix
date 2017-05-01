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

#ifndef MKNIXBODY_H
#define MKNIXBODY_H

// #include <string>
#include "common.h"
#include "LMX/lmx.h"
#include "gmm/gmm_matrix.h"
#include "boundarygroup.h"

#include <core/cellboundary.h>
#include <core/node.h>
#include <gpu/assembly_cpu.h>
#include "gpu/functions_cpu.h"
#ifdef HAVE_CUDA
//#include <gpu/cuda_helper.h>
#include <gpu/assembly_kernels.h>

#endif
//MULTICPU
#define MAXTHREADS 4

namespace mknix {

class Point;

// class Node;
class Cell;

// class BoundaryGroup;
class LoadThermalBody;

class LoadThermalBoundary1D;

/**
	@author AUTHORS <MAILS>
*/
class Body
{

public:
    Body();

    Body(std::string);

    virtual ~Body();

    virtual std::string getType() = 0;

    virtual void initialize();

    virtual void setMaterialTable(MaterialTable* mt_ptr);

    virtual void setTemperatureVector(lmx::Vector<data_type>& q);

    virtual void calcFactors();

    virtual void calcCapacityMatrix();

    virtual void calcConductivityMatrix();

    virtual void calcExternalHeat();

    virtual void assembleCapacityMatrix(lmx::Matrix<data_type>&);
    virtual void mapConectivityCapacityMatrix();

    virtual void assembleConductivityMatrix(lmx::Matrix<data_type>&);

    virtual void assembleExternalHeat(lmx::Vector<data_type>&);

    virtual void calcMassMatrix() = 0;

    virtual void calcExternalForces() = 0;

    virtual void assembleMassMatrix(lmx::Matrix<data_type>&) = 0;

    virtual void assembleExternalForces(lmx::Vector<data_type>&) = 0;

    void setTemperature(double);

    virtual void setMechanical()
    {
        this->isThermal = 0;
    }

    virtual void setOutput(std::string) = 0;

    virtual void outputStep
            (const lmx::Vector<data_type>&, const lmx::Vector<data_type>&) = 0;

    virtual void outputStep(const lmx::Vector<data_type>&) = 0;

    virtual void outputStep();

    virtual void outputToFile(std::ofstream *);

    virtual void addNode(Node * node_in)
    {
        this->nodes.push_back(node_in);
    }

    void addNodes(std::vector<Node *>& nodes_in)
    {
        this->bondedBodyNodes.insert(bondedBodyNodes.end(), nodes_in.begin(), nodes_in.end());
    }

    virtual int getNodesSize()
    {
        return this->nodes.size();
    }

    std::vector<Node *>& getNodes()
    {
        return this->nodes;
    }

    virtual Node * getNode(int node_number)
    {
        return this->nodes[node_number];
    }

    virtual Node * getLastNode()
    {
        return lastNode;
    }

    virtual void addBoundaryLine(Point * node1, Point * node2)
    {
        this->linearBoundary[node1] = node2;
//       cout << "linearBoundary: " << node1->getNumber()
//      << " " << node2->getNumber() << endl;
    }

    virtual void addBoundaryConnectivity(std::vector<int> connectivity_in);

//     void addBoundaryGroup( std::string boundary_name )
//     {
//       this->boundaryGroups[boundary_name];
//     }
//
//     BoundaryGroup& getBoundaryGroup( std::string boundary_name )
//     {
//       return this->boundaryGroups[boundary_name];
//     }

    virtual int getBoundarySize()
    {
        return this->linearBoundary.size();
    }

    virtual Point * getBoundaryFirstNode()
    {
        return this->linearBoundary.begin()->first;
    }

    virtual Point * getBoundaryNextNode(Point * node_in)
    {
        return this->linearBoundary[node_in];
    }

    void addBoundaryGroup(std::string boundaryName_in)
    {
        this->boundaryGroups[boundaryName_in] = new BoundaryGroup;
    }

    void addNodeToBoundaryGroup(int nodeNumber_in, std::string boundaryName_in)
    {
        Node * temp = nodes[nodeNumber_in];
        this->boundaryGroups[boundaryName_in]->addNode(temp);
    }

    void addCellToBoundaryGroup(CellBoundary * cell_in, std::string boundaryName_in)
    {
        this->boundaryGroups[boundaryName_in]->addCell(cell_in);
    }

    void setLoadThermalInBoundaryGroup(LoadThermalBoundary1D * load_in, std::string boundaryName_in)
    {
        this->boundaryGroups[boundaryName_in]->setLoadThermal(load_in);
    }

    std::map<int, Cell *>& getCells()
    {
        return this->cells;
    }

    int getCellLastNumber()
    {
        return this->cells.end()->first;
    }

    virtual void addCell(int num, Cell * cell_in)
    {
        this->cells[num] = cell_in;
    }

    virtual void writeBodyInfo(std::ofstream *) = 0;

    virtual void writeBoundaryNodes(std::vector<Point *>&) = 0;

    virtual void writeBoundaryConnectivity(std::vector<std::vector<Point *> >&) = 0;

    virtual void translate(double, double, double);

    // Temporary, should be a pointer to a load class
    virtual void setLoadThermal(LoadThermalBody * theLoad)
    {
        loadThermalBody = theLoad;
    }

    virtual MaterialTable* getMaterialTablePointer(){
      return this->_h_materials;
    }


protected:
    std::string title;
    Node * lastNode;
    std::vector<Node *> nodes;
    std::vector<Node *> bondedBodyNodes;
    std::map<std::string, BoundaryGroup *> boundaryGroups;
    std::vector<std::vector<int> > boundaryConnectivity;
    std::map<Point *, Point *> linearBoundary;
    /**< Map of linear boundary. */
    std::map<int, Cell *> cells;
    /**< Map of integration cells. */
    bool computeEnergy;
    bool isThermal;
    std::vector<lmx::Vector<data_type> *> temperature;
    LoadThermalBody * loadThermalBody;

    //map
    std::vector<uint> _full_map_cap;//also used for cpu new assembly functions
    std::vector<uint> _vec_ind_cap;
    std::vector<uint> _cvec_ptr_cap;
    std::vector<uint> _full_map_cond;
    std::vector<uint> _vec_ind_cond;
    std::vector<uint> _cvec_ptr_cond;
    std::vector<int> _locaThermalNumbers;

    //GPU related
    bool _use_gpu;
    data_type *_d_globalCapacity;
    data_type *_d_globalConductivity;
    data_type *_h_globalCapacity;
    data_type *_h_globalCapacity2;
    data_type *_h_globalConductivity;
    data_type *_h_globalConductivity2;
    data_type *_d_globalCapacityf;
    data_type *_d_globalConductivityf;
    data_type *_d_localCapacityf;
    data_type *_h_localCapacityf;
    data_type *_d_localConductivityf;
    data_type *_h_localConductivityf;

    data_type *_h_local_capacity_factor;
    data_type *_h_local_conductivity_factor;
    data_type *_h_local_temperatures_array;
    data_type *_h_local_jacobian_array;
    data_type *_h_local_weight_array;

    data_type *_d_local_capacity_factor;
    data_type *_d_local_conductivity_factor;
    data_type *_d_local_temperatures_array;
    data_type *_d_local_shapeFun_phis;
    data_type *_d_local_jacobian_array;
    data_type *_d_local_weight_array;

    MaterialTable *_h_materials;
    int *_h_materials_ids;
    ShapeFunctionTable *_h_shapeFunctionTable;

    int       *_d_capacity_map;
    int       *_h_presence_matrix;
    int       _number_nodes;
    int       _number_points;
    int       _support_node_size;
    int       _sparse_matrix_size;
    //MULTI CPU part
    pthread_t _threads[MAXTHREADS];
    p_struct* _param_array_capacity;
    p_struct* _param_array_conductivity;
    //measurement related
    std::vector<double> microCPU1, microCPU1b, microCPU1c,
                        microCPU2, microCPU2b, microCPU2c, microCPU3;
    std::vector<double> microGPU1, microGPU2, microGPU3;
    std::vector<double> microCPU_old_capacity, microCPU_single_capacity, microCPU_multi_capacity;
    std::vector<double> microCPU_old_conductivity, microCPU_single_conductivity, microCPU_multi_conductivity;
    std::vector<double> microCPU_single_factors, microCPU_multi_factors;
    std::vector<double> microGPU_capacity, microGPU_conductivity, microGPU_factors;
    //testing gmm casting
    gmm::csr_matrix<data_type> gmm_sp_globalCapacity;

};

}

#endif
