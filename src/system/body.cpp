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

#include "body.h"

#include <core/cell.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <map>
#include <sstream>
#include <iomanip>

#ifdef HAVE_VTK
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkCellType.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#endif 

namespace mknix
{

/**
 * @brief Default constructor. Initializes the body with thermal simulation enabled and energy computation disabled.
 */
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

/**
 * @brief Destructor. Releases memory for temperature vectors, cells, boundary groups and volumetric heat sources.
 */
Body::~Body()
{
    for (auto& temp : temperature)
    {
        delete temp;
    }
    for (auto& cell : cells)
    {
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
    for (auto& group : boundaryGroups)
    {
        delete group.second;
    }
    for (auto p : volumetricHeatSources)
    {
        delete p;
    } 
    volumetricHeatSources.clear();
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
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->initialize(this->nodes);
    }

//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i)
    {
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
    for (auto i = 0u; i < end_int; ++i)
    {
        if (this->nodes[i]->getShapeFunType() == "RBF" ||
                this->nodes[i]->getShapeFunType() == "MLS")
        {
            this->nodes[i]->findSupportNodes(this->nodes);
        }
    }

//     #pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i)
    {
        if (this->nodes[i]->getShapeFunType() == "RBF")
        {
            this->nodes[i]->shapeFunSolve("RBF", 1.03);
        }
        if (this->nodes[i]->getShapeFunType() == "MLS")
        {
            this->nodes[i]->shapeFunSolve("MLS", 1.03);
        }
    }
    std::map<std::string, BoundaryGroup *>::iterator it_boundaryGroups;
    for (auto& group : boundaryGroups)
    {
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
//#pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i)
    {
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
//#pragma omp parallel for
    for (auto i = 0u; i < end_int; ++i)
    {
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
    for (auto i = 0u; i < end_int; ++i){
            this->cells[i]->computeQextGaussPoints(this->volumetricHeatSources);
        }
    for (auto group : boundaryGroups){
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
    for (auto i = 0u; i < end_int; ++i)
    {
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
    for (auto i = 0u; i < end_int; ++i)
    {
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
    for (auto i = 0u; i < end_int; ++i)
    {
        this->cells[i]->assembleQextGaussPoints(globalExternalHeat);
    }
    for (auto group : boundaryGroups)
    {
        group.second->assembleExternalHeat(globalExternalHeat);
    }
//     cout << globalExternalHeat << endl;
}

/**
 * @brief Sets the temperature of every node in the body.
 * @param temp_in Temperature value to assign to all nodes.
 */
void Body::setTemperature(double temp_in)
{
    for (auto& node : nodes)
    {
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
    if (isThermal && nodes.size() != 0)
    {
        temperature.push_back(new lmx::Vector<data_type>(nodes.size())); //temperature
        for (auto i = 0u; i < nodes.size(); ++i)
        {
            temperature.back()->writeElement(nodes[i]->getqt(), i);
        }
    }

    if (computeEnergy)   // TODO: store thermal energy
    {
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

    if (boundaryConnectivity.size() > 0)
    {
        *outFile << "BOUNDARY " << title << " " << boundaryConnectivity.size() << endl;
        for (auto& boundary : boundaryConnectivity)
        {
            for (auto& segment : boundary)
            {
                *outFile << segment << " ";
            }
            *outFile << endl;
        }
    }

    if (temperature.size() != 0)
    {
        *outFile << "TEMPERATURE " << title << endl;
        for (auto& temp : temperature)
        {
            auto  vectorSize = temp->size();
            for (auto i = 0u; i < vectorSize; ++i)
            {
                *outFile << temp->readElement(i) << " ";
            }
            *outFile << endl;
        }
    }
    outputVTK();
}

/**
 * @brief Prepare for saving VTK files.
 *
 * @return void
 **/
void Body::outputVTK( )
{
#ifdef HAVE_VTK
    std::string outputDir = "./" + this->title + "/";
    
    // Create output directory
    struct stat info;
    if (stat(outputDir.c_str(), &info) != 0) {
        if (mkdir(outputDir.c_str(), 0777) != 0) {
            std::cerr << "Error: Cannot create directory " << outputDir << std::endl;
            return;
        }
    }

    cout << "VTK OUTPUT: Writing " << this->temperature.size() << " timesteps for " << this->title << endl;
    
    // Build mapping from global node number to local body index
    std::map<int, int> globalToLocalNodeMap;
    for (size_t i = 0; i < this->nodes.size(); ++i) {
        int globalNum = this->nodes[i]->getNumber();
        globalToLocalNodeMap[globalNum] = i;
    }
    
    std::vector<std::pair<int, std::string>> pvd_entries;
    
    // Write each timestep
    for (size_t t = 0; t < this->temperature.size(); ++t) {
        vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        
        // Create points from body nodes
        vtkSmartPointer<vtkPoints> vtkpoints = vtkSmartPointer<vtkPoints>::New();
        for (auto node : this->nodes) {
            vtkpoints->InsertNextPoint(node->getConf(0), node->getConf(1), node->getConf(2));
        }
        uGrid->SetPoints(vtkpoints);
        
        // Add cells with proper node mapping
        size_t numCells = 0;
        for (auto& cellPair : this->cells) {
            Cell* cell = cellPair.second;
            std::vector<int> globalNodeNumbers = cell->getNodeNumbers();
            
            // Map global node numbers to local body indices
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            for (int globalNum : globalNodeNumbers) {
                if (globalToLocalNodeMap.find(globalNum) != globalToLocalNodeMap.end()) {
                    idList->InsertNextId(globalToLocalNodeMap[globalNum]);
                }
            }
            
            // Only add cell if all nodes were found in mapping
            if ((int)idList->GetNumberOfIds() == (int)globalNodeNumbers.size()) {
                int cellType = VTK_EMPTY_CELL;
                if (globalNodeNumbers.size() == 4) cellType = VTK_TETRA;
                else if (globalNodeNumbers.size() == 3) cellType = VTK_TRIANGLE;
                else if (globalNodeNumbers.size() == 8) cellType = VTK_HEXAHEDRON;
                else if (globalNodeNumbers.size() == 6) cellType = VTK_WEDGE;
                
                if (cellType != VTK_EMPTY_CELL) {
                    uGrid->InsertNextCell(cellType, idList);
                    numCells++;
                }
            }
        }
        
        // Add temperature data
        vtkSmartPointer<vtkFloatArray> tempArray = vtkSmartPointer<vtkFloatArray>::New();
        tempArray->SetName("Temperature");
        tempArray->SetNumberOfComponents(1);
        tempArray->SetNumberOfTuples(this->nodes.size());
        
        for (size_t i = 0; i < this->nodes.size(); ++i) {
            float value = 0.0f;
            if (i < this->temperature[t]->size()) {
                value = this->temperature[t]->readElement(i);
            }
            tempArray->SetValue(i, value);
        }
        uGrid->GetPointData()->SetScalars(tempArray);
        
        // Write VTU file
        std::ostringstream filename;
        filename << "step_" << std::setw(3) << std::setfill('0') << t << ".vtu";
        std::string filepath = outputDir + filename.str();
        
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(filepath.c_str());
        writer->SetInputData(uGrid);
        writer->SetDataModeToAppended();
        writer->SetCompressorTypeToZLib();
        writer->EncodeAppendedDataOff();
        
        if (!writer->Write()) {
            std::cerr << "ERROR writing: " << filepath << std::endl;
        }
        
        pvd_entries.push_back({(int)t, filename.str()});
    }
    
    // Write PVD file
    std::string pvd_filename = outputDir + "output.pvd";
    std::ofstream pvd_file(pvd_filename.c_str());
    pvd_file << "<?xml version=\"1.0\"?>\n";
    pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd_file << "  <Collection>\n";
    for (auto& entry : pvd_entries) {
        pvd_file << "    <DataSet timestep=\"" << entry.first 
                << "\" group=\"\" part=\"0\" file=\"" << entry.second << "\"/>\n";
    }
    pvd_file << "  </Collection>\n";
    pvd_file << "</VTKFile>\n";
    pvd_file.close();
    
    cout << "VTK OUTPUT completed for " << this->title << endl;

#endif
}

/**
 * @brief Appends a connectivity record (list of node indices) to the boundary connectivity table.
 * @param connectivity_in Vector of node indices defining one boundary segment.
 */
void Body::addBoundaryConnectivity(std::vector<int> connectivity_in)
{
    this->boundaryConnectivity.push_back(std::vector<int>(connectivity_in));
}

/**
 * @brief Translates all nodes of the body by the given displacement vector.
 * @param x_in Displacement along the X axis.
 * @param y_in Displacement along the Y axis.
 * @param z_in Displacement along the Z axis.
 */
void Body::translate(double x_in, double y_in, double z_in)
{
    for (auto node : nodes)
    {
        node->setX(node->getX() + x_in);
        node->setY(node->getY() + y_in);
        node->setZ(node->getZ() + z_in);
    }
}

constexpr double deg2rad(double deg)
{
    return deg * M_PI / 180.0;
}

/**
 * @brief Rotates all nodes of the body using Euler angles (ZYX convention).
 * @param phi   Rotation angle about the X axis (degrees).
 * @param theta Rotation angle about the Y axis (degrees).
 * @param psi   Rotation angle about the Z axis (degrees).
 */
void Body::rotate(double phi, double theta, double psi)
{
    phi = deg2rad(phi);
    theta = deg2rad(theta);
    psi = deg2rad(psi);

    lmx::DenseMatrix<double> r_x(3, 3);
    r_x.writeElement(1, 0, 0);
    r_x.writeElement(cos(phi), 1, 1);
    r_x.writeElement(-sin(phi), 2, 1);
    r_x.writeElement(sin(phi), 1, 2);
    r_x.writeElement(cos(phi), 2, 2);

    lmx::DenseMatrix<double> r_y(3, 3);
    r_y.writeElement(cos(theta), 0, 0);
    r_y.writeElement(sin(theta), 2, 0);
    r_y.writeElement(1, 1, 1);
    r_y.writeElement(-sin(theta), 0, 2);
    r_y.writeElement(cos(theta), 2, 2);

    lmx::DenseMatrix<double> r_z(3, 3);
    r_z.writeElement(cos(psi), 0, 0);
    r_z.writeElement(-sin(psi), 1, 0);
    r_z.writeElement(sin(psi), 0, 1);
    r_z.writeElement(cos(psi), 1, 1);
    r_z.writeElement(1, 2, 2);

    auto r = r_z * r_y * r_x;

    lmx::Vector<double> v(3);

    for (auto node : nodes)
    {
        v.writeElement(node->getX(), 0);
        v.writeElement(node->getY(), 1);
        v.writeElement(node->getZ(), 2);
        auto v2 = r * v;
        node->setX(v2(0));
        node->setY(v2(1));
        node->setZ(v2(2));
    }
}

}
