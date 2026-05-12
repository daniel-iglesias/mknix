/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   https://github.com/daniel-iglesias/mknix                                        *
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

#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_VTK
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#endif 

namespace mknix
{

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
    // Create VTK file for visualization
    std::stringstream ss;
    ss << this->title << ".vtu";
    auto outFileNameVTK = ss.str();

    std::string outputDir = "./" + this->title + "/";

   // Check and create output directory if needed
    struct stat info;
    if (stat(outputDir.c_str(), &info) != 0) {
        if (mkdir(outputDir.c_str(), 0777) != 0) {
            std::cerr << "Error: Cannot create directory " << outputDir << ": " << std::strerror(errno) << std::endl;
            return;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "Error: " << outputDir << " exists but is not a directory." << std::endl;
        return;
    }


    // Save nodes in VTK variables 
    vtkPoints * vtkpoints = vtkPoints::New();
    
    // Create an unstructured grid and add points
    vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create an the writer
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    cout << "VTK OUTPUT: Initial configuration." << endl;
    for (auto& point : nodes)
    {
        vtkpoints->InsertNextPoint(point->getConf(0),
                                point->getConf(1),
                                point->getConf(2));
    }

    // TBD : Save mesh in VTK cells variables
    // baseSystem->writeFlexBodies(&outFile);

    uGrid->SetPoints(vtkpoints);

    for (auto& cell : cells) {
        std::vector<int> nodeNumbers = cell.second->getNodeNumbers();
        // Create a VTK cell (e.g., tetrahedron, hexahedron) based on nodeNumbers
        // This part depends on the type of cell and needs to be implemented accordingly
        // We assume it's a TETRA cell for the moment:
        if ( nodeNumbers.size() == 4 ) {
            vtkSmartPointer<vtkTetra> quad = vtkSmartPointer<vtkTetra>::New();
            for (vtkIdType i = 0; i < nodeNumbers.size(); ++i) {
                quad->GetPointIds()->SetId(i, nodeNumbers[i]);
            }
            uGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
        }
        else if ( nodeNumbers.size() == 3 ) {
            // Handle triangle case
            // Similar to tetrahedron, but using vtkTriangle
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
            for (vtkIdType i = 0; i < nodeNumbers.size(); ++i) {
                triangle->GetPointIds()->SetId(i, nodeNumbers[i]);
            }
            uGrid->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
        }
    }
    std::vector<std::pair<int, std::string>> pvd_entries;
    int t=0;
    for (auto& temp : temperature)
    {
        // Create temperature array
        vtkSmartPointer<vtkFloatArray> temperature = vtkSmartPointer<vtkFloatArray>::New();
        temperature->SetName("Temperature");
        temperature->SetNumberOfComponents(1);
        auto  vectorSize = temp->size();
        temperature->SetNumberOfTuples(vectorSize);
        for (auto i = 0u; i < vectorSize; ++i) {
            temperature->SetValue(i, temp->readElement(i));
        }

        uGrid->GetPointData()->SetScalars(temperature);

        // Format filename
        std::ostringstream filename;
        filename << "step_" << std::setw(3) << std::setfill('0') << t << ".vtu";
        std::string filepath = "./" + this->title + "/" + filename.str();

        // Write .vtu file
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(filepath.c_str());
        writer->SetInputData(uGrid);
        writer->Write();

        pvd_entries.emplace_back(t, filename.str());
        
        ++t;
    }

    // Write to VTU file
    writer->SetFileName(outFileNameVTK.c_str());
    writer->SetInputData(uGrid);
    writer->Write();

    // Write .pvd file
    std::string pvd_filename = outputDir + "/output.pvd";
    std::ofstream pvd_file(pvd_filename.c_str());
    if (!pvd_file) {
        std::cerr << "Error: Could not open " << pvd_filename << " for writing." << std::endl;
        return;
    }

    pvd_file << "<?xml version=\"1.0\"?>\n";
    pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd_file << "  <Collection>\n";

    for (size_t i = 0; i < pvd_entries.size(); ++i) {
        pvd_file << "    <DataSet timestep=\"" << pvd_entries[i].first
                << "\" group=\"\" part=\"0\" file=\"" << pvd_entries[i].second << "\"/>\n";
    }

    pvd_file << "  </Collection>\n";
    pvd_file << "</VTKFile>\n";
    pvd_file.close();


    
#endif
}

void Body::addBoundaryConnectivity(std::vector<int> connectivity_in)
{
    this->boundaryConnectivity.push_back(std::vector<int>(connectivity_in));
}


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
