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
#include "readerflex.h"

#include <simulation/simulation.h>
#include <system/bodyflexglobalgalerkin.h>
#include <core/cellrect.h>
#include <core/celltriang.h>
#include <core/celltetrahedron.h>
#include <core/cellboundarylinear.h>
#include <core/elemtriangle.h>
#include <core/elemtetrahedron.h>
#include <system/system.h>


namespace mknix {

ReaderFlex::ReaderFlex()
        : theSimulation(nullptr)
        , output(nullptr)
        , input(nullptr)
{
}

ReaderFlex::ReaderFlex(Simulation * simulation_in,
                       std::ofstream& output_in,
                       std::ifstream& input_in
)
        : theSimulation(simulation_in)
        , output(&output_in)
        , input(&input_in)
{
}


ReaderFlex::~ReaderFlex()
{
}


} // namespace mknix


void mknix::ReaderFlex::readFlexBodies(System * system_in)
{
    std::string keyword, bodyType;
    std::string boundaryType("CLOCKWISE");
    std::string stringKeyword, flexTitle;
    std::string meshfreeFormulation("RPIM");

    while (*input >> bodyType) {
        if (bodyType == "ENDFLEXBODIES") {
            break;
        } else if (bodyType == "SHARENODES") {
            std::string bodyTitleA, bodyTitleB;
            *input >> bodyTitleA >> bodyTitleB;
            *output << "SHARENODES: " << bodyTitleA << " " << bodyTitleB << endl;
            if (system_in->flexBodies.find(bodyTitleA) !=
                system_in->flexBodies.end()) {
                if (system_in->flexBodies.find(bodyTitleB) !=
                    system_in->flexBodies.end()) {
                    system_in->flexBodies[bodyTitleA]->addNodes(system_in->flexBodies[bodyTitleB]->getNodes());
                    system_in->flexBodies[bodyTitleB]->addNodes(system_in->flexBodies[bodyTitleA]->getNodes());
                }
                else {
                    cerr << "ERROR: Body " << bodyTitleB << " not found!!" << std::endl;
                    return;
                }
            }
            else {
                cerr << "ERROR: Body " << bodyTitleA << " not found!!" << std::endl;
                return;
            }
        }
        else if (bodyType == "MESHFREE" || bodyType == "FEMESH") {
            char a;
            *input >> flexTitle;
            *output << "FLEXBODY: "
            << system_in->getTitle()
            << "."
            << flexTitle << std::endl;
            system_in->flexBodies[flexTitle]
                    = new FlexGlobalGalerkin(flexTitle);
            system_in->flexBodies[flexTitle]->setType(bodyType);
            system_in->thermalBodies[flexTitle] = system_in->flexBodies[flexTitle];

            while (*input >> keyword) {
// 	      cout << "keyword: " << keyword << endl;
                if (keyword == "ENDMESHFREE" || keyword == "ENDFEMESH") {
                    break;
                } else if (keyword == "OUTPUT") {
                    *input >> stringKeyword;
                    if (stringKeyword == "STRESS") {
                        system_in->flexBodies[flexTitle]->setOutput(stringKeyword);
                        *output << "OUTPUT: " << stringKeyword << endl;
                    }
                    if (stringKeyword == "ENERGY") {
                        system_in->flexBodies[flexTitle]->setOutput(stringKeyword);
                        *output << "OUTPUT: " << stringKeyword << endl;
                    }
                }
                else if (keyword == "FORMULATION") {
                    *input >> stringKeyword;
                    system_in->flexBodies[flexTitle]->setFormulation(stringKeyword);
                    *output << "FORMULATION: " << stringKeyword << endl;
                    if (bodyType == "MESHFREE") {
                    }
                }
                else if (keyword == "INITIALTEMPERATURE") {
                    double init_temp;
                    *input >> init_temp; // Type of formulation: EFG or RPIM
                    *output << "INITIAL TEMPERATURE SET TO: "
                    << init_temp
                    << std::endl;
                    // BUG: No effect as it's overriden by Simulation class setting:
                    system_in->flexBodies[flexTitle]->setTemperature(init_temp);
                }
                else if (keyword == "METHOD") {
                    *input >> meshfreeFormulation; // Type of formulation: EFG or RPIM
                    *output << "METHOD SET TO: "
                    << meshfreeFormulation
                    << std::endl;
                }
                else if (keyword == "LAYER") { // Needs to be called after TRIANGLES or TETRAHEDRONS
                    int newMaterial, counter(0);
                    double thickness;
                    *input >> newMaterial >> thickness;
                    *output << "Layer: Material = " << newMaterial << '\t Thickness = ' << thickness << endl;
                    for(auto& p_cell : system_in->flexBodies[flexTitle]->getCells() ){
                        counter += p_cell.second->setMaterialIfLayer(theSimulation->materials[newMaterial], thickness);
                    }
                    *output << '\t' << "Number of cells changed: " << counter << endl;                        
                }
                else if (keyword == "BOUNDARYGROUP") {
                    std::string boundaryFormulation;
                    std::string boundaryName;
                    double nGPs, alpha;
                    *input >> boundaryName;
                    system_in->flexBodies[flexTitle]->addBoundaryGroup(boundaryName);
                    *output << "BOUNDARYGROUP: " << keyword << endl;
                    while (*input >> keyword) {
                        if (keyword == "ENDBOUNDARYGROUP") {
                            break;
                        } else if (keyword == "METHOD") {
                            *input >> boundaryFormulation;
                            *output << "METHOD: " << boundaryFormulation << endl;
                            *input >> nGPs >> alpha;
                            *output << "Number of GPs per cell: " << nGPs << endl;
                            *output << "Influence domain factor (alpha): " << alpha << endl;
                        }
                        else if (keyword == "FILE") {
                            int totalNodes, totalCells;
                            *input >> keyword;
                            *output << "FILE: " << keyword << endl;
                            std::ifstream meshfile(keyword.c_str()); // file to read points from
                            int nodeNumberInFile; // read number but do not use it
                            double x, y, z;
                            meshfile >> totalNodes >> totalCells;
                            *output << '\t' << "Nodes: " << totalNodes << endl;
                            *output << "\tNode\tx\ty\tz \n";
                            for (int i = 0; i < totalNodes; ++i) {
                                meshfile >> nodeNumberInFile >> x >> y >> z;
                                system_in->flexBodies[flexTitle]->addNodeToBoundaryGroup(nodeNumberInFile - 1,
                                                                                         boundaryName);
                                *output
                                << '\t' << system_in->flexBodies[flexTitle]->getNode(nodeNumberInFile - 1)->getNumber()
                                << std::endl;
                            }
                            int elementType, node1, node2, cellCount;
                            *output << '\t' << "Cells: " << totalCells << endl;
                            *output << "\tCell\tnode1\tnode2\tnode3 \n";
                            //             old_size = theSimulation->cells.end()->first;
                            cellCount = 0;
                            for (int i = 0; i < totalCells; ++i) {
                                meshfile >> nodeNumberInFile >> elementType;
                                if (elementType == 102) {
                                    meshfile >> node1 >> node2;
                                    system_in->flexBodies[flexTitle]
                                            ->addCellToBoundaryGroup(new CellBoundaryLinear
                                                                             (boundaryFormulation,
                                                                              alpha,
                                                                              nGPs,
                                                                              system_in->flexBodies[flexTitle]->getNode(
                                                                                      node1 - 1),
                                                                              system_in->flexBodies[flexTitle]->getNode(
                                                                                      node2 - 1)
                                                                             ), boundaryName
                                            );
                                    *output
                                    << '\t'
                                    << /*old_size +*/ cellCount
                                    << "\t"
                                    << system_in->flexBodies[flexTitle]->getNode(node1 - 1)->getNumber()
                                    << "\t"
                                    << system_in->flexBodies[flexTitle]->getNode(node2 - 1)->getNumber()
                                    << std::endl;
                                }
                            }

                        }
                    }
                }
                else if (keyword == "NODES") {
                    do {
                        input->get(a);
                    }
                    while (a != '\n');
                    double x, y;
                    double z = 0.0; // Initialized for 2-D problems.
                    *output << "NODES:" << endl;
                    *input >> keyword;
                    if (keyword == "RECTANGULAR") {
                        int nx, ny; // number of cells in x and y direction
                        double x1, y1, x2, y2; // begin and end of rectangle
                        double Ax, Ay;
                        int old_size = theSimulation->nodes.end()->first;
                        int thermal_old_size = theSimulation->thermalNodes.end()->first;
                        int newNumber, newThermalNumber;
                        *input >> nx >> ny >> x1 >> y1 >> x2 >> y2;
                        *output << '\t' << "Rectangular patch: " << nx << "x" << ny << endl;
                        Ax = (x2 - x1) / nx;
                        Ay = (y2 - y1) / ny;
                        *output << "\tNode\tx\ty\tz \n";
                        for (int j = 0; j <= ny; ++j) {
                            for (int i = 0; i <= nx; ++i) {
                                newNumber = old_size + (ny + 1) * i + j;
                                newThermalNumber = thermal_old_size + (ny + 1) * i + j;
                                theSimulation->nodes[newNumber] =
                                        new Node(newNumber, x1 + i * Ax, y1 + j * Ay, z);
                                system_in->flexBodies[flexTitle]
                                        ->addNode(theSimulation->nodes[newNumber]);
                                theSimulation->nodes[newNumber]
                                        ->setThermalNumber(newThermalNumber);
                                theSimulation->thermalNodes[newThermalNumber]
                                        = theSimulation->nodes[newNumber];
                                *output
                                << '\t'
                                << theSimulation->nodes[newNumber]->getNumber()
                                << "\t"
                                << theSimulation->nodes[newNumber]->getX()
                                << "\t"
                                << theSimulation->nodes[newNumber]->getY()
                                << "\t"
                                << theSimulation->nodes[newNumber]->getZ()
                                << std::endl;
                            }
                        }
                        // build the boundary lines:
                        for (int i = 0; i < nx; ++i) { //bottom line:
                            system_in->flexBodies[flexTitle]
                                    ->addBoundaryLine(system_in->flexBodies[flexTitle]
                                                              ->getNode((ny + 1) * i),
                                                      system_in->flexBodies[flexTitle]
                                                              ->getNode((ny + 1) * (i + 1))
                                    );
                        }
                        for (int i = 0; i < ny; ++i) { //right line:
                            system_in->flexBodies[flexTitle]
                                    ->addBoundaryLine(system_in->flexBodies[flexTitle]
                                                              ->getNode(nx * (ny + 1) + i),
                                                      system_in->flexBodies[flexTitle]
                                                              ->getNode(nx * (ny + 1) + i + 1)
                                    );
                        }
                        for (int i = nx; i > 0; --i) { //top line:
                            system_in->flexBodies[flexTitle]
                                    ->addBoundaryLine(system_in->flexBodies[flexTitle]
                                                              ->getNode((ny + 1) * (i + 1) - 1),
                                                      system_in->flexBodies[flexTitle]
                                                              ->getNode((ny + 1) * (i) - 1)
                                    );
                        }
                        for (int i = ny; i > 0; --i) { //left line:
                            system_in->flexBodies[flexTitle]
                                    ->addBoundaryLine(system_in->flexBodies[flexTitle]
                                                              ->getNode(i),
                                                      system_in->flexBodies[flexTitle]
                                                              ->getNode(i - 1)
                                    );
                        }

                        do {
                            input->get(a);
                        }
                        while (a != '\n');
                    }
                    else if (keyword == "BLOCK") {
                        //           int nx, ny, nz; // number of cells in x and y direction
                        //           double x1,y1, x2,y2, z1, z2; // begin and end of rectangle
                        //           double Ax, Ay, Az;
                        //           int old_size = nodes.end()->first;
                        //           *input >> nx >> ny >> nz >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
                        //           *output << '\t' << "Block patch: " << nx << " x " << ny << " x "
                        //            << nz << endl;
                        //           Ax = (x2 - x1) / nx;
                        //           Ay = (y2 - y1) / ny;
                        //           Az = (z2 - z1) / nz;
                        //           *output << "\tNode\tx\ty\tz \n";
                        //           for (int i=0; i <= nx; ++i){
                        //             for (int j=0; j <= ny; ++j){
                        //               for (int k=0; k <= nz; ++k){
                        //                 nodes[old_size + (ny+1)*(nz+1)*i + (nz+1)*j + k] =
                        //                   *(new Node(old_size + (ny+1)*(nz+1)*i + (nz+1)*j + k,
                        //                     x1 + i*Ax , y1 + j*Ay, z1 + k*Az) );
                        //                 *output
                        //                   << '\t'
                        //                   << nodes[old_size + (ny+1)*(nz+1)*i + (nz+1)*j + k]
                        //                      .nodeNumber()
                        //                   << "\t"
                        //                   << nodes[old_size + (ny+1)*(nz+1)*i + (nz+1)*j + k]
                        //                      .getx()
                        //                   << "\t"
                        //                   << nodes[old_size + (ny+1)*(nz+1)*i + (nz+1)*j + k]
                        //                      .gety()
                        //                   << "\t"
                        //                   << nodes[old_size + (ny+1)*(nz+1)*i + (nz+1)*j + k]
                        //                      .getz()
                        //                   << std::endl;
                        //               }
                        //             }
                        //           }
                        //           do {input->get(a);} while(a!='\n');
                    }
                    else if (keyword == "GRID") {
                        char a;
                        do {
                            input->get(a);
                        }
                        while (a != '\n');
                        int totalNodes, totalCells;
                        //          int nodeNumber;
                        //          double x, y;
                        //          double z = 0.0; // Initialized for 2-D problems.
                        *output << "NODES GRID:" << endl;
                        *input >> keyword;
                        if (keyword == "TRIANGLES") {
                            *input >> keyword;
                            *output << "FILE: " << keyword << endl;
                            std::ifstream meshfile(keyword.c_str()); // file to read points from
                            int nodeNumberInFile; // read number but do not use it
// 			  double x, y, z;
                            int old_size = theSimulation->nodes.end()->first;
                            int thermal_old_size = theSimulation->thermalNodes.end()->first;
                            meshfile >> totalNodes >> totalCells;
                            *output << '\t' << "Nodes: " << totalNodes << endl;
                            *output << "\tNode\tx\ty\tz \n";
                            for (int i = 0; i < totalNodes; ++i) {
                                meshfile >> nodeNumberInFile >> x >> y >> z;
                                theSimulation->nodes[old_size + i] =
                                        new Node(old_size + i, x, y, z);
                                theSimulation->nodes[old_size + i]->setThermalNumber(thermal_old_size + i);
                                theSimulation->thermalNodes[thermal_old_size + i]
                                        = theSimulation->nodes[old_size + i];
                                system_in->flexBodies[flexTitle]
                                        ->addNode(theSimulation->nodes[old_size + i]);
                                *output
                                << '\t' << theSimulation->nodes[old_size + i]->getNumber()
                                << "\t" << theSimulation->nodes[old_size + i]->getX()
                                << "\t" << theSimulation->nodes[old_size + i]->getY()
                                << "\t" << theSimulation->nodes[old_size + i]->getZ()
                                << std::endl;
                            }
                        }
                        else if (keyword == "TETRAHEDRONS") {
                            *input >> keyword;
                            *output << "FILE: " << keyword << endl;
                            std::ifstream meshfile(keyword.c_str()); // file to read points from
                            int nodeNumberInFile; // read number but do not use it
// 			  double x, y, z;
                            int old_size = theSimulation->nodes.end()->first;
                            int thermal_old_size = theSimulation->thermalNodes.end()->first;
                            std::map<int, int> nodesFile;
                            meshfile >> totalNodes >> totalCells;
                            *output << '\t' << "Nodes: " << totalNodes << endl;
                            *output << "\tNode\tx\ty\tz \n";
                            for (int i = 0; i < totalNodes; ++i) {
                                meshfile >> nodeNumberInFile >> x >> y >> z;
                                nodesFile[nodeNumberInFile] = old_size + i;
                                theSimulation->nodes[old_size + i] =
                                        new Node(old_size + i, x, y, z);
                                theSimulation->nodes[old_size + i]->setThermalNumber(thermal_old_size + i);
                                theSimulation->thermalNodes[thermal_old_size + i]
                                        = theSimulation->nodes[old_size + i];
                                theSimulation->outputPoints[old_size + i]
                                        = theSimulation->nodes[old_size + i];
                                system_in->flexBodies[flexTitle]
                                        ->addNode(theSimulation->nodes[old_size + i]);
                                *output
                                << '\t' << theSimulation->nodes[old_size + i]->getNumber()
                                << "\t" << theSimulation->nodes[old_size + i]->getX()
                                << "\t" << theSimulation->nodes[old_size + i]->getY()
                                << "\t" << theSimulation->nodes[old_size + i]->getZ()
                                << std::endl;
                            }
                        }
                    }
                }
                else if (keyword == "CELLS") {
                    char a;
                    do {
                        input->get(a);
                    }
                    while (a != '\n');
                    int material(1), nGPs(1);
                    double alpha(2.5), dc;
                    int totalPoints, totalCells;
//          int nodeNumber;
//          double x, y;
//          double z = 0.0; // Initialized for 2-D problems.
                    *output << "CELLS:" << endl;
                    *input >> material >> nGPs >> alpha >> dc;
                    *output << "Material: " << material << endl;
                    *output << "Number of GPs per cell: " << nGPs << endl;
                    *output << "Influence domain factor (alpha): " << alpha << endl;
                    *output << "Average node distance in domain (dc): " << dc << endl;
                    do {
                        input->get(a);
                    }
                    while (a != '\n');

                    *input >> keyword;
                    if (keyword == "RECTANGULAR") {
                        int nx, ny; // number of cells in x and y direction
                        double x1, y1, x2, y2; // begin and end of rectangle
                        double Ax, Ay, dcx, dcy;
                        int material_in;
                        //           int old_size = cells.size();
                        *input >> material_in
                        >> dcx >> dcy
                        >> nx >> ny
                        >> x1 >> y1
                        >> x2 >> y2;
                        *output << '\t' << "Rectangular patch: " << "material: "
                        << material_in << " dimension: " << nx << " x " << ny
                        << endl;
                        Ax = (x2 - x1) / nx;
                        Ay = (y2 - y1) / ny;
                        for (int i = 0; i < nx; ++i) {
                            for (int j = 0; j < ny; ++j) {
                                system_in->flexBodies[flexTitle]
                                        ->addCell(ny * i + j,
                                                  new CellRect(
                                                          theSimulation->materials[material_in],
                                                          meshfreeFormulation,
                                                          alpha,
                                                          nGPs,
                                                          x1 + i * Ax, y1 + j * Ay,
                                                          x1 + (i + 1) * Ax, y1 + (j) * Ay,
                                                          x1 + (i + 1) * Ax, y1 + (j + 1) * Ay,
                                                          x1 + (i) * Ax, y1 + (j + 1) * Ay,
                                                          dcx, dcy,
                                                          x1, y1,
                                                          x2, y2
                                                  )
                                        );
                                *output << "\t\t" << "Cell:" << /*old_size +*/ ny * i + j
                                << endl
                                << "\t\t\t(" << x1 + i * Ax << ", " << y1 + j * Ay << ")"
                                << endl
                                << "\t\t\t(" << x1 + (i + 1) * Ax << ", " << y1 + j * Ay << ")"
                                << endl
                                << "\t\t\t(" << x1 + (i + 1) * Ax << ", " << y1 + (j + 1) * Ay << ")"
                                << endl
                                << "\t\t\t(" << x1 + i * Ax << ", " << y1 + (j + 1) * Ay << ")"
                                << endl;
                            }
                        }
                        do {
                            input->get(a);
                        }
                        while (a != '\n');
                    }
                    else if (keyword == "BLOCK") {
                        //           int nx, ny, nz; // number of cells in x, y and z direction
                        //           double x1,y1, x2,y2, z1, z2; // begin and end of patch
                        //           double Ax, Ay, Az, dcx, dcy, dcz;
                        //           int material_in;
                        //           int old_size = cells.size();
                        //           *input >> material_in
                        //                 >> dcx >> dcy >> dcz
                        //                 >> nx >> ny >> nz
                        //                 >> x1 >> y1 >> z1
                        //                 >> x2 >> y2 >> z2;
                        //           *output << '\t' << "Block patch: " << "material: " << material_in
                        //                  << " dimension: " << nx << " x " << ny <<  " x " << nz
                        //                  << endl;
                        //           Ax = (x2 - x1) / nx;
                        //           Ay = (y2 - y1) / ny;
                        //           Az = (z2 - z1) / nz;
                        //           for (int i=0; i < nx; ++i){
                        //             for (int j=0; j < ny; ++j){
                        //               for (int k=0; k < nz; ++k){
                        //                 cells[old_size + (ny+nz)*i + nz*j + k] =
                        //                   new CellRect( materials[material_in],
                        //                                 alpha,
                        //                                 gPoints,
                        //                                 x1 + (i)*Ax ,   y1 + (j)*Ay,   z1 + (k)*Az,
                        //                                 x1 + (i+1)*Ax , y1 + (j)*Ay,   z1 + (k)*Az,
                        //                                 x1 + (i+1)*Ax , y1 + (j+1)*Ay, z1 + (k)*Az,
                        //                                 x1 + (i)*Ax ,   y1 + (j+1)*Ay, z1 + (k)*Az,
                        //                                 x1 + (i)*Ax ,   y1 + (j)*Ay,   z1 + (k+1)*Az,
                        //                                 x1 + (i+1)*Ax , y1 + (j)*Ay,   z1 + (k+1)*Az,
                        //                                 x1 + (i+1)*Ax , y1 + (j+1)*Ay, z1 + (k+1)*Az,
                        //                                 x1 + (i)*Ax ,   y1 + (j+1)*Ay, z1 + (k+1)*Az,
                        //                                 dcx, dcy, dcz,
                        //                                 x1, y1, z1,
                        //                                 x2, y2, z2
                        //                               );
                        //               *output << "\t\t" << "Cell:" << old_size + ny*i + j
                        //                      << endl << "\t\t\t("
                        //                      << x1 + i*Ax << ", " << y1 + j*Ay << ")"
                        //                      << endl << "\t\t\t("
                        //                      << x1 + (i+1)*Ax << ", " << y1 + j*Ay << ")"
                        //                      << endl << "\t\t\t("
                        //                      << x1 + (i+1)*Ax << ", " << y1 + (j+1)*Ay << ")"
                        //                      << endl << "\t\t\t("
                        //                      << x1 + i*Ax << ", " << y1 + (j+1)*Ay << ")"
                        //                      << endl;
                        //               }
                        //             }
                        //           }
                        //           do {input->get(a);} while(a!='\n');
                    }
                    else if (keyword == "GRID") {
                        *input >> keyword;
                        if (keyword == "TRIANGLES") {
                            *input >> keyword;
                            *output << "FILE: " << keyword << endl;
                            std::ifstream meshfile(keyword.c_str()); // file to read points from
                            int pointNumberInFile; // read number but do not use it
                            double x, y, z;
                            int old_size = theSimulation->outputPoints.end()->first;
                            std::map<int, int> pointsFile;
                            meshfile >> totalPoints >> totalCells;
                            *output << '\t' << "Points: " << totalPoints << endl;
                            *output << "\tPoint\tx\ty\tz \n";
                            for (int i = 0; i < totalPoints; ++i) {
                                meshfile >> pointNumberInFile >> x >> y >> z;
                                pointsFile[pointNumberInFile] = old_size + i;
                                theSimulation->outputPoints[old_size + i] =
                                        new Point(Simulation::getDim(), old_size + i, x, y, z, alpha,
                                                  dc); // dc is set to 1.03 for the moment
                                system_in->flexBodies[flexTitle]
                                        ->addBodyPoint(theSimulation->outputPoints[old_size + i], meshfreeFormulation);
                                *output
                                << '\t' << theSimulation->outputPoints[old_size + i]->getNumber()
                                << "\t" << theSimulation->outputPoints[old_size + i]->getX()
                                << "\t" << theSimulation->outputPoints[old_size + i]->getY()
                                << "\t" << theSimulation->outputPoints[old_size + i]->getZ()
                                << std::endl;
                            }
                            int elementType, node1, node2, node3, cellCount;
                            *output << '\t' << "Cells: " << totalCells << endl;
                            *output << "\tCell\tnode1\tnode2\tnode3 \n";
                            //             old_size = theSimulation->cells.end()->first;
                            cellCount = system_in->flexBodies[flexTitle]->getCellLastNumber();
                            for (int i = 0; i < totalCells; ++i) {
                                meshfile >> pointNumberInFile >> elementType;
                                if (elementType == 203) {
                                    meshfile >> node1 >> node2 >> node3;
                                    if (bodyType == "MESHFREE") {
                                        system_in->flexBodies[flexTitle]
                                                ->addCell(/*old_size +*/ cellCount,
                                                                         new CellTriang
                                                                                 (theSimulation->materials[material],
                                                                                  meshfreeFormulation,
                                                                                  alpha,
                                                                                  nGPs,
                                                                                  theSimulation->outputPoints[pointsFile[node1]],
                                                                                  theSimulation->outputPoints[pointsFile[node2]],
                                                                                  theSimulation->outputPoints[pointsFile[node3]],
                                                                                  dc
                                                                                 )
                                                );
                                        *output
                                        << '\t'
                                        << /*old_size +*/ cellCount
                                        << "\t"
                                        << theSimulation->outputPoints[pointsFile[node1]]->getNumber()
                                        << "\t"
                                        << theSimulation->outputPoints[pointsFile[node2]]->getNumber()
                                        << "\t"
                                        << theSimulation->outputPoints[pointsFile[node3]]->getNumber()
                                        << std::endl;
                                    }
                                    else if (bodyType == "FEMESH") {
                                        std::cout << "ERROR: NOT POSSIBLE TO CREATE A FEMESH OF THIS TYPE" << std::endl;
                                    }

                                    ++cellCount;
                                }
                                else if (elementType == 102) {
                                    meshfile >> node1 >> node2;
                                    if (boundaryType == "CLOCKWISE") {
                                        system_in->flexBodies[flexTitle]->addBoundaryLine
                                                (system_in->flexBodies[flexTitle]->getBodyPoint(pointsFile[node1]),
                                                 system_in->flexBodies[flexTitle]->getBodyPoint(pointsFile[node2])
                                                );
                                    }
                                    else if (boundaryType == "COUNTERCLOCKWISE") {
                                        system_in->flexBodies[flexTitle]->addBoundaryLine
                                                (system_in->flexBodies[flexTitle]->getBodyPoint(pointsFile[node2]),
                                                 system_in->flexBodies[flexTitle]->getBodyPoint(pointsFile[node1])
                                                );
                                    }
                                    else {
                                        cerr << ":::ERROR: BOUNDARY ORIENTATION UNKNOWN:::" << endl;
                                    }
                                    *output << "\t Bound"
                                    << "\t" << node1
                                    << "\t" << node2
                                    << std::endl;
                                }
                                else {
                                    do {
                                        meshfile.get(a);
                                    }
                                    while (a != '\n');
                                }
                            }
                            do {
                                input->get(a);
                            }
                            while (a != '\n');
                        }
                        else if (keyword == "TETRAHEDRONS") {
                            *input >> keyword;
                            *output << "FILE: " << keyword << endl;
                            std::ifstream meshfile(keyword.c_str()); // file to read points from
                            int pointNumberInFile; // read number but do not use it
                            double x, y, z;
                            int old_size = theSimulation->outputPoints.end()->first;
                            std::map<int, int> pointsFile;
                            meshfile >> totalPoints >> totalCells;
                            *output << '\t' << "Nodes: " << totalPoints << endl;
                            *output << "\tNode\tx\ty\tz \n";
                            for (int i = 0; i < totalPoints; ++i) {
                                meshfile >> pointNumberInFile >> x >> y >> z;
                                pointsFile[pointNumberInFile] = old_size + i;
                                theSimulation->outputPoints[old_size + i] =
                                        new Point(Simulation::getDim(), old_size + i, x, y, z, alpha,
                                                  dc); // dc is set to 1.03 for the moment
                                system_in->flexBodies[flexTitle]
                                        ->addBodyPoint(theSimulation->outputPoints[old_size + i], meshfreeFormulation);
                                *output
                                << '\t' << theSimulation->outputPoints[old_size + i]->getNumber()
                                << "\t" << theSimulation->outputPoints[old_size + i]->getX()
                                << "\t" << theSimulation->outputPoints[old_size + i]->getY()
                                << "\t" << theSimulation->outputPoints[old_size + i]->getZ()
                                << std::endl;
                            }
                            int elementType, node1, node2, node3, node4, cellCount;
                            std::vector<int> nodes(3);
                            *output << '\t' << "Cells: " << totalCells << endl;
                            *output << "\tCell\tnode1\tnode2\tnode3 \n";
                            //             old_size = theSimulation->cells.end()->first;
                            cellCount = system_in->flexBodies[flexTitle]->getCellLastNumber();
                            for (int i = 0; i < totalCells; ++i) {
                                meshfile >> pointNumberInFile >> elementType;
                                if (elementType == 203) {
                                    meshfile >> node1 >> node2 >> node3;
                                    nodes[0] = --node1;
                                    nodes[1] = --node2;
                                    nodes[2] = --node3;
                                    system_in->flexBodies[flexTitle]->addBoundaryConnectivity(nodes);
                                    *output << "\t" << i << "\t"
                                    << nodes[0] << "\t"
                                    << nodes[1] << "\t"
                                    << nodes[2] << endl;
                                }
                                else if (elementType == 102) { //TODO: read linear boundary
                                    // at the moment we dump the information
                                    meshfile >> node1 >> node2;
                                }
                                else if (elementType == 304) {
                                    meshfile >> node1 >> node2 >> node3 >> node4;
                                    if (bodyType == "MESHFREE") {
                                        system_in->flexBodies[flexTitle]
                                                ->addCell(/*old_size +*/ cellCount,
                                                                         new CellTetrahedron
                                                                                 (theSimulation->materials[material],
                                                                                  meshfreeFormulation,
                                                                                  alpha,
                                                                                  nGPs,
                                                                                  theSimulation->outputPoints[pointsFile[node1]],
                                                                                  theSimulation->outputPoints[pointsFile[node2]],
                                                                                  theSimulation->outputPoints[pointsFile[node3]],
                                                                                  theSimulation->outputPoints[pointsFile[node4]]
                                                                                 )
                                                );
                                    }
                                    else if (bodyType == "FEMESH") {
                                        std::cout << "ERROR: NOT POSSIBLE TO CREATE A FEMESH OF THIS TYPE" << std::endl;
                                    }
                                    *output
                                    << '\t'
                                    << /*old_size +*/ cellCount
                                    << "\t"
                                    << theSimulation->outputPoints[pointsFile[node1]]->getNumber()
                                    << "\t"
                                    << theSimulation->outputPoints[pointsFile[node2]]->getNumber()
                                    << "\t"
                                    << theSimulation->outputPoints[pointsFile[node3]]->getNumber()
                                    << "\t"
                                    << theSimulation->outputPoints[pointsFile[node4]]->getNumber()
                                    << std::endl;
                                    ++cellCount;
                                }
                            }
                        }
                    }
                }
                else if (keyword == "MESH") {
                    cout << "bodyType " << bodyType << endl;
                    char a;
                    do {
                        input->get(a);
                    }
                    while (a != '\n');
                    int material, nGPs;
                    double alpha;
                    int totalNodes, totalCells;
//          int nodeNumber;
//          double x, y;
//          double z = 0.0; // Initialized for 2-D problems.
                    *output << "MESH:" << endl;
                    *input >> material >> nGPs >> alpha;
                    *output << "Material: " << material << endl;
                    *output << "Number of GPs per cell: " << nGPs << endl;
                    *output << "Influence domain factor (alpha): " << alpha << endl;
                    do { input->get(a); } while (a != '\n');
                    *input >> keyword;
                    if (keyword == "TRIANGLES") {
                        *input >> keyword;
                        *output << "FILE: " << keyword << endl;
                        std::ifstream meshfile(keyword.c_str()); // file to read points from
                        int nodeNumberInFile; // read number but do not use it
                        double x, y, z;
                        int old_size = theSimulation->nodes.end()->first;
                        int thermal_old_size = theSimulation->thermalNodes.end()->first;
                        std::map<int, int> nodesFile;
                        meshfile >> totalNodes >> totalCells;
                        *output << '\t' << "Nodes: " << totalNodes << endl;
                        *output << "\tNode\tx\ty\tz \n";
                        for (int i = 0; i < totalNodes; ++i) {
                            meshfile >> nodeNumberInFile >> x >> y >> z;
                            nodesFile[nodeNumberInFile] = old_size + i;
                            theSimulation->nodes[old_size + i] =
                                    new Node(old_size + i, x, y, z);
                            theSimulation->nodes[old_size + i]->setThermalNumber(thermal_old_size + i);
                            theSimulation->thermalNodes[thermal_old_size + i]
                                    = theSimulation->nodes[old_size + i];
                            theSimulation->outputPoints[old_size + i]
                                    = theSimulation->nodes[old_size + i];
                            system_in->flexBodies[flexTitle]
                                    ->addNode(theSimulation->nodes[old_size + i]);
                            *output
                            << '\t' << theSimulation->nodes[old_size + i]->getNumber()
                            << "\t" << theSimulation->nodes[old_size + i]->getX()
                            << "\t" << theSimulation->nodes[old_size + i]->getY()
                            << "\t" << theSimulation->nodes[old_size + i]->getZ()
                            << std::endl;
                        }
                        int elementType, node1, node2, node3, cellCount;
                        *output << '\t' << "Cells: " << totalCells << endl;
                        *output << "\tCell\tnode1\tnode2\tnode3 \n";
//             old_size = theSimulation->cells.end()->first;
                        cellCount = 0;
                        for (int i = 0; i < totalCells; ++i) {
                            meshfile >> nodeNumberInFile >> elementType;
                            if (elementType == 203) {
                                meshfile >> node1 >> node2 >> node3;
// 		  cout << "bodyType " << bodyType << endl;
                                if (bodyType == "MESHFREE") {
                                    system_in->flexBodies[flexTitle]
                                            ->addCell(/*old_size +*/ cellCount,
                                                                     new CellTriang
                                                                             (theSimulation->materials[material],
                                                                              meshfreeFormulation,
                                                                              alpha,
                                                                              nGPs,
                                                                              theSimulation->nodes[nodesFile[node1]],
                                                                              theSimulation->nodes[nodesFile[node2]],
                                                                              theSimulation->nodes[nodesFile[node3]]
                                                                             )
                                            );
                                    *output
                                    << '\t'
                                    << /*old_size +*/ cellCount
                                    << "\t"
                                    << theSimulation->nodes[nodesFile[node1]]->getNumber()
                                    << "\t"
                                    << theSimulation->nodes[nodesFile[node2]]->getNumber()
                                    << "\t"
                                    << theSimulation->nodes[nodesFile[node3]]->getNumber()
                                    << std::endl;
                                }
                                else if (bodyType == "FEMESH") {
                                    system_in->flexBodies[flexTitle]
                                            ->addCell( /*old_size +*/ cellCount,
                                                                      new ElemTriangle
                                                                              (theSimulation->materials[material],
                                                                               alpha,
                                                                               nGPs,
                                                                               theSimulation->nodes[nodesFile[node1]],
                                                                               theSimulation->nodes[nodesFile[node2]],
                                                                               theSimulation->nodes[nodesFile[node3]]
                                                                              )
                                            );
                                    *output
                                    << '\t'
                                    << /*old_size +*/ cellCount
                                    << "\t"
                                    << theSimulation->nodes[nodesFile[node1]]->getNumber()
                                    << "\t"
                                    << theSimulation->nodes[nodesFile[node2]]->getNumber()
                                    << "\t"
                                    << theSimulation->nodes[nodesFile[node3]]->getNumber()
                                    << std::endl;
                                }

                                ++cellCount;
                            }
                            else if (elementType == 102) {
                                std::vector<int> nodes(2);
                                meshfile >> node1 >> node2;
                                nodes[0] = --node1;
                                nodes[1] = --node2;
                                system_in->flexBodies[flexTitle]->addBoundaryConnectivity(nodes);
// 				*output << "\t" << i << "\t"
// 				<< nodes[0] << "\t"
// 				<< nodes[1] << "\t"
// 				<< nodes[2] << endl;
//                                 if(boundaryType == "CLOCKWISE") {
//                                     system_in->flexBodies[flexTitle]->addBoundaryLine
//                                     ( system_in->flexBodies[flexTitle]->getNode(node1-1),
//                                       system_in->flexBodies[flexTitle]->getNode(node2-1)
//                                     );
//                                 }
//                                 else if(boundaryType == "COUNTERCLOCKWISE") {
//                                     system_in->flexBodies[flexTitle]->addBoundaryLine
//                                     ( system_in->flexBodies[flexTitle]->getNode(node2-1),
//                                       system_in->flexBodies[flexTitle]->getNode(node1-1)
//                                     );
//                                 }
//                                 else
//                                     cerr << ":::ERROR: BOUNDARY ORIENTATION UNKNOWN:::" << endl;
                                *output << "\t Bound"
                                << "\t" << node1
                                << "\t" << node2
                                << std::endl;
                            }
                            else {
                                do {
                                    meshfile.get(a);
                                }
                                while (a != '\n');
                            }
                        }
                        do {
                            input->get(a);
                        }
                        while (a != '\n');
                    }
                    else if (keyword == "TETRAHEDRONS") {
                        *input >> keyword;
                        *output << "FILE: " << keyword << endl;
                        std::ifstream meshfile(keyword.c_str()); // file to read points from
                        int nodeNumberInFile; // read number but do not use it
                        double x, y, z;
                        int old_size = theSimulation->nodes.end()->first;
                        int thermal_old_size = theSimulation->thermalNodes.end()->first;
                        std::map<int, int> nodesFile;
                        meshfile >> totalNodes >> totalCells;
                        *output << '\t' << "Nodes: " << totalNodes << endl;
                        *output << "\tNode\tx\ty\tz \n";
                        for (int i = 0; i < totalNodes; ++i) {
                            meshfile >> nodeNumberInFile >> x >> y >> z;
                            nodesFile[nodeNumberInFile] = old_size + i;
                            theSimulation->nodes[old_size + i] =
                                    new Node(old_size + i, x, y, z);
                            theSimulation->nodes[old_size + i]->setThermalNumber(thermal_old_size + i);
                            theSimulation->thermalNodes[thermal_old_size + i]
                                    = theSimulation->nodes[old_size + i];
                            theSimulation->outputPoints[old_size + i]
                                    = theSimulation->nodes[old_size + i];
                            system_in->flexBodies[flexTitle]
                                    ->addNode(theSimulation->nodes[old_size + i]);
                            *output
                            << '\t' << theSimulation->nodes[old_size + i]->getNumber()
                            << "\t" << theSimulation->nodes[old_size + i]->getX()
                            << "\t" << theSimulation->nodes[old_size + i]->getY()
                            << "\t" << theSimulation->nodes[old_size + i]->getZ()
                            << std::endl;
                        }
                        int elementType, node1, node2, node3, node4, cellCount;
                        std::vector<int> nodes(3);
                        *output << '\t' << "Cells: " << totalCells << endl;
                        *output << "\tCell\tnode1\tnode2\tnode3 \n";
                        //             old_size = theSimulation->cells.end()->first;
                        cellCount = 0;
                        for (int i = 0; i < totalCells; ++i) {
                            meshfile >> nodeNumberInFile >> elementType;
                            //                             cout << elementType << endl;
                            if (elementType == 203) { //TODO: read facet boundary
                                // at the moment we dump the information
                                meshfile >> node1 >> node2 >> node3;
                                nodes[0] = --node1;
                                nodes[1] = --node2;
                                nodes[2] = --node3;
                                system_in->flexBodies[flexTitle]->addBoundaryConnectivity(nodes);
                                *output << "\t" << i << "\t"
                                << nodes[0] << "\t"
                                << nodes[1] << "\t"
                                << nodes[2] << endl;
                            }
                            else if (elementType == 102) { //TODO: read linear boundary
                                // at the moment we dump the information
                                meshfile >> node1 >> node2;
                            }
                            else if (elementType == 304) {
                                meshfile >> node1 >> node2 >> node3 >> node4;
                                if (bodyType == "MESHFREE") {
                                    system_in->flexBodies[flexTitle]
                                            ->addCell(/*old_size +*/ cellCount,
                                                                     new CellTetrahedron
                                                                             (theSimulation->materials[material],
                                                                              meshfreeFormulation,
                                                                              alpha,
                                                                              nGPs,
                                                                              theSimulation->nodes[nodesFile[node1]],
                                                                              theSimulation->nodes[nodesFile[node2]],
                                                                              theSimulation->nodes[nodesFile[node3]],
                                                                              theSimulation->nodes[nodesFile[node4]]
                                                                             )
                                            );
                                }
                                else if (bodyType == "FEMESH") {
                                    system_in->flexBodies[flexTitle]
                                            ->addCell( /*old_size +*/ cellCount,
                                                                      new ElemTetrahedron
                                                                              (theSimulation->materials[material],
                                                                               alpha,
                                                                               nGPs,
                                                                               theSimulation->nodes[nodesFile[node1]],
                                                                               theSimulation->nodes[nodesFile[node2]],
                                                                               theSimulation->nodes[nodesFile[node3]],
                                                                               theSimulation->nodes[nodesFile[node4]]
                                                                              )
                                            );
                                }
                                *output
                                << '\t'
                                << /*old_size +*/ cellCount
                                << "\t"
                                << theSimulation->nodes[nodesFile[node1]]->getNumber()
                                << "\t"
                                << theSimulation->nodes[nodesFile[node2]]->getNumber()
                                << "\t"
                                << theSimulation->nodes[nodesFile[node3]]->getNumber()
                                << "\t"
                                << theSimulation->nodes[nodesFile[node4]]->getNumber()
                                << std::endl;
                                ++cellCount;
                            }
                        }
                    }
                }
            }
        } // MESHFREE
    }
    system_in->initFlexBodies();
}
