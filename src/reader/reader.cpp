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

#include "reader.h"
#include "readerrigid.h"
#include "readerflex.h"
#include "readerconstraints.h"

#include <simulation/simulation.h>
#include <simulation/analysisdynamic.h>
#include <simulation/analysisstatic.h>
#include <simulation/analysisthermaldynamic.h>
#include <simulation/analysisthermalstatic.h>
#include <simulation/analysisthermomechanicaldynamic.h>
#include <system/bodyflex.h>
#include <system/bodyrigid.h>
#include <system/bodythermal.h>
#include <system/force.h>
#include <system/loadradiation.h>
#include <system/loadthermal.h>
#include <system/motion.h>
#include <system/loadthermalbody.h>
#include <system/loadthermalboundary1D.h>

namespace {
constexpr auto pathSep = "/";

std::string dirName(const std::string& path)
{
    auto found = path.find_last_of(pathSep);
    return path.substr(0, found);
}
}

namespace mknix {

Reader::Reader()
        : theSimulation(nullptr)
        , theReaderRigid(nullptr)
{
    output.open("output.reader");
}

Reader::Reader(Simulation * simulation_in)
        : theSimulation(simulation_in)
        , theReaderRigid(nullptr)
{
    output.open("output.reader");
}


Reader::~Reader()
{
}


} // namespace mknix

void mknix::Reader::inputFromFile(const std::string& fileIn)
{
  std::cout << "Opening " << fileIn << std::endl;
    input.open(fileIn, std::ifstream::in);
    if(!input.good()){
      std::cout << "Error! Could not read the file: " << fileIn << std::endl;
      throw std::runtime_error("Could not read the file: " + fileIn);
    }

    char keyword[100];

    // Change to directory of input file so all mesh file paths are relative to here
    auto dir = dirName(fileIn);
    chdir(dir.c_str());

    while (input >> keyword) {

        if (!strcmp(keyword, "//")) {
            char a;
            do {
                input.get(a);
            }
            while (a != '\n');
        }
        else if (!strcmp(keyword, "TITLE")) {
            input >> theSimulation->title;
            output << "TITLE: " << theSimulation->title << std::endl;
        }
        else if (!strcmp(keyword, "WORKINGDIR")) {
            std::string working_dir;
            input >> keyword;
            chdir(keyword);
            working_dir = getcwd(keyword, 100);
            output << "WORKINGDIR: " << working_dir << std::endl;
        }
        else if (!strcmp(keyword, "GRAVITY")) {
            double value;
            input >> value;
            Simulation::gravity(0) = value;
            input >> value;
            Simulation::gravity(1) = value;
            input >> value;
            Simulation::gravity(2) = value;
            output << "GRAVITY: "
            << Simulation::getGravity(0) << ", "
            << Simulation::getGravity(1) << ", "
            << Simulation::getGravity(2) << std::endl;
        }
        else if (!strcmp(keyword, "DIMENSION")) {
            int value;

            input >> value;
            Simulation::dimension = value;

            output << "DIMENSION: "
            << Simulation::dimension << std::endl;
        }
        else if (!strcmp(keyword, "MATERIALS")) {
            char a;
            do {
                input.get(a);
            }
            while (a != '\n');
            output << "MATERIALS: " << std::endl;
            while (input >> keyword) {
                if (!strcmp(keyword, "ENDMATERIALS")) {
                    break;
                    // TODO Add PLSTRESS and 3D
                } else if (!strcmp(keyword, "PLSTRAIN")) {
                    int num_mat;
                    double young, poisson, density;
                    input >> num_mat >> young >> poisson >> density;

                    if (theSimulation->materials.count(num_mat) == 0) {
                        theSimulation->materials[num_mat];
                    }

                    theSimulation->materials.at(num_mat).setMechanicalProps(Simulation::dimension, young, poisson,
                                                                            density);
                    theSimulation->materials.at(num_mat).setMaterialId(num_mat);
                    output << "MATERIAL: " << keyword
                    << ", number = " << num_mat << ",\n\t E = " << young
                    << ", mu = " << poisson << ", density = " << density << std::endl;
                    do {
                        input.get(a);
                    }
                    while (a != '\n');
                }
                else if (!strcmp(keyword, "THERMAL")) {
                    int num_mat;
                    double capacity, kappa, beta, density; // Capacity, Conductivity, Thermal expansion, Density
                    input >> num_mat >> capacity >> kappa >> beta >> density;

                    if (theSimulation->materials.count(num_mat) == 0) {
                        theSimulation->materials[num_mat];
                    }
// 		      theSimulation->materials.insert(std::pair<int,Material>( num_mat, Material() ) );

                    theSimulation->materials.at(num_mat).setThermalProps(capacity, kappa, beta, density);
                    theSimulation->materials.at(num_mat).setMaterialId(num_mat);
                    //std::cout << "Material " << num_mat << " set" << std::endl;

                    output << "MATERIAL: " << keyword
                    << ", number = " << num_mat << ",\n\t Cp = " << capacity
                    << ", kappa = " << kappa << ", density = " << density << std::endl;
                    do {
                        input.get(a);
                    }
                    while (a != '\n');
                }
                else if (!strcmp(keyword, "FILES")) {
                    while (input >> keyword) {
                        if (!strcmp(keyword, "ENDFILES")) {
                            break;
                        } else if (!strcmp(keyword, "CAPACITY")) {
                            int num_mat;
                            double temperature, capacity;
                            char a;
                            input >> num_mat >> keyword;
                            if (theSimulation->materials.count(num_mat) == 0) {
                                theSimulation->materials[num_mat];
                            }
                            output << "THERMALFILE CAPACITY: " << keyword
                            << " for mat # " << num_mat << endl;
                            std::ifstream thermalfile(keyword); // file to read points from
                            while (thermalfile >> temperature) {
                                thermalfile >> capacity;
                                theSimulation->materials.at(num_mat).addThermalCapacity(temperature, capacity);
                                output << "TEMP:" << temperature << ", \t" << capacity << endl;
                            }
                            do {
                                input.get(a);
                            }
                            while (a != '\n');
                        }
                        else if (!strcmp(keyword, "CONDUCTIVITY")) {
                            int num_mat;
                            double temperature, conductivity;
                            char a;
                            input >> num_mat >> keyword;
                            if (theSimulation->materials.count(num_mat) == 0) {
                                theSimulation->materials[num_mat];
                            }
                            output << "THERMALFILE CONDUCTIVITY: " << keyword
                            << " for mat # " << num_mat << endl;
                            std::ifstream thermalfile(keyword); // file to read points from
                            while (thermalfile >> temperature) {
                                thermalfile >> conductivity;
                                theSimulation->materials.at(num_mat).addThermalConductivity(temperature, conductivity);
                                output << "TEMP:" << temperature << ", \t" << conductivity << endl;
                            }
                            do {
                                input.get(a);
                            }
                            while (a != '\n');
                        }
                    }
                }
            }
            //
          theSimulation->myMaterialTable = (MaterialTable*)malloc(1*sizeof(MaterialTable));
          setupMaterialTables(&theSimulation->myMaterialTable,theSimulation->materials);
            //
        }
        else if (!strcmp(keyword, "CONTACT")) {
            std::string type;
            input >> type; // Valid types: GLOBAL, NONE
            Simulation::contact = type;
            output << "CONTACT: "
            << Simulation::contact << std::endl;
        }
        else if (!strcmp(keyword, "VISUALIZATION")) {
            std::string type;
            input >> type;
            if (type == "ON") Simulation::visualization = 1;
            if (type == "OFF")Simulation::visualization = 0;
            output << "VISUALIZATION: "
            << Simulation::contact << std::endl;
        }
        else if (!strcmp(keyword, "OUTPUT")) {
            std::string type;
            input >> type;
            if (type == "MATRICES") {
                Simulation::outputMatrices = 1;
                output << "OUTPUT: MATRICES"
                << Simulation::outputMatrices << std::endl;
            }
//      else if(type == "OTHER")
        }
        else if (!strcmp(keyword, "SMOOTHING")) {
            std::string type;
            input >> type;
            if (type == "OFF") {
                Simulation::smoothingType = "OFF";
                output << keyword << " "
                << Simulation::smoothingType << std::endl;
            }
            if (type == "LOCAL") {
                Simulation::smoothingType = "LOCAL";
                output << keyword << " "
                << Simulation::smoothingType << std::endl;
            }
            if (type == "CONSTANT") {
                Simulation::smoothingType = "CONSTANT";
                output << keyword << " "
                << Simulation::smoothingType << std::endl;
            }
            if (type == "GLOBAL") {
                Simulation::smoothingType = "GLOBAL";
                output << keyword << " "
                << Simulation::smoothingType << std::endl;
            }
//      else if(type == "OTHER")
        }
        else if (!strcmp(keyword, "INITIALTEMPERATURE")) {
            double init_temp;
            input >> init_temp; // Type of formulation: EFG or RPIM
            output << "SIMULATION INITIAL TEMPERATURE SET TO: "
            << init_temp
            << std::endl;
            theSimulation->setInitialTemperatures(init_temp);
        }
        else if (!strcmp(keyword, "SYSTEM")) {
            readSystem(theSimulation->baseSystem.get());
        }
        else if (!strcmp(keyword, "ANALYSIS")) {
            readAnalysis();
        }
    }

}

void mknix::Reader::readSystem(System * system_in)
{
    std::string sysTitle;
    char keyword[20];
    input >> sysTitle;
    output << "SYSTEM: "
    << system_in->getTitle()
    << "."
    << sysTitle << std::endl;

    system_in->subSystems[sysTitle] = new System(sysTitle);
    while (input >> keyword) {
        if (!strcmp(keyword, "ENDSYSTEM")) {
            return;

        } else if (!strcmp(keyword, "RIGIDBODIES")) {
            theReaderRigid = new ReaderRigid(theSimulation, output, input);
            theReaderRigid->readRigidBodies(system_in->subSystems[sysTitle]);
            delete theReaderRigid;
            theReaderRigid = nullptr;
        }
        else if (!strcmp(keyword, "FLEXBODIES")) {
            theReaderFlex = new ReaderFlex(theSimulation, output, input);
            theReaderFlex->readFlexBodies(system_in->subSystems[sysTitle]);
            delete theReaderFlex;
            theReaderFlex = nullptr;
        }
        else if (!strcmp(keyword, "BODYPOINTS")) {
            this->readBodyPoints(system_in->subSystems[sysTitle]);
        }
        else if (!strcmp(keyword, "JOINTS")) {
            theReaderConstraints
                    = new ReaderConstraints(theSimulation, output, input);
            theReaderConstraints->readConstraints(system_in->subSystems[sysTitle]);
            delete theReaderConstraints;
            theReaderConstraints = nullptr;
        }
        else if (!strcmp(keyword, "LOADS")) {
            this->readLoads(system_in->subSystems[sysTitle]);
        }
        else if (!strcmp(keyword, "ENVIRONMENT")) {
            this->readEnvironment(system_in->subSystems[sysTitle]);
        }
        else if (!strcmp(keyword, "MOTION")) {
            this->readMotion(system_in->subSystems[sysTitle]);
        }
        else if (!strcmp(keyword, "SCALE")) {
            double temp, xValue, yValue, zValue;

            input >> xValue >> yValue >> zValue;

            for (auto& node : theSimulation->nodes) {
                temp = node.second->getX();
                node.second->setX(temp * xValue);
                temp = node.second->getY();
                node.second->setY(temp * yValue);
                temp = node.second->getZ();
                node.second->setZ(temp * zValue);
            }
        }
        else if (!strcmp(keyword, "MIRROR")) {
            double temp;
            std::string axis;

            input >> axis;

            if (axis == "x" || axis == "X") {
                for (auto& node : theSimulation->nodes) {
                    temp = node.second->getX();
                    node.second->setX(-temp);
                }
            }
            else if (axis == "y" || axis == "Y") {
                for (auto& node : theSimulation->nodes) {
                    temp = node.second->getY();
                    node.second->setY(-temp);
                }
            }
            else if (axis == "z" || axis == "Z") {
                for (auto& node : theSimulation->nodes) {
                    temp = node.second->getZ();
                    node.second->setZ(-temp);
                }
            }
        }
        else if (!strcmp(keyword, "SHIFT")) {
            double temp, xValue, yValue, zValue;

            input >> xValue >> yValue >> zValue;

            for (auto& node : theSimulation->nodes) {
                temp = node.second->getX();
                node.second->setX(temp + xValue);
                temp = node.second->getY();
                node.second->setY(temp + yValue);
                temp = node.second->getZ();
                node.second->setZ(temp + zValue);
            }
        }
    }
}

void mknix::Reader::readBodyPoints(System * system_in)
{
    char keyword[20];
//   std::string loadTitle;

    Node * pNode = 0;
    double x, y, z;

    while (input >> keyword) {
        if (!strcmp(keyword, "ENDBODYPOINTS")) return;

//       if (a == '.'){ // we read the node...
//         while(input.get(a)){
//           if (a == '\n') break;
//           else if (a==' ') break;
//           else{
//             sNode.push_back( a );
//           }
//         }
//       }
        input >> x >> y >> z;

        int nodeNumber(0);

        if (system_in->rigidBodies.find(keyword) !=
            system_in->rigidBodies.end()) //if the body is a rigidbody
        {
            if (system_in->rigidBodies[keyword]->getNodesSize() > 0) {
                nodeNumber = system_in->rigidBodies[keyword]->getLastNode()->getNumber();
                ++nodeNumber;
            }

            pNode = new Node(nodeNumber, x, y, z);
            system_in->rigidBodies[keyword]->addNode(pNode);
            cout << "BODYPOINTS in " << keyword << ", Number " << nodeNumber << endl;
        }
        else { //the body is a flexbody
            double alpha, dc;
            input >> alpha >> dc;
            nodeNumber = system_in->flexBodies[keyword]->getNumberOfPoints();

            system_in->flexBodies[keyword]->addPoint(nodeNumber, x, y, z, alpha, dc);
        }
    }
}


void mknix::Reader::readLoads(System * system_in)
{
    char keyword[20];
//   std::string loadTitle;

    while (input >> keyword) {
        if (!strcmp(keyword, "ENDLOADS")) {
            return;

        } else if (!strcmp(keyword, "FORCE")) {
            std::string sBody;
            std::string sNode;
            Node * pNode = 0;
            char a;
            double fx, fy, fz;

            while (input.get(a)) { // we read the body...
                if (a == '.') {
                    break;
                } else if (a == '\n') {
                    break;
                } else if (a == ' ') { //get blank space
                } else {
                    sBody.push_back(a);
                }
            }
            if (a == '.') { // we read the node...
                while (input.get(a)) {
                    if (a == '\n') {
                        break;
                    } else if (a == ' ') {
                        break;
                    } else {
                        sNode.push_back(a);
                    }
                }
            }
            if (system_in->rigidBodies.find(sBody) !=
                system_in->rigidBodies.end()) { //if the body is a rigidbody
                    pNode = system_in->rigidBodies[sBody]
    // 		  ->getDomainNode( sNode );
                            ->getNode(atoi(sNode.c_str()));
                } else { //the body is a flexbody
                    pNode = system_in->flexBodies[sBody]->getNode(atoi(sNode.c_str()));
            }

            input >> fx >> fy >> fz;

            system_in->loads.push_back(new Force(pNode, fx, fy, fz));
        }

        else if (!strcmp(keyword, "THERMALFLUENCE")) {
            std::string sBody;
            std::string sNode;
            Node * pNode = 0;
            char a;
            double fluence;

            while (input.get(a)) { // we read the body...
                if (a == '.') {
                    break;
                } else if (a == '\n') {
                    break;
                } else if (a == ' ') { //get blank space
                } else {
                    sBody.push_back(a);
                }
            }
            if (a == '.') { // we read the node...
                while (input.get(a)) {
                    if (a == '\n') {
                        break;
                    } else if (a == ' ') {
                        break;
                    } else {
                        sNode.push_back(a);
                    }
                }
            }
            if (system_in->rigidBodies.find(sBody) !=
                system_in->rigidBodies.end()) { //if the body is a rigidbody
                    pNode = system_in->rigidBodies[sBody]
    // 		  ->getDomainNode( sNode );
                            ->getNode(atoi(sNode.c_str()));
                } else { //the body is a flexbody
                    pNode = system_in->flexBodies[sBody]->getNode(atoi(sNode.c_str()));
            }

            input >> fluence;
            output << "THERMALFLUENCE " << pNode->getNumber() << " " << fluence << endl;

            system_in->loadsThermal.push_back(new LoadThermal(pNode, fluence));
        }
        else if (!strcmp(keyword, "THERMALOUTPUT")) {
            std::string sBody;
            std::string sNode;
            Node * pNode = 0;
            char a;

            while (input.get(a)) { // we read the body...
                if (a == '.') {
                    break;
                } else if (a == '\n') {
                    break;
                } else if (a == ' ') { //get blank space
                } else {
                    sBody.push_back(a);
                }
            }
            if (a == '.') { // we read the node...
                while (input.get(a)) {
                    if (a == '\n') {
                        break;
                    } else if (a == ' ') {
                        break;
                    } else {
                        sNode.push_back(a);
                    }
                }
            }
            if (system_in->rigidBodies.find(sBody) !=
                system_in->rigidBodies.end()) { //if the body is a rigidbody
                    pNode = system_in->rigidBodies[sBody]
    // 		  ->getDomainNode( sNode );
                            ->getNode(atoi(sNode.c_str()));
                } else if (sBody == "MAX_INTERFACE_TEMP") { // avoid else
            }
            else { //the body is a flexbody
                pNode = system_in->flexBodies[sBody]->getNode(atoi(sNode.c_str()));
            }

            if (sBody == "MAX_INTERFACE_TEMP") {
                system_in->outputMaxInterfaceTemp = true;
            }
            else {
                output << "THERMALOUTPUT " << pNode->getNumber() << endl;

                system_in->outputSignalThermal.push_back(pNode);
            }
        }

        else if (!strcmp(keyword, "THERMALBODY")) {
            std::string sBody;
            char a;
            while (input.get(a)) { // we read the body...
                if (a == '.') {
                    break;
                } else if (a == '\n') {
                    break;
                } else if (a == ' ') { //get blank space
                } else {
                    sBody.push_back(a);
                }
            }
            system_in->thermalBodies[sBody]->setLoadThermal(new LoadThermalBody());
        }
        else if (!strcmp(keyword, "THERMALFLUX1D")) { //Improvement from above
            std::string sBody, sBoundary;
            char a;
            while (input.get(a)) { // we read the body...
                if (a == '.') {
                    while (input.get(a)) { // we read the boundary...
                        if (a == '.') {
                            break;
                        } else if (a == '\n') {
                            break;
                        } else if (a == ' ') { //get blank space
                        } else {
                            sBoundary.push_back(a);
                        }
                    }
                    break;
                }

                else if (a == '\n') {
                    break;
                } else if (a == ' ') { //get blank space
                } else {
                    sBody.push_back(a);
                }
                output << "THERMALFLUX1D " << sBody << "." << sBoundary << endl;

            }
            LoadThermalBoundary1D * temp = new LoadThermalBoundary1D();
            system_in->thermalBodies[sBody]->setLoadThermalInBoundaryGroup(temp, sBoundary);
            while (input >> keyword) { // Options for definition: 2D, 3D, rad map in file
                if (!strcmp(keyword, "ENDTHERMALFLUX1D")) {
                    return;
                } else if (!strcmp(keyword, "FILE")) {
                    input >> keyword;
                    temp->loadFile(keyword);
                }
                else if (!strcmp(keyword, "TIMEFILE")) {
                    input >> keyword;
                    temp->loadTimeFile(keyword);
                }
                else if (!strcmp(keyword, "SCALE")) {
                    double factor;
                    input >> factor;
                    temp->scaleLoad(factor);
                }
            }
        }
        else if (!strcmp(keyword, "RADIATION")) {
            int mapDim = 3, lineSkip(0);
            double lengthFactor(1.), doseFactor(1.);
            double scaleX(1.), scaleY(1.), scaleZ(1.);
            Radiation * theRad = new Radiation;

            system_in->loads.push_back(theRad);

            output << "ENVIRONMENT "
            << system_in->getTitle()
            << "."
            << "RADIATION" << std::endl;

            while (input >> keyword) { // Options for definition: 2D, 3D, rad map in file
                if (!strcmp(keyword, "ENDRADIATION")) {
                    return;
                } else if (!strcmp(keyword, "STATIC3D")) {
                    mapDim = 3;
                }
                else if (!strcmp(keyword, "STATIC2D")) {
                    mapDim = 2;
                }
                else if (!strcmp(keyword, "SKIPLINES")) {
                    input >> lineSkip;
                }
                else if (!strcmp(keyword, "LENGTHFACTOR")) {
                    input >> lengthFactor;
                }
                else if (!strcmp(keyword, "SCALEAXIS")) {
                    input >> scaleX >> scaleY >> scaleZ;
                }
                else if (!strcmp(keyword, "DOSEFACTOR")) {
                    input >> doseFactor;
                }
                else if (!strcmp(keyword, "MAPFILE")) {
                    input >> keyword;
                    output << "FILE: " << keyword << endl;
                    std::ifstream mapFile(keyword); // file to read data from
                    cout << "RADFILE " << keyword << " is open? " << mapFile.is_open() << endl;
                    double x, y, z, value;
                    std::string junk;
                    for (int i = 0; i < lineSkip; ++i) {
                        std::getline(mapFile, junk); // thrash lines
                        cout << "Thrasing header from RADFILE:\n \t" << junk << endl;
                    }
                    while (mapFile >> x) {
                        mapFile >> y;
                        if (mapDim == 3) {
                            mapFile >> z;
                        }
                        else {
                            z = 0;
                        }
                        mapFile >> value; // input in mm and muSv
                        theRad->addVoxel(x * lengthFactor * scaleX,
                                         y * lengthFactor * scaleY,
                                         z * lengthFactor * scaleZ,
                                         value * doseFactor);
                    }
                }
            }
        }
    }
}

void mknix::Reader::readEnvironment(System * system_in)
{
    // Prepared to read radiation and convection.
    // At this moment they are implemented as loads
    char keyword[20];

    while (input >> keyword) {
        if (!strcmp(keyword, "ENDENVIRONMENT")) return;
    }
}

void mknix::Reader::readMotion(System * system_in)
{
    char keyword[20];
    std::string groundNodeNumber;
    Node * pNode = 0;
    std::map<double, double> timex, timey, timez;
    double time, ux, uy, uz;

    input >> groundNodeNumber; // Ground node to move
    pNode = system_in->getNode(atoi(groundNodeNumber.c_str()));
    output << "MOTION "
    << system_in->getTitle()
    << "."
    << groundNodeNumber << std::endl;

    system_in->motions.push_back(new Motion(pNode));

    while (input >> keyword) {
        if (!strcmp(keyword, "ENDMOTION")) {
            system_in->motions.back()->setTimeUx(timex);
            system_in->motions.back()->setTimeUy(timey);
            system_in->motions.back()->setTimeUz(timez);
            return;
        }
        else if (!strcmp(keyword, "TIMECONF")) {
            input >> time >> ux >> uy >> uz; // time and movement
            timex[time] = ux;
            timey[time] = uy;
            timez[time] = uz;
            output << "\t"
            << /*groundNodeNumber <<*/ "TIMECONF: "
            << time << " = "
            << ux << ", " << uy << ", " << uz << ", "
            << endl;
        }
    }
}

void mknix::Reader::readAnalysis()
{
    char keyword[20];

    output << "ANALYSYS: " << std::endl;

    while (input >> keyword) {
        if (!strcmp(keyword, "ENDANALYSIS")) {
            return;
        } else if (!strcmp(keyword, "STATIC")) {
            double time;
            output << "\t" << keyword << ":"
            << std::endl;
            while (input >> keyword) {
                if (!strcmp(keyword, "ENDSTATIC")) {
                    break;
                } else if (!strcmp(keyword, "EPSILON")) {
                    input >> time; //just to not create another variable
                    Simulation::epsilon = time;
                    output << "\t\t"
                    << "EPSILON: "
                    << Simulation::epsilon
                    << endl;
                }
                else if (!strcmp(keyword, "TIME")) {
                    input >> time;
                    output << "\t\t"
                    << "TIME: " << time << std::endl;
                }
            }
            this->theSimulation->analyses.push_back
                    (make_unique<AnalysisStatic>(theSimulation, time));

        }
        else if (!strcmp(keyword, "THERMALSTATIC")) {
            double time;
            output << "\t" << keyword << ":"
            << std::endl;
            while (input >> keyword) {
                if (!strcmp(keyword, "ENDTHERMALSTATIC")) {
                    break;
                } else if (!strcmp(keyword, "EPSILON")) {
                    input >> time; //just to not create another variable
                    Simulation::epsilon = time;
                    output << "\t\t"
                    << "EPSILON: "
                    << Simulation::epsilon
                    << endl;
                }
                else if (!strcmp(keyword, "TIME")) {
                    input >> time;
                    output << "\t\t"
                    << "TIME: " << time << std::endl;
                }
            }
            this->theSimulation->analyses.push_back
                    (make_unique<AnalysisThermalStatic>(theSimulation, time));

        }
        else if (!strcmp(keyword, "THERMALDYNAMIC")) {
            char integratorType[20];
            double to, tf, At;
            output << "\t" << keyword << ":"
            << std::endl;
            while (input >> keyword) {
                if (!strcmp(keyword, "ENDTHERMALDYNAMIC")) {
                    break;
                } else if (!strcmp(keyword, "EPSILON")) {
                    input >> to; //just to not create another variable
                    Simulation::epsilon = to;
                    output << "\t\t"
                    << "EPSILON: "
                    << Simulation::epsilon
                    << endl;
                }
                else if (!strcmp(keyword, "INTEGRATOR")) {
                    input >> integratorType;
                    output << "\t\t"
                    << "INTEGRATOR: " << integratorType
                    << std::endl;
                }
                else if (!strcmp(keyword, "TIME")) {
                    input >> to
                    >> tf
                    >> At;
                    output << "\t\t"
                    << "TIME: "
                    << to << ", "
                    << tf << ", "
                    << At << std::endl;
                }
            }
            this->theSimulation->analyses.push_back
                    (make_unique<AnalysisThermalDynamic>(theSimulation, to, tf, At, integratorType));
        }
        else if (!strcmp(keyword, "THERMOMECHANICALDYNAMIC")) {
            char integratorType[20];
            double to, tf, At;
            output << "\t" << keyword << ":"
            << std::endl;
            while (input >> keyword) {
                if (!strcmp(keyword, "ENDTHERMOMECHANICALDYNAMIC")) {
                    break;
                } else if (!strcmp(keyword, "EPSILON")) {
                    input >> to; //just to not create another variable
                    Simulation::epsilon = to;
                    output << "\t\t"
                    << "EPSILON: "
                    << Simulation::epsilon
                    << endl;
                }
                else if (!strcmp(keyword, "INTEGRATOR")) {
                    input >> integratorType;
                    output << "\t\t"
                    << "INTEGRATOR: " << integratorType
                    << std::endl;
                }
                else if (!strcmp(keyword, "TIME")) {
                    input >> to
                    >> tf
                    >> At;
                    output << "\t\t"
                    << "TIME: "
                    << to << ", "
                    << tf << ", "
                    << At << std::endl;
                }
            }
            this->theSimulation->analyses.push_back
                    (make_unique<AnalysisThermoMechanicalDynamic>(theSimulation, to, tf, At, integratorType));
        }
        else if (!strcmp(keyword, "DYNAMIC")) {
            char integratorType[20];
            double to, tf, At;
            double par1 = -1.;
            double par2 = -1.;
            double par3 = -1.;
            output << "\t" << keyword << ":"
            << std::endl;
            while (input >> keyword) {
                if (!strcmp(keyword, "ENDDYNAMIC")) {
                    break;
                } else if (!strcmp(keyword, "EPSILON")) {
                    input >> to; //just to not create another variable
                    Simulation::epsilon = to;
                    output << "\t\t"
                    << "EPSILON: "
                    << Simulation::epsilon
                    << endl;
                }
                else if (!strcmp(keyword, "INTEGRATOR")) {
                    input >> integratorType;
                    if (!strcmp(integratorType, "NEWMARK")) {
                        input >> par1 >> par2;
                        output << "\t\t"
                        << "INTEGRATOR: "
                        << integratorType
                        << " " << par1 // beta
                        << " " << par2 // gamma
                        << std::endl;
                    }
                    else if (!strcmp(integratorType, "NEWMARK-ALPHA")) {
                        strcpy(integratorType, "NEWMARK");
                        input >> par1 >> par2;
                        output << "\t\t"
                        << "INTEGRATOR: "
                        << integratorType
                        << " " << par1
                        << " " << par2
                        << " " << par3
                        << std::endl;
                    }
                    else if (!strcmp(integratorType, "HHT-SIMPLE")) {
                        strcpy(integratorType, "ALPHA");
                        input >> par1;
                        output << "\t\t"
                        << "INTEGRATOR: " << integratorType << " " << par1
                        << std::endl;
                    }
                    else if (!strcmp(integratorType, "HHT-GENERALIZED")) {
                        strcpy(integratorType, "ALPHA");
                        input >> par1;
                        output << "\t\t"
                        << "INTEGRATOR: "
                        << integratorType
                        << " " << par1
                        << " " << par2
                        << " " << par3
                        << std::endl;
                    }
                    else {
                        output << "\t\t"
                        << "INTEGRATOR: " << integratorType
                        << std::endl;
                    }
                }
                else if (!strcmp(keyword, "TIME")) {
                    input >> to
                    >> tf
                    >> At;
                    output << "\t\t"
                    << "TIME: "
                    << to << ", "
                    << tf << ", "
                    << At << std::endl;
                }
            }
            this->theSimulation->analyses.push_back
                    (make_unique<AnalysisDynamic>(theSimulation, to, tf, At, integratorType, par1, par2, par3));
        }
        else if (!strcmp(keyword, "OTRO")) {
        }
    }
}
