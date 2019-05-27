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
#include "readerrigid.h"

#include <core/node.h>
#include <simulation/simulation.h>
#include <system/bodyrigid0D.h>
#include <system/bodyrigid1D.h>
#include <system/bodyrigid2D.h>
#include <system/bodyrigid3D.h>
#include <system/constraintdistance.h>
#include <system/system.h>
#include <system/systemchain.h>

namespace mknix
{

ReaderRigid::ReaderRigid()
    : theSimulation(nullptr)
    , output(nullptr)
    , input(nullptr)
{
}

ReaderRigid::ReaderRigid(Simulation * simulation_in,
                         std::ofstream& output_in,
                         std::ifstream& input_in)
    : theSimulation(simulation_in)
    , output(&output_in)
    , input(&input_in)
{
}


ReaderRigid::~ReaderRigid()
{
}


} // namespace mknix

void mknix::ReaderRigid::readRigidBodies(System * system_in)
{
    std::string keyword;

    while (*input >> keyword)
    {
        if (keyword == "ENDRIGIDBODIES")
        {
            return;
        }
        /*Set Penalty method*/
        else if (keyword == "PENALTY")
        {
            Simulation::constraintMethod = "PENALTY";
            *output << "PENALTY set" << endl;
        }
        /*Set Augmented Lagrange method*/
        else if (keyword == "AUGMENTED")
        {
            Simulation::constraintMethod = "AUGMENTED";
            *output << "AUGMENTED set" << endl;
        }
        /*Set Alpha parameter*/
        else if (keyword == "ALPHA")
        {
            *input >> Simulation::alpha;
            *output << "ALPHA: " << Simulation::getAlpha() << endl;
        }
        /*Defined RigidBody:0D:MassPoint*/
        else if (keyword == "MASSPOINT")
        {
            this->readRigidBody0DMassPoint(system_in);
        }
        /*Defined RigidBody:1D:Bar*/
        else if (keyword == "BAR")
        {
            this->readRigidBody1DBar(system_in);
        }
        else if (keyword == "CHAIN")
        {
            this->readRigidBody1DChain(system_in);
        }
        /*Defined RigidBody:2D:Generic*/
        else if (keyword == "GENERIC2D")
        {
            this->readRigidBody2DMesh(system_in);
        }
        /*Defined RigidBody:3D:Generic*/
        else if (keyword == "GENERIC3D")
        {
            this->readRigidBody3DGeneric(system_in);
        }
        else if (keyword == "OTHER")
        {
        }
    }
}

void mknix::ReaderRigid::readRigidBody0DMassPoint(System * system_in)
{
    std::string keyword;

    std::string rigidTitle;
    std::string nameNode;

    double xA, yA, zA;
    double mass = 0.0;

    /*Read Rigid Body Title*/
    *input >> rigidTitle;
    *output << "MASSPOINT: "
            << system_in->getTitle()
            << "."
            << rigidTitle << std::endl;

    while (*input >> keyword)
    {
        if (keyword == "ENDMASSPOINT")
        {
            break;
        }
        /*Read NodeA*/
        else if (keyword == "NODEA")
        {
            nameNode = rigidTitle + ".NODEA";
            this->readNode(xA, yA, zA, nameNode);
        }
        /*Read Mass*/
        else if (keyword == "MASS")
        {
            *input >> mass;
            *output << "\t"
                    << rigidTitle << ".MASS: " << mass << endl;
        }
    }

    /*New Object: MassPoint*/
    int nodeA;
    nodeA = theSimulation->nodes.end()->first;

    theSimulation->nodes[nodeA] = new Node(nodeA, xA, yA, zA);

    system_in->rigidBodies[rigidTitle] = new RigidBodyMassPoint(rigidTitle,
            theSimulation->nodes[nodeA],
            mass);
}

void mknix::ReaderRigid::readRigidBody1DBar(System * system_in)
{
    std::string keyword;

    std::string rigidTitle;
    std::string nameNode;
    std::string energyKeyword; //chapuza

    double xA, yA, zA;
    double xB, yB, zB;
    double mass;

    /*Read Rigid Body Title*/
    *input >> rigidTitle;
    *output << "BAR: "
            << system_in->getTitle()
            << "."
            << rigidTitle << std::endl;

    while (*input >> keyword)
    {
        if (keyword == "ENDBAR")
        {
            break;
        }
        else if (keyword == "OUTPUT")
        {
            *input >> energyKeyword;
            if (energyKeyword == "ENERGY")
            {
                *output << "OUTPUT: " << energyKeyword << endl;
            }
        }
        /*Read NodeA*/
        else if (keyword == "NODEA")
        {
            nameNode = rigidTitle + ".NODEA";
            this->readNode(xA, yA, zA, nameNode);
        }
        /*Read NodeB*/
        else if (keyword == "NODEB")
        {
            nameNode = rigidTitle + ".NODEB";
            this->readNode(xB, yB, zB, nameNode);
        }
        /*Read Mass*/
        else if (keyword == "MASS")
        {
            *input >> mass;
            *output << "\t"
                    << rigidTitle << ".MASS: " << mass << endl;
        }
    }

    int lastNode;
    lastNode = theSimulation->nodes.end()->first;
    cout << "LASTNODE: " << lastNode << endl;

    theSimulation->nodes[lastNode] = new Node(lastNode, xA, yA, zA);
    theSimulation->nodes[lastNode + 1] = new Node(lastNode + 1, xB, yB, zB);
    theSimulation->outputPoints[lastNode] = theSimulation->nodes[lastNode];
    theSimulation->outputPoints[lastNode + 1] = theSimulation->nodes[lastNode + 1];

    system_in->rigidBodies[rigidTitle] = new RigidBar(rigidTitle,
            theSimulation->nodes[lastNode + 0],
            theSimulation->nodes[lastNode + 1],
            mass
                                                     );

    system_in->constraints[rigidTitle] = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 0],
        theSimulation->nodes[lastNode + 1],
        Simulation::alpha,
        Simulation::constraintMethod
    );

    system_in->rigidBodies[rigidTitle]->setOutput(energyKeyword); //chapuza

}

void mknix::ReaderRigid::readRigidBody1DChain(System * system_in)
{
    std::string keyword;

    std::string rigidTitle;
    std::string nameNode;
    std::string energyKeyword; //chapuza

    double xA, yA, zA;
    double xB, yB, zB;
    double mass, length, temp1, temp2;
    std::map<double, double> timelengths;
    int segments;

    /*Read Rigid Body Title*/
    *input >> rigidTitle;
    *output << "CHAIN: "
            << system_in->getTitle()
            << "."
            << rigidTitle << std::endl;

    while (*input >> keyword)
    {
        if (keyword == "ENDCHAIN")
        {
            break;
        }
        else if (keyword == "OUTPUT")
        {
            *input >> energyKeyword;
            if (energyKeyword == "ENERGY")
            {
                *output << "OUTPUT: " << energyKeyword << endl;
            }
        }
        /*Read NodeA*/
        else if (keyword == "NODEA")
        {
            nameNode = rigidTitle + ".NODEA";
            this->readNode(xA, yA, zA, nameNode);
        }
        /*Read NodeB*/
        else if (keyword == "NODEB")
        {
            nameNode = rigidTitle + ".NODEB";
            this->readNode(xB, yB, zB, nameNode);
        }
        /*Read Mass*/
        else if (keyword == "MASS")
        {
            *input >> mass;
            *output << "\t"
                    << rigidTitle << ".MASS: " << mass << endl;
        }
        else if (keyword == "SEGMENTS")
        {
            *input >> segments;
            *output << "\t"
                    << rigidTitle << ".SEGMENTS: " << segments << endl;
        }
        else if (keyword == "LENGTH")
        {
            *input >> length;
            *output << "\t"
                    << rigidTitle << ".LENGTH: " << length << endl;
        }
        else if (keyword == "TIMELENGTH")
        {
            *input >> temp1 >> temp2; // time and length
            timelengths[temp1] = temp2;
            *output << "\t"
                    << rigidTitle << ".TIMELENGTH: "
                    << temp1 << ", "
                    << timelengths[temp1] << endl;
        }
    }

    // Create system in system_in
    SystemChain * theChain = new SystemChain(rigidTitle.c_str());
    system_in->subSystems[rigidTitle] = theChain;
    // Define the interface nodes
    theChain->setInterfaceNodeA(xA, yA, zA);
    theChain->setInterfaceNodeB(xB, yB, zB);
    // Define total mass
    theChain->setMass(mass);
    // Tell the system to populate nodes telling the number of the lastNode in theSimulation
    // and the number of segments
    theChain->setProperties(segments, length);
    theChain->populate(theSimulation, energyKeyword);
    // Define the temporal behaviour
    theChain->setTimeLengths(timelengths);
}

void mknix::ReaderRigid::readRigidBody2DMesh(System * system_in)
{
    std::string rigidTitle;
    std::string keyword;

    std::string stringKeyword;
    std::string meshfreeFormulation("RPIM");
    std::string nameNode;

    double aValue;
    std::vector<double> someValues;

    /*Read Rigid Body Title*/
    *input >> rigidTitle;
    *output << "GENERIC2D: "
            << system_in->getTitle()
            << "."
            << rigidTitle << std::endl;

    int lastNode(0);
    if (theSimulation->nodes.size() != 0)
    {
        lastNode = theSimulation->nodes.rbegin()->first + 1;
    }
//     lastNode = theSimulation->nodes.end()->first ;
    theSimulation->nodes[lastNode] = new Node(lastNode, 0., 0., 1.);
    theSimulation->nodes[lastNode + 1] = new Node(lastNode + 1, 1., 0., 1.);
    theSimulation->nodes[lastNode + 2] = new Node(lastNode + 2, 0., 1., 1.);

    theSimulation->outputPoints[lastNode] = theSimulation->nodes[lastNode];
    theSimulation->outputPoints[lastNode + 1] = theSimulation->nodes[lastNode + 1];
    theSimulation->outputPoints[lastNode + 2] = theSimulation->nodes[lastNode + 2];

    system_in->rigidBodies[rigidTitle] = new RigidBody2D(rigidTitle,
            theSimulation->nodes[lastNode + 0],
            theSimulation->nodes[lastNode + 1],
            theSimulation->nodes[lastNode + 2]
                                                        );

    system_in->constraints[rigidTitle + std::string("A")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 0],
        theSimulation->nodes[lastNode + 1],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("B")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 0],
        theSimulation->nodes[lastNode + 2],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("C")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 1],
        theSimulation->nodes[lastNode + 2],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    /*Read Rigid Body Attributes*/
    while (*input >> keyword)
    {
        if (keyword == "ENDGENERIC2D")
        {
            break;
        }
        else if (keyword == "DENSITY")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".DENSITYFACTOR: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setDensityFactor(aValue);
        }
        else if (keyword == "OUTPUT")
        {
            *input >> stringKeyword;
            if (stringKeyword == "ENERGY")
            {
                system_in->rigidBodies[rigidTitle]->setOutput(stringKeyword);
                *output << "OUTPUT: " << stringKeyword << endl;
            }
        }
        else if (keyword == "MASS")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".MASS: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setMass(aValue);
        }
        else if (keyword == "IXX")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IXX: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 0); // 0 is xx-axis
        }
        else if (keyword == "IYY")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IYY: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 1); // 1 is yy-axis
        }
        else if (keyword == "IXY")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IXY: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 2); // 2 is xy-axis
        }
        else if (keyword == "POSITION")
        {
            someValues.clear();
            // Read three values: CoG (x, y) and rotation angle (Theta)
            *input >> aValue;
            someValues.push_back(aValue);
            *output << "\t"
                    << rigidTitle << ".POSITION: CoG (" << aValue;
            *input >> aValue;
            someValues.push_back(aValue);
            *output << ", " << aValue;
            *input >> aValue;
            someValues.push_back(aValue);
            *output << "), Rotation (" << aValue << ")" << endl;
            system_in->rigidBodies[rigidTitle]->setPosition(someValues);
        }
        else if (keyword == "TRIANGLES")
        {
            char a;
            std::string boundaryType("CLOCKWISE");
            int totalNodes, totalCells;
            *input >> keyword;
            *output << "FILE: " << keyword << endl;
            std::ifstream meshfile(keyword); // file to read points from
            int number; // read number and use it. Differs from FlexBody
            // Instead of using global numbering, the node is created inside of the body
            // Global number is unknown and is not part of the system formulation
            double x, y, z;
//       int old_size = theSimulation->nodes.end()->first;
//       std::map<int,int> nodesFile;
            meshfile >> totalNodes >> totalCells;
            *output << '\t' << "Nodes: " << totalNodes << endl;
            *output << "\tNode\tx\ty\tz \n";
            for (int i = 0; i < totalNodes; ++i)
            {
                meshfile >> number >> x >> y >> z;
                --number;
                system_in->rigidBodies[rigidTitle]
                ->addNode(new Node(number, x, y, z));
                *output
                        << '\t' << system_in->rigidBodies[rigidTitle]->getNode(number)->getNumber()
                        << "\t" << system_in->rigidBodies[rigidTitle]->getNode(number)->getX()
                        << "\t" << system_in->rigidBodies[rigidTitle]->getNode(number)->getY()
                        << "\t" << system_in->rigidBodies[rigidTitle]->getNode(number)->getZ()
                        << std::endl;
            }

            // Reading connectivity of the mesh
            int elementType, node1, node2, node3;
            std::vector<int> nodes(3);
            *output << '\t' << "Cells: " << totalCells << endl;
            *output << "\tCell\tnode1\tnode2\tnode3 \n";
//             old_size = theSimulation->cells.end()->first;
            for (int i = 0; i < totalCells; ++i)
            {
                meshfile >> number >> elementType;
                if (elementType == 203)
                {
                    meshfile >> node1 >> node2 >> node3;
                    nodes[0] = --node1;
                    nodes[1] = --node2;
                    nodes[2] = --node3;
                    system_in->rigidBodies[rigidTitle]->addBoundaryConnectivity(nodes);
                    *output << "\t" << i << "\t"
                            << nodes[0] << "\t"
                            << nodes[1] << "\t"
                            << nodes[2] << endl;
                }
                else if (elementType == 102)
                {
                    meshfile >> node1 >> node2;
                    if (boundaryType == "CLOCKWISE")
                    {
                        system_in->rigidBodies[rigidTitle]->addBoundaryLine
                        (system_in->rigidBodies[rigidTitle]->getNode(node1 - 1),
                         system_in->rigidBodies[rigidTitle]->getNode(node2 - 1)
                        );
                    }
                    else if (boundaryType == "COUNTERCLOCKWISE")
                    {
                        system_in->rigidBodies[rigidTitle]->addBoundaryLine
                        (system_in->rigidBodies[rigidTitle]->getNode(node2 - 1),
                         system_in->rigidBodies[rigidTitle]->getNode(node1 - 1)
                        );
                    }
                    else
                    {
                        cerr << ":::ERROR: BOUNDARY ORIENTATION UNKNOWN:::" << endl;
                    }
                    *output << "\t Bound"
                            << "\t" << node1
                            << "\t" << node2
                            << std::endl;
                }
                else
                {
                    do
                    {
                        meshfile.get(a);
                    }
                    while (a != '\n');
                }
            }
            do
            {
                input->get(a);
            }
            while (a != '\n');
        }

    }
}


void mknix::ReaderRigid::readRigidBody3DGeneric(System * system_in)
{
    std::string rigidTitle;
    std::string keyword;

    std::string stringKeyword;
    std::string meshfreeFormulation("RPIM");
    std::string nameNode;

    double aValue;
    std::vector<double> someValues;

    /*Read Rigid Body Title*/
    *input >> rigidTitle;
    *output << "GENERIC3D: "
            << system_in->getTitle()
            << "."
            << rigidTitle << std::endl;

    int lastNode(0);
    if (theSimulation->nodes.size() != 0)
    {
        lastNode = theSimulation->nodes.rbegin()->first + 1;
    }
//     lastNode = theSimulation->nodes.end()->first ;
    cout << "lastNode = " << lastNode << endl;

    theSimulation->nodes[lastNode] = new Node(lastNode, 0., 0., 0.);
    theSimulation->nodes[lastNode + 1] = new Node(lastNode + 1, 1., 0., 0.);
    theSimulation->nodes[lastNode + 2] = new Node(lastNode + 2, 0., 1., 0.);
    theSimulation->nodes[lastNode + 3] = new Node(lastNode + 3, 0., 0., 1.);

    theSimulation->outputPoints[lastNode] = theSimulation->nodes[lastNode];
    theSimulation->outputPoints[lastNode + 1] = theSimulation->nodes[lastNode + 1];
    theSimulation->outputPoints[lastNode + 2] = theSimulation->nodes[lastNode + 2];
    theSimulation->outputPoints[lastNode + 3] = theSimulation->nodes[lastNode + 3];

    system_in->rigidBodies[rigidTitle] = new RigidBody3D(rigidTitle,
            theSimulation->nodes[lastNode + 0],
            theSimulation->nodes[lastNode + 1],
            theSimulation->nodes[lastNode + 2],
            theSimulation->nodes[lastNode + 3]
                                                        );

    system_in->constraints[rigidTitle + std::string("01")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 0],
        theSimulation->nodes[lastNode + 1],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("02")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 0],
        theSimulation->nodes[lastNode + 2],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("03")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 0],
        theSimulation->nodes[lastNode + 3],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("12")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 1],
        theSimulation->nodes[lastNode + 2],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("23")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 2],
        theSimulation->nodes[lastNode + 3],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    system_in->constraints[rigidTitle + std::string("13")]
        = new ConstraintDistance( // Can be changed to a constraint in state space
        theSimulation->nodes[lastNode + 1],
        theSimulation->nodes[lastNode + 3],
        Simulation::alpha,
        Simulation::constraintMethod
    );
    /*Read Rigid Body Attributes*/
    while (*input >> keyword)
    {
        if (keyword == "ENDGENERIC3D")
        {
            break;
        }
        else if (keyword == "DENSITYFACTOR")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".DENSITYFACTOR: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setDensityFactor(aValue);
        }
        else if (keyword == "OUTPUT")
        {
            *input >> stringKeyword;
            if (stringKeyword == "ENERGY")
            {
                system_in->rigidBodies[rigidTitle]->setOutput(stringKeyword);
                *output << "OUTPUT: " << stringKeyword << endl;
            }
        }
        else if (keyword == "MASS")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".MASS: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setMass(aValue);
        }
        else if (keyword == "IXX")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IXX: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 0); // 0 is xx-axis
        }
        else if (keyword == "IYY")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IYY: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 1); // 1 is yy-axis
        }
        else if (keyword == "IZZ")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IZZ: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 2); // 2 is zz-axis
        }
        else if (keyword == "IXY")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IXY: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 3); // 3 is xy-axis
        }
        else if (keyword == "IYZ")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IYZ: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 4); // 4 is yz-axis
        }
        else if (keyword == "IXZ")
        {
            *input >> aValue;
            *output << "\t"
                    << rigidTitle << ".IXZ: " << aValue << endl;
            system_in->rigidBodies[rigidTitle]->setInertia(aValue, 5); // 5 is yz-axis
        }
        else if (keyword == "POSITION")
        {
            someValues.clear();
            // Read three values: CoG (x, y, z). No rotation angles yet
            *input >> aValue;
            someValues.push_back(aValue);
            *output << "\t"
                    << rigidTitle << ".POSITION: CoG (" << aValue;
            *input >> aValue;
            someValues.push_back(aValue);
            *output << ", " << aValue;
            *input >> aValue;
            someValues.push_back(aValue);
            *output << ", " << aValue << ")" << endl;
            system_in->rigidBodies[rigidTitle]->setPosition(someValues);
        }
        else if (keyword == "TETRAHEDRONS")
        {
            int totalNodes, totalCells;
            *input >> keyword;
            *output << "FILE: " << keyword << endl;
            std::ifstream meshfile(keyword); // file to read points from
            int number; // read number and use it. Differs from FlexBody
            // Instead of using global numbering, the node is created inside of the body
            // Global number is unknown and is not part of the system formulation
            double x, y, z;
            meshfile >> totalNodes >> totalCells;
            *output << '\t' << "Nodes: " << totalNodes << endl;
            *output << "\tNode\tx\ty\tz \n";
            for (int i = 0; i < totalNodes; ++i)
            {
                meshfile >> number >> x >> y >> z;
                --number;
                system_in->rigidBodies[rigidTitle]
                ->addNode(new Node(number, x, y, z));
                *output
                        << '\t' << system_in->rigidBodies[rigidTitle]->getNode(number)->getNumber()
                        << "\t" << system_in->rigidBodies[rigidTitle]->getNode(number)->getX()
                        << "\t" << system_in->rigidBodies[rigidTitle]->getNode(number)->getY()
                        << "\t" << system_in->rigidBodies[rigidTitle]->getNode(number)->getZ()
                        << std::endl;
            }
            int elementType, node1, node2, node3, node4/*, cellCount*/;
            std::vector<int> nodes(3);
            *output << '\t' << "Cells: " << totalCells << endl;
            *output << "\tCell\tnode1\tnode2\tnode3 \n";
//             old_size = theSimulation->cells.end()->first;
            for (int i = 0; i < totalCells; ++i)
            {
                meshfile >> number >> elementType;
//         cout << elementType << endl;
                if (elementType == 203)
                {
                    meshfile >> node1 >> node2 >> node3;
                    nodes[0] = --node1;
                    nodes[1] = --node2;
                    nodes[2] = --node3;
                    system_in->rigidBodies[rigidTitle]->addBoundaryConnectivity(nodes);
                    *output << "\t" << i << "\t"
                            << nodes[0] << "\t"
                            << nodes[1] << "\t"
                            << nodes[2] << endl;
                }
                else if (elementType == 102)   //TODO: read linear boundary
                {
                    // at the moment we dump the information
                    meshfile >> node1 >> node2;
                }
                else if (elementType == 304)   //TODO: read linear boundary
                {
                    // at the moment we dump the information
                    meshfile >> node1 >> node2 >> node3 >> node4;
                }
                else   //elementType is now the first node of the tetrahedral
                {
                    meshfile >> node1 >> node2; //these are the 3, & 4th
                    // TODO: Read the mesh and create cells for thermal formulation
//           if( !strcmp(bodyType,"MESHFREE") ){
//             system_in->flexBodies[flexTitle]
//               ->addCell(/*old_size +*/ cellCount,
//                         new CellTetrahedron
//                         ( theSimulation->materials[material],
//                           meshfreeFormulation,
//                           alpha,
//                           nGPs,
//                           theSimulation->nodes[nodesFile[nodeNumberInFile]],
//                           theSimulation->nodes[nodesFile[elementType]],
//                           theSimulation->nodes[nodesFile[node1]],
//                           theSimulation->nodes[nodesFile[node2]]
//                         )
//                         );
//           }
//           else if( !strcmp(bodyType,"FEMESH") ){
//             system_in->flexBodies[flexTitle]
//               ->addCell( /*old_size +*/ cellCount,
//                           new ElemTetrahedron
//                           ( theSimulation->materials[material],
//                             alpha,
//                             nGPs,
//                             theSimulation->nodes[nodesFile[nodeNumberInFile]],
//                             theSimulation->nodes[nodesFile[elementType]],
//                             theSimulation->nodes[nodesFile[node1]],
//                             theSimulation->nodes[nodesFile[node2]]
//                           )
//                         );
//           }
//           *output
//               << '\t'
//               << /*old_size +*/ cellCount
//               << "\t"
//               << theSimulation->nodes[nodesFile[nodeNumberInFile]]->getNumber()
//               << "\t"
//               << theSimulation->nodes[nodesFile[elementType]]->getNumber()
//               << "\t"
//               << theSimulation->nodes[nodesFile[node1]]->getNumber()
//               << "\t"
//               << theSimulation->nodes[nodesFile[node2]]->getNumber()
//               << std::endl;
//           ++cellCount;
                }
            }
        }
    }
}


void mknix::ReaderRigid::readNode(double& x_in,
                                  double& y_in,
                                  double& z_in,
                                  std::string name_in)
{
    *input >> x_in >> y_in >> z_in;
    *output << "\t"
            << name_in << " : "
            << x_in << ", "
            << y_in << ", "
            << z_in << std::endl;
}
