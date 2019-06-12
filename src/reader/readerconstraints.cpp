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
#include "readerconstraints.h"

#include <core/node.h>
#include <simulation/simulation.h>
#include <system/bodyflex.h>
#include <system/bodyrigid.h>
#include <system/constraintclearance.h>
#include <system/constraintdistance.h>
#include <system/constraintfixedaxis.h>
#include <system/constraintfixedcoordinates.h>
#include <system/constraintthermalfixed.h>
#include <system/system.h>

mknix::ReaderConstraints::ReaderConstraints()
    : theSimulation(0)
    , output(0)
    , input(0)
    , p_nodeA(0)
    , p_nodeB(0)
{
}

mknix::ReaderConstraints::ReaderConstraints(Simulation* simulation_in,
        std::ofstream& output_in,
        std::ifstream& input_in)
    : theSimulation(simulation_in)
    , output(&output_in)
    , input(&input_in)
    , p_nodeA(0)
    , p_nodeB(0)
{
}

mknix::ReaderConstraints::~ReaderConstraints()
{
}

void mknix::ReaderConstraints::readConstraints(System* system_in)
{
    std::string keyword;
    std::string consTitle;

    while (*input >> keyword)
    {
        if (keyword == "ENDJOINTS")
        {
            return;
        }
        else if (keyword == "PENALTY")
        {
            Simulation::constraintMethod = "PENALTY";
            *output << "PENALTY set" << endl;
        }
        else if (keyword == "AUGMENTED")
        {
            Simulation::constraintMethod = "AUGMENTED";
            *output << "AUGMENTED set" << endl;
        }
        else if (keyword == "ALPHA")
        {
            *input >> Simulation::alpha;
            *output << "ALPHA: "
                    << Simulation::getAlpha()
                    << endl;
        }
        else if (keyword == "SPHERICAL")
        {
            /* Igual a una restriccion de distancia constante */
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA = 0;
            p_nodeB = 0;

            *input >> consTitle;
            *output << "SPHERICAL: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;

            while (*input >> keyword)
            {
                if (keyword == "ENDSPHERICAL")
                {
                    break;
                }
                else if (keyword == "NODEA")
                {
                    this->readNodeName(bodyTitleA, nodeA);
                }
                else if (keyword == "NODEB")
                {
                    this->readNodeName(bodyTitleB, nodeB);
                }
            }
            this->assignConstraintNodes(system_in, consTitle, bodyTitleA, nodeA, bodyTitleB, nodeB);
            system_in->constraints[consTitle]
                = new ConstraintFixedCoordinates(p_nodeA, p_nodeB, Simulation::alpha, Simulation::constraintMethod);
            system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode(system_in, consTitle, "NODEA", bodyTitleA, nodeA, 0);
            this->outputConstraintNode(system_in, consTitle, "NODEB", bodyTitleB, nodeB, 1);
        }

        else if (keyword == "DISTANCE")
        {
            /* Igual a una restriccion de distancia constante */
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA = 0;
            p_nodeB = 0;

            *input >> consTitle;
            *output << "DISTANCE: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;

            while (*input >> keyword)
            {
                if (keyword == "ENDDISTANCE")
                {
                    break;
                }
                else if (keyword == "NODEA")
                {
                    this->readNodeName(bodyTitleA, nodeA);
                }
                else if (keyword == "NODEB")
                {
                    this->readNodeName(bodyTitleB, nodeB);
                }
            }
            this->assignConstraintNodes(system_in, consTitle, bodyTitleA, nodeA, bodyTitleB, nodeB);
            system_in->constraints[consTitle]
                = new ConstraintDistance(p_nodeA, p_nodeB, Simulation::alpha, Simulation::constraintMethod);
            system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode(system_in, consTitle, "NODEA", bodyTitleA, nodeA, 0);
            this->outputConstraintNode(system_in, consTitle, "NODEB", bodyTitleB, nodeB, 1);
        }

        else if (keyword == "AXIS")
        {
            // Igual a una restricción de distancia constante
            std::string axisName;
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA = 0;
            p_nodeB = 0;

            *input >> consTitle;
            *output << "AXIS: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;
            while (*input >> keyword)
            {
                if (keyword == "ENDAXIS")
                {
                    break;
                }
                else if (keyword == "DIRECTION")
                {
                    *input >> axisName;
                    *output << "DIRECTION: " << axisName << "-axis" << endl;
                }
                else if (keyword == "NODEA")
                {
                    this->readNodeName(bodyTitleA, nodeA);
                }
                else if (keyword == "NODEB")
                {
                    this->readNodeName(bodyTitleB, nodeB);
                }
            }
            this->assignConstraintNodes(system_in, consTitle, bodyTitleA, nodeA, bodyTitleB, nodeB);
            system_in->constraints[consTitle]
                = new ConstraintFixedAxis(p_nodeA, p_nodeB, axisName, Simulation::alpha,
                                          Simulation::constraintMethod);
            system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode(system_in, consTitle, "NODEA", bodyTitleA, nodeA, 0);
            this->outputConstraintNode(system_in, consTitle, "NODEB", bodyTitleB, nodeB, 1);
        }
        else if (keyword == "CLEARANCE")
        {
            // Igual a una restricción de distancia constante
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;
            double tolerance = 0.;

            p_nodeA = 0;
            p_nodeB = 0;

            *input >> consTitle;
            *output << "CLEARANCE: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;
            while (*input >> keyword)
            {
                if (keyword == "ENDCLEARANCE")
                {
                    break;
                }
                else if (keyword == "TOLERANCE")
                {
                    *input >> tolerance;
                    *output << "TOLERANCE: " << tolerance << endl;
                }
                else if (keyword == "NODEA")
                {
                    this->readNodeName(bodyTitleA, nodeA);
                }
                else if (keyword == "NODEB")
                {
                    this->readNodeName(bodyTitleB, nodeB);
                }
            }
            this->assignConstraintNodes(system_in, consTitle, bodyTitleA, nodeA, bodyTitleB, nodeB);
            system_in->constraints[consTitle]
                = new ConstraintClearance(p_nodeA, p_nodeB, tolerance, Simulation::alpha,
                                          Simulation::constraintMethod);
            system_in->constraints[consTitle]->setTitle(consTitle);

            this->outputConstraintNode(system_in, consTitle, "NODEA", bodyTitleA, nodeA, 0);
            this->outputConstraintNode(system_in, consTitle, "NODEB", bodyTitleB, nodeB, 1);
        }
        else if (keyword == "THERMALSPHERICAL")
        {
            /* Igual a una restriccion de distancia constante */
            std::string bodyTitleA, bodyTitleB;
            std::string nodeA, nodeB;

            p_nodeA = 0;
            p_nodeB = 0;

            *input >> consTitle;
            *output << "THERMALSPHERICAL: "
                    << system_in->getTitle()
                    << "."
                    << consTitle << std::endl;

            while (*input >> keyword)
            {
                if (keyword == "ENDTHERMALSPHERICAL")
                {
                    break;
                }
                else if (keyword == "NODEA")
                {
                    this->readNodeName(bodyTitleA, nodeA);
                }
                else if (keyword == "NODEB")
                {
                    this->readNodeName(bodyTitleB, nodeB);
                }
            }
            this->assignConstraintNodes(system_in, consTitle, bodyTitleA, nodeA, bodyTitleB, nodeB);
            system_in->constraintsThermal[consTitle]
                = new ConstraintThermalFixed(p_nodeA, p_nodeB, Simulation::alpha, Simulation::constraintMethod);
            system_in->constraintsThermal[consTitle]->setTitle(consTitle);

            *output << "THERMALSPHERICAL: "
                    << system_in->getTitle()
                    << "."
                    << system_in->constraintsThermal[consTitle]->getTitle() << std::endl;
            this->outputConstraintThermalNode(system_in, consTitle, "NODEA", bodyTitleA, nodeA, 0);
            this->outputConstraintThermalNode(system_in, consTitle, "NODEB", bodyTitleB, nodeB, 1);
        }

        else if (keyword == "OTRO")
        {
        }
    }
}


void mknix::ReaderConstraints::readNodeName(std::string& bodyTitle, std::string& node)
{
    char a;

    input->get(a); // get blank space
    while (input->get(a))
    {
        if (a == '.')
        {
            break;
        }
        else if (a == '\n')
        {
            break;
        }
        else
        {
            bodyTitle.push_back(a);
        }
    }
    if (a != '\n')
    {
        while (input->get(a))
        {
            if (a == '.')
            {
                break;
            }
            else if (a == '\n')
            {
                break;
            }
            else
            {
                node.push_back(a);;
            }
        }
    }
    cout << "NODE read: " << bodyTitle << "." << node << endl;
}


void mknix::ReaderConstraints::assignConstraintNodes(System* system_in,
        const std::string& consName,
        const std::string& bodyTitleA,
        const std::string& nodeA,
        const std::string& bodyTitleB,
        const std::string& nodeB)
{
    if (bodyTitleA == "GROUND")
    {
        Node* node = nullptr;

        if (system_in->subSystems.find(bodyTitleB) != system_in->subSystems.end())
        {
            /* If the body is a system */
            node = new Node(*system_in->subSystems[bodyTitleB]->getNode(atoi(nodeB.c_str())));
        }
        else if (system_in->rigidBodies.find(bodyTitleB) != system_in->rigidBodies.end())
        {
            /* If the body is a rigidbody */
            node = new Node(*system_in->rigidBodies[bodyTitleB]->getNode(atoi(nodeB.c_str())));
        }
        else
        {
            /* The body is a flexbody*/
            node = new Node(*system_in->flexBodies[bodyTitleB]->getNode(atoi(nodeB.c_str())));
        }

        system_in->groundNodes.push_back(node);
        system_in->groundNodesMap[consName] = node;
        p_nodeA = node;
        p_nodeA->setNumber(-1);
        p_nodeA->setThermalNumber(-1);
    }

    if (bodyTitleB == "GROUND")
    {
        Node* node = nullptr;

        if (system_in->subSystems.find(bodyTitleA) != system_in->subSystems.end())
        {
            /* If the body is a system */
            node = new Node(*system_in->subSystems[bodyTitleA]->getNode(atoi(nodeA.c_str())));
        }
        else if (system_in->rigidBodies.find(bodyTitleA) != system_in->rigidBodies.end())
        {
            /* If the body is a rigidbody */
            node = new Node(system_in->rigidBodies[bodyTitleA]->getNode(atoi(nodeA.c_str())));
        }
        else
        {
            /* the body is a flexbody */
            node = new Node(*system_in->flexBodies[bodyTitleA]->getNode(atoi(nodeA.c_str())));
        }

        system_in->groundNodes.push_back(node);
        system_in->groundNodesMap[consName] = node;
        p_nodeB = node;
        p_nodeB->setNumber(-1);
        p_nodeB->setThermalNumber(-1);
    }

    /* if it's not grounded */
    if (p_nodeA == nullptr)
    {
        if (system_in->subSystems.find(bodyTitleA) != system_in->subSystems.end())
        {
            p_nodeA = system_in->subSystems[bodyTitleA]->getNode(atoi(nodeA.c_str()));
        }
        else if (system_in->rigidBodies.find(bodyTitleA) != system_in->rigidBodies.end())
        {
            //if the body is a rigidbody
            p_nodeA = system_in->rigidBodies[bodyTitleA]->getNode(atoi(nodeA.c_str()));
        }
        else
        {
            //the body is a flexbody
            p_nodeA = system_in->flexBodies[bodyTitleA]->getNode(atoi(nodeA.c_str()));
        }
    }

    /* if it's not grounded */
    if (p_nodeB == nullptr)
    {
        if (system_in->subSystems.find(bodyTitleB) != system_in->subSystems.end())
        {
            p_nodeB = system_in->subSystems[bodyTitleB]->getNode(atoi(nodeB.c_str()));
        }
        else if (system_in->rigidBodies.find(bodyTitleB) != system_in->rigidBodies.end())
        {
            //if the body is a rigidbody
            p_nodeB = system_in->rigidBodies[bodyTitleB]->getNode(atoi(nodeB.c_str()));
        }
        else
        {
            //the body is a flexbody
            p_nodeB = system_in->flexBodies[bodyTitleB]->getNode(atoi(nodeB.c_str()));
        }
    }
}


void mknix::ReaderConstraints::outputConstraintNode(System* system_in,
        const std::string& consTitle,
        const std::string& nodeName,
        const std::string& bodyTitle,
        const std::string& node,
        std::size_t i)
{
    *output << "\t"
            << consTitle << "." << nodeName << ": "
            << bodyTitle << "." << node
            << ": " << system_in->constraints[consTitle]->getNode(i)->getNumber()
            << ": " << system_in->constraints[consTitle]->getNode(i)->getX()
            << ", " << system_in->constraints[consTitle]->getNode(i)->getY()
            << ", " << system_in->constraints[consTitle]->getNode(i)->getZ()
            << std::endl;
}


void mknix::ReaderConstraints::outputConstraintThermalNode(System* system_in,
        const std::string& consTitle,
        const std::string& nodeName,
        const std::string& bodyTitle,
        const std::string& node,
        std::size_t i)
{
    *output << "\t"
            << consTitle << "." << nodeName << ": "
            << bodyTitle << "." << node
            << ": " << system_in->constraintsThermal[consTitle]->getNode(i)->getNumber()
            << ": " << system_in->constraintsThermal[consTitle]->getNode(i)->getX()
            << ", " << system_in->constraintsThermal[consTitle]->getNode(i)->getY()
            << ", " << system_in->constraintsThermal[consTitle]->getNode(i)->getZ()
            << std::endl;
}

