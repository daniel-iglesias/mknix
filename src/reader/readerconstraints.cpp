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

#include <sstream>
#include <stdexcept>
 
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
            double temperature = theSimulation->getInitialTemperature();
            // bool hasTemperature = false;
 
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
                else if (keyword == "TEMPERATURE")
                {
                    *input >> temperature;
                    // hasTemperature = true;
                    *output << "TEMPERATURE: " << temperature << std::endl;
                }
            }
            this->assignConstraintNodes(system_in, consTitle, bodyTitleA, nodeA, bodyTitleB, nodeB);

            theSimulation->setThermalNodeInitialTemperature(p_nodeA, temperature);
            theSimulation->setThermalNodeInitialTemperature(p_nodeB, temperature);

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
 
    bodyTitle.clear();
    node.clear();
 
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

    if (bodyTitle == "GROUND")
    {
        cout << "NODE read: " << bodyTitle << endl;
        return;
    }

    /* Node id is the next whitespace-delimited token (avoids swallowing NODEB on one line). */
    if (a != '\n' && *input)
    {
        *input >> node;
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
    auto bodyInventory = [&]() {
        std::ostringstream stream;
        stream << "Available bodies with node counts:";

        if (system_in->subSystems.empty()
            && system_in->rigidBodies.empty()
            && system_in->flexBodies.empty())
        {
            stream << " none";
            return stream.str();
        }

        for (const auto& body : system_in->subSystems)
        {
            stream << "\n  SYSTEM " << body.first
                   << " (" << body.second->getNumberOfNodes() << " nodes)";
        }

        for (const auto& body : system_in->rigidBodies)
        {
            stream << "\n  RIGIDBODY " << body.first
                   << " (" << body.second->getNodesSize() << " nodes)";
        }

        for (const auto& body : system_in->flexBodies)
        {
            stream << "\n  FLEXBODY " << body.first
                   << " (" << body.second->getNodesSize() << " nodes)";
        }

        return stream.str();
    };

    auto throwLookupError = [&](const std::string& message) -> void {
        std::ostringstream stream;
        stream << message << '\n' << bodyInventory();
        throw std::runtime_error(stream.str());
    };

    auto parseNodeIndex = [&](const std::string& nodeName,
                              const std::string& bodyTitle,
                              const std::string& role) -> int {
        try
        {
            std::size_t processed = 0;
            const int index = std::stoi(nodeName, &processed);
            if (processed != nodeName.size())
            {
                throw std::invalid_argument(nodeName);
            }
            return index;
        }
        catch (const std::exception&)
        {
            throwLookupError("ERROR: invalid node " + role + " '" + nodeName
                             + "' for body '" + bodyTitle + "' in constraint '" + consName + "'.");
        }

        return 0;
    };

    auto resolveBodyNode = [&](const std::string& bodyTitle,
                               const std::string& nodeName,
                               const std::string& role) -> Node* {
        const int nodeIndex = parseNodeIndex(nodeName, bodyTitle, role);

        if (system_in->subSystems.find(bodyTitle) != system_in->subSystems.end())
        {
            auto* body = system_in->subSystems.at(bodyTitle);
            if (nodeIndex < 0 || static_cast<std::size_t>(nodeIndex) >= body->getNumberOfNodes())
            {
                throwLookupError("ERROR: node " + role + " index " + nodeName
                                 + " is out of range for system body '" + bodyTitle
                                 + "' in constraint '" + consName + "'.");
            }
            return body->getNode(static_cast<std::size_t>(nodeIndex));
        }

        if (system_in->rigidBodies.find(bodyTitle) != system_in->rigidBodies.end())
        {
            auto* body = system_in->rigidBodies.at(bodyTitle);
            if (nodeIndex < 0 || static_cast<std::size_t>(nodeIndex) >= static_cast<std::size_t>(body->getNodesSize()))
            {
                throwLookupError("ERROR: node " + role + " index " + nodeName
                                 + " is out of range for rigid body '" + bodyTitle
                                 + "' in constraint '" + consName + "'.");
            }
            return body->getNode(nodeIndex);
        }

        if (system_in->flexBodies.find(bodyTitle) != system_in->flexBodies.end())
        {
            auto* body = system_in->flexBodies.at(bodyTitle);
            if (nodeIndex < 0 || static_cast<std::size_t>(nodeIndex) >= static_cast<std::size_t>(body->getNodesSize()))
            {
                throwLookupError("ERROR: node " + role + " index " + nodeName
                                 + " is out of range for flex body '" + bodyTitle
                                 + "' in constraint '" + consName + "'.");
            }
            return body->getNode(nodeIndex);
        }

        throwLookupError("ERROR: body '" + bodyTitle + "' not found for " + role
                         + " '" + nodeName + "' in constraint '" + consName + "'.");
    };

    if (bodyTitleA == "GROUND")
    {
        Node* node = new Node(*resolveBodyNode(bodyTitleB, nodeB, "NODEB"));
 
        system_in->groundNodes.push_back(node);
        system_in->groundNodesMap[consName] = node;
        p_nodeA = node;
        p_nodeA->setNumber(-1);
        p_nodeA->setThermalNumber(-1);
    }
 
    if (bodyTitleB == "GROUND")
    {
        Node* node = new Node(*resolveBodyNode(bodyTitleA, nodeA, "NODEA"));
 
        system_in->groundNodes.push_back(node);
        system_in->groundNodesMap[consName] = node;
        p_nodeB = node;
        p_nodeB->setNumber(-1);
        p_nodeB->setThermalNumber(-1);
    }
 
    /* if it's not grounded */
    if (p_nodeA == nullptr)
    {
        p_nodeA = resolveBodyNode(bodyTitleA, nodeA, "NODEA");
    }
 
    /* if it's not grounded */
    if (p_nodeB == nullptr)
    {
        p_nodeB = resolveBodyNode(bodyTitleB, nodeB, "NODEB");
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
            << bodyTitle;

    if (!node.empty())
    {
        *output << "." << node;
    }

    *output << ": " << system_in->constraints[consTitle]->getNode(i)->getNumber()
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
            << bodyTitle;

    if (!node.empty())
    {
        *output << "." << node;
    }

    *output << ": " << system_in->constraintsThermal[consTitle]->getNode(i)->getNumber()
            << ": " << system_in->constraintsThermal[consTitle]->getNode(i)->getX()
            << ", " << system_in->constraintsThermal[consTitle]->getNode(i)->getY()
            << ", " << system_in->constraintsThermal[consTitle]->getNode(i)->getZ()
            << std::endl;
}