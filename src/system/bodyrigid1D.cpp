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

#include "bodyrigid1D.h"

#include <simulation/simulation.h>

namespace mknix
{

/**
 * @brief Default constructor. Initializes the bar with zero moment of inertia.
 */
RigidBar::RigidBar()
    : RigidBody()
    , Jo(0)
{
}

/**
 * @brief Constructor with title, end nodes and mass.
 * @param title_in  Name identifier for this bar body.
 * @param nodeA_in  Pointer to node A (first endpoint).
 * @param nodeB_in  Pointer to node B (second endpoint).
 * @param mass_in   Total mass of the bar.
 */
RigidBar::RigidBar(std::string title_in,
                   Node * nodeA_in,
                   Node * nodeB_in,
                   double mass_in
                  )
    : RigidBody(title_in)
{
    this->mass = mass_in;
    // Read frameNodes from theReader:
    // It is assumed that correspond to extremes of bar and
    frameNodes.push_back(nodeA_in);
    frameNodes.push_back(nodeB_in);
    // We calculate bar properties and reposition frameNodes to the normalized position
    lenght = std::sqrt(std::pow(nodeB_in->getX() - nodeA_in->getX(), 2) +
                       std::pow(nodeB_in->getY() - nodeA_in->getY(), 2) +
                       std::pow(nodeB_in->getZ() - nodeA_in->getZ(), 2));
    rho = mass_in / lenght;

    this->localMassMatrix.resize(2 * dim, 2 * dim);
    this->externalForces.resize(2 * dim);

}

/**
 * @brief Destructor.
 */
RigidBar::~RigidBar()
{
}

/**
 * @brief Sets the rotational inertia about axis 0 (only one axis supported).
 * @param inertia_in Inertia value to assign.
 * @param axis       Axis index (only 0 is valid).
 */
void RigidBar::setInertia(double inertia_in, int axis)
{
    if (axis == 0)
    {
        Jo = inertia_in;
    }
    else
    {
        cerr << endl << "ERROR: Trying to set inertia out of bounds in RigidBar" << endl;
    }
}

/**
 * @brief Sets initial end-node positions from a centre-of-mass position and rotation angle.
 * @param position Vector [CoG_x, CoG_y, z, angle] defining the bar pose.
 */
void RigidBar::setPosition(std::vector<double>& position)
{
    // TODO: check vector size. Should have 3 elements: CoG_x, CoG_y, rotation angle
    this->frameNodes[0]->setX(position[0] - lenght / 2 * std::cos(position[3]));
    this->frameNodes[0]->setY(position[1] - lenght / 2 * std::sin(position[3]));
    this->frameNodes[1]->setX(position[0] + lenght / 2 * std::cos(position[3]));
    this->frameNodes[1]->setY(position[1] + lenght / 2 * std::sin(position[3]));
}

/**
 * @brief Builds the consistent local mass matrix for the bar (2-node linear element).
 */
void RigidBar::calcMassMatrix()
{
    this->localMassMatrix(0, 0) = 1. / 3. * mass;
    this->localMassMatrix(1, 1) = 1. / 3. * mass;
    this->localMassMatrix(Simulation::getDim(), Simulation::getDim()) = 1. / 3. * mass;
    this->localMassMatrix(Simulation::getDim() + 1, Simulation::getDim() + 1) = 1. / 3. * mass;

    this->localMassMatrix(0, Simulation::getDim()) = 1. / 6. * mass;
    this->localMassMatrix(Simulation::getDim(), 0) = 1. / 6. * mass;
    this->localMassMatrix(1, Simulation::getDim() + 1) = 1. / 6. * mass;
    this->localMassMatrix(Simulation::getDim() + 1, 1) = 1. / 6. * mass;

    if (Simulation::getDim() == 3)
    {
        this->localMassMatrix(2, 2) = 1. / 3. * mass;
        this->localMassMatrix(Simulation::getDim() + 2, Simulation::getDim() + 2) = 1. / 3. * mass;
        this->localMassMatrix(2, Simulation::getDim() + 2) = 1. / 6. * mass;
        this->localMassMatrix(Simulation::getDim() + 2, 2) = 1. / 6. * mass;
    }
}

/**
 * @brief Computes the gravitational load vector: half of -m*g at each end node.
 */
void RigidBar::calcExternalForces()
{
    this->externalForces(0) = -0.5 * mass * Simulation::getGravity(0);
    this->externalForces(1) = -0.5 * mass * Simulation::getGravity(1);
    if (Simulation::getDim() == 3)
    {
        this->externalForces(2) = -0.5 * mass * Simulation::getGravity(2);
    }
    this->externalForces(Simulation::getDim()) = -0.5 * mass * Simulation::getGravity(0);
    this->externalForces(Simulation::getDim() + 1) = -0.5 * mass * Simulation::getGravity(1);
    if (Simulation::getDim() == 3)
    {
        this->externalForces(Simulation::getDim() + 2) = -0.5 * mass * Simulation::getGravity(2);
    }
}


/**
 * @brief Adds a domain node to the bar, assigning the two frame nodes as support nodes
 *        and solving the 1D shape function.
 * @param node_in Pointer to the domain node to add.
 */
void RigidBar::addNode(Node * node_in)
{
    mknix::Body::addNode(node_in); // adds node_in to node vector
    this->nodes.back()->addSupportNode(this->frameNodes[0]);
    this->nodes.back()->addSupportNode(this->frameNodes[1]);
    this->nodes.back()->setShapeFunType("1D");
    this->nodes.back()->shapeFunSolve("1D", 1.);
}


// Node * RigidBar::getDomainNode( std::string name_in )
// {
//   if( name_in == "NODEA" ) return this->frameNodes[0];
//   else if( name_in == "NODEB" ) return this->frameNodes[1];
//   else cerr << "ERROR: NO NODE WITH THAT NAME IN BAR" << endl;
// }

/**
 * @brief Collects the frame nodes (bar endpoints) as boundary nodes.
 * @param boundary_nodes Vector to which the frame node pointers are appended.
 */
void RigidBar::writeBoundaryNodes(std::vector<Node *>& boundary_nodes)
{
    // We specialize only for the bar, as its frameNodes are its own boundary
    for (auto& node : this->frameNodes)
    {
        boundary_nodes.push_back(node);
    }
}

/**
 * @brief Builds the boundary connectivity using the frame nodes (bar endpoints).
 * @param connectivity_nodes Vector of node chains to which the bar boundary is appended.
 */
void RigidBar::writeBoundaryConnectivity(std::vector<std::vector<Node *> >& connectivity_nodes)
{
    // We specialize only for the bar, as its frameNodes are its own boundary
    connectivity_nodes.push_back(std::vector<Node *>());
    for (auto& node : this->frameNodes)
    {
        connectivity_nodes.back().push_back(node);
    }
}

}
