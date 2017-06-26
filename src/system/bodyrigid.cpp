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
#include "bodyrigid.h"

#include <simulation/simulation.h>

namespace mknix {

RigidBody::RigidBody()
        : Body()
        , computeEnergy(false)
        , dim(Simulation::getDim())
        , mass(0)
        , densityFactor(1)
{
}

RigidBody::RigidBody(std::string title_in)
        : Body(title_in)
        , computeEnergy(false)
        , dim(Simulation::getDim())
        , mass(0)
        , densityFactor(1)
{
}

RigidBody::~RigidBody()
{
}


void RigidBody::assembleMassMatrix(lmx::Matrix<data_type>& globalMass)
{
    class Body;

    int frameNodesSize = frameNodes.size();
    int i, j, m, n;

    for (i = 0; i < frameNodesSize; ++i) {
        for (j = 0; j < frameNodesSize; ++j) {
            for (m = 0; m < Simulation::getDim(); ++m) {
                for (n = 0; n < Simulation::getDim(); ++n) {
                    globalMass(Simulation::getDim() * frameNodes[i]->getNumber() + m,
                               Simulation::getDim() * frameNodes[j]->getNumber() + n)
                            += localMassMatrix.readElement(Simulation::getDim() * i + m,
                                                           Simulation::getDim() * j + n);
                }
            }
        }
    }
}

void RigidBody::assembleExternalForces
        (lmx::Vector<data_type>& globalExternalForces)
{
    int frameNodesSize = frameNodes.size();
    int i, m;

    for (i = 0; i < frameNodesSize; ++i) {
        for (m = 0; m < Simulation::getDim(); ++m) {
            globalExternalForces(Simulation::getDim() * frameNodes[i]->getNumber() + m)
                    += externalForces.readElement(Simulation::getDim() * i + m);
        }
    }
}

// Node* RigidBody::getDomainNode(std::string node_in )
// {
//   if( node_in == "CoG" ) return this->getNode(0);
//   else if( node_in == "x1" ) return this->getNode(1); // for testing proposes, returns first orientation node
//   else return this->getNode( atoi(node_in.c_str()) );
// }


void RigidBody::setOutput(std::string outputType_in)
{
    if (outputType_in == "ENERGY") computeEnergy = 1;
}

Node * RigidBody::getNode(int node_number)
{
    if (node_number == -1) {
        return this->frameNodes[0];
    }
    else if (node_number == -2) {
        return this->frameNodes[1];
    }
    else if (node_number == -3 && this->frameNodes.size() > 2) {
        return this->frameNodes[2];
    }
    else if (node_number == -4 && this->frameNodes.size() > 3) {
        return this->frameNodes[3];
    }
    else if (node_number >=0 && (size_t)node_number < this->nodes.size()) {
        return this->nodes[node_number];
    }
    else {
        cerr << "ERROR: NODE " << node_number << " not found!!!" << endl;
        return nullptr;
    }
}


void RigidBody::outputStep
        (const lmx::Vector<data_type>& q, const lmx::Vector<data_type>& qdot)
{
    Body::outputStep();

    auto nodesSize = nodes.size();
    if (nodesSize > 0) {
        domainConf.push_back(new lmx::Vector<data_type>(dim * nodesSize));
        for (auto i = 0u; i < nodesSize; ++i) {
            for (auto j = 0u; j < (size_t)dim; ++j) {
                domainConf.back()->writeElement(nodes[i]->getConf(j), dim * i + j);
            }
        }
    }
    if (computeEnergy) {
        // TODO: compute energy for each kind of rigid body
        /*potential, kinetic, elastic, total*/
//     energy.push_back( new lmx::Vector<data_type>( 4 ) );
//     energy.back()->fillIdentity( 0. );
//
//     frameNodesSize = frameNodes.size();
//
//     for( i=0; i<frameNodesSize; ++i )
//     {
//       for( n=0; n<Simulation::getDim() ; ++n)
//       {
//         /*potential*/
//         energy.back()->operator()(0)
//             -= q.readElement( Simulation::getDim()*frameNodes[i]->getNumber() + n )
//             * externalForces.readElement(Simulation::getDim()*i + n);
//
//         /*kinetic*/
//         for (j=0; j<frameNodesSize; ++j)
//         {
//           for (m=0; m<Simulation::getDim(); ++m)
//           {
//             energy.back()->operator()(1)
//             += 0.5
//             * qdot.readElement( Simulation::getDim()*frameNodes[i]->getNumber() + n )
//             * localMassMatrix.readElement( Simulation::getDim()*i + n, Simulation::getDim()*j + m)
//             * qdot.readElement( Simulation::getDim()*frameNodes[j]->getNumber() + m );
//           }
//         }
//       }
//     }
//     //no Elastic
//     /*Total*/
//     energy.back()->operator()(3) = energy.back()->readElement(0)
//                                  + energy.back()->readElement(1)
//                                  + energy.back()->readElement(2);
    }
}

void RigidBody::outputStep (const VectorX<data_type>& q, const VectorX<data_type>& qdot)
{//TODO
  /*  Body::outputStep();

    auto nodesSize = nodes.size();
    if (nodesSize > 0) {
        domainConf.push_back(new lmx::Vector<data_type>(dim * nodesSize));
        for (auto i = 0u; i < nodesSize; ++i) {
            for (auto j = 0u; j < (size_t)dim; ++j) {
                domainConf.back()->writeElement(nodes[i]->getConf(j), dim * i + j);
            }
        }
    }
    if (computeEnergy) {
        // TODO: compute energy for each kind of rigid body

    }*/
}


void RigidBody::outputStep(const lmx::Vector<data_type>& q)
{
    Body::outputStep();

    int nodesSize = nodes.size();
    int i, j;
    if (nodesSize > 0) {
        domainConf.push_back(new lmx::Vector<data_type>(dim * nodesSize));
        for (i = 0; i < nodesSize; ++i) {
            for (j = 0; j < dim; ++j) {
                domainConf.back()->writeElement(nodes[i]->getConf(j), dim * i + j);
            }
        }
    }
    if (computeEnergy) {
        // TODO: compute energy for each kind of rigid body
//     energy.push_back( new lmx::Vector<data_type>( 4 ) ); //potential, kinetic, elastic
//
//     energy.back()->fillIdentity( 0. );
//
//     nodesSize = nodes.size();
//     for (i=0; i<nodesSize; ++i){
//       for (j=0; j<dim; ++j){
//         energy.back()->operator()(0) //potential
//             -= q.readElement( Simulation::getDim()*nodes[i]->getNumber() + j )
//             * externalForces.readElement(Simulation::getDim()*i + j);
//       }
//     }
//     //no kinetic, no Elastic
//     energy.back()->operator()(3) = energy.back()->readElement(0)
// 				 + energy.back()->readElement(1)
// 				 + energy.back()->readElement(2); //total
//
    }
}

void RigidBody::outputStep(const VectorX<data_type>& q)
{//TODO complete
    Body::outputStep();
/*
    int nodesSize = nodes.size();
    int i, j;
    if (nodesSize > 0) {
        domainConf.push_back(new VectorX<data_type>(dim * nodesSize));
        for (i = 0; i < nodesSize; ++i) {
            for (j = 0; j < dim; ++j) {
                *(domainConf.back())[ dim * i + j] = nodes[i]->getConf(j);
            }
        }
    }
    if (computeEnergy) {
        // TODO: compute energy for each kind of rigid body
//
    }*/
}

void RigidBody::outputToFile(std::ofstream * outFile)
{
    Body::outputToFile(outFile);

    std::vector<lmx::Vector<data_type> *>::iterator itDomainNodes;

    if (domainConf.size() > 0) {
        *outFile << "DOMAIN " << title << " " << nodes.size() << endl;
        for (auto& domainNode : domainConf) {
//       vectorSize = (*itDomainNodes)->size();
//       for( i=0; i<vectorSize; ++i){
            for (auto i = 0u; i < nodes.size(); ++i) {
                for (auto j = 0u; j < static_cast<size_t>(dim); ++j) {
//         *outFile << (*itDomainNodes)->readElement(i) << " ";
                    *outFile << domainNode->readElement(dim * i + j) << " ";
                }
                *outFile << endl;
            }
            *outFile << endl;
        }
    }
    if (energy.size() > 0) {
        *outFile << "ENERGY " << title << endl;
        for (auto& el : energy) {
            auto vectorSize = el->size();
            for (auto i = 0u; i < vectorSize; ++i) {
                *outFile << el->readElement(i) << " ";
            }
            *outFile << endl;
        }
    }
}

void RigidBody::writeBodyInfo(std::ofstream * outFile)
{
    *outFile << "\t" << this->getType() << " "
    << this->title << " "
    << "0 "; // should give material information
    std::vector<Node *>::iterator it_frameNodes;
    for (it_frameNodes = frameNodes.begin();
         it_frameNodes != frameNodes.end();
         ++it_frameNodes) {
        *outFile << (*it_frameNodes)->getNumber() << " ";
    }
    *outFile << endl;
}


void RigidBody::writeBoundaryNodes(std::vector<Point *>& boundary_nodes)
{
//   std::vector<Node*>::iterator it_nodes;
//   for( it_nodes = this->nodes.begin();
//        it_nodes!= nodes.end();
//        ++it_nodes )
//   {
//     boundary_nodes.push_back( *it_nodes );
//   }
}

void RigidBody::writeBoundaryConnectivity(std::vector<std::vector<Point *> >& connectivity_nodes)
{
//   std::vector<Node*>::iterator it_nodes;
//   connectivity_nodes.push_back( std::vector< Node* >() );
//   for( it_nodes = this->nodes.begin();
//        it_nodes!= nodes.end();
//        ++it_nodes )
//   {
//     connectivity_nodes.back().push_back( *it_nodes );
//   }
}


}
