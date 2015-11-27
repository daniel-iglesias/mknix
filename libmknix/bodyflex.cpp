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
#include "bodyflex.h"
#include "cell.h"
#include "simulation.h"

namespace mknix {

FlexBody::FlexBody()
        : Body()
        , formulation("NONLINEAR")
        , computeStress(0)
        , computeEnergy(0)
{
}


FlexBody::FlexBody(std::string title_in)
        : Body(title_in)
        , formulation("NONLINEAR")
        , computeStress(0)
        , computeEnergy(0)
{
}


FlexBody::~FlexBody()
{
}

void FlexBody::initialize()
{
    mknix::Body::initialize();

    for (auto& point : points) {
        point->findSupportNodes(this->nodes);
        point->shapeFunSolve(1.03);
    }
    for (auto& point : bodyPoints) {
        point->findSupportNodes(this->nodes);
        point->shapeFunSolve(1.03);
    }

}


void FlexBody::addBodyPoint( /*const*/ Point * point_in, std::string method_in)
{   // Creates or redefines a point of the integration domain:
//       cout << "addBodyPoint, with points(size) = " << this->points.size() << endl;
    this->bodyPoints.push_back(point_in);
//     cout << "POINT ADDED in FLEXBODY:" << points.back()->getNumber()
//          << ", Support nodes: " << points.back()->getSupportSize(0) << endl;
    if (method_in == "RPIM" || method_in == "RBF") {
        method_in = "RBF";
    } else if (method_in == "EFG" || method_in == "MLS") {
        method_in = "MLS";
    } else { cerr << "ERROR: UNKNOWN METHOD FOR BODY POINT. " << method_in << " WAS DEFINED."; }
    bodyPoints.back()->setShapeFunType(method_in);
}

void FlexBody::addPoint( /*const*/ Node * node_in) // problems with munmap_chunk(). Use the later function
{   // Creates or redefines a point, not part of the calculation domain:
//       cout << "addPoint, with points(size) = " << this->points.size() << endl;
    this->points.push_back(new Node(node_in));
    cout << "POINT ADDED in FLEXBODY:" << points.back()->getNumber()
    << ", Support nodes: " << points.back()->getSupportSize(0) << endl;
    points.back()->setShapeFunType("MLS");
}


void FlexBody::addPoint(int nodeNumber, double x, double y, double z, double alpha, double dc)
{   // Creates or redefines a point, not part of the calculation domain:
//       cout << "addPoint, with points(size) = " << this->points.size() << endl;
    this->points.push_back(new Node(nodeNumber, x, y, z));
    points.back()->setAlphai(alpha);
    points.back()->setDc(dc);

    points.back()->findSupportNodes(this->nodes);
    cout << "POINT ADDED in FLEXBODY:" << points.back()->getNumber()
    << ", Support nodes: " << points.back()->getSupportSize(0) << endl;
    points.back()->setShapeFunType("MLS");
    points.back()->shapeFunSolve("MLS", 1.03);
}


/**
 * @brief Activates a flag for output data at the end of the analysis.
 *
 * @see outputToFile()
 * @see outputStep()
 * @param outputType_in Keyword of the flag. Options are: [STRESS, ENERGY]
 * @return void
 **/
void FlexBody::setOutput(std::string outputType_in)
{
    if (outputType_in == "STRESS") {
        computeStress = 1;
    } else if (outputType_in == "ENERGY") {
        computeEnergy = 1;
    }
}


/**
 * @brief Streams the data stored during the analysis to a file.
 *
 * @param outFile Output files
 * @return void
 **/
void FlexBody::outputToFile(std::ofstream * outFile)
{
    Body::outputToFile(outFile);

    if (computeStress) {
        *outFile << "STRESS " << title << endl;
        for (auto& stress: stresses) {
            auto vectorSize = stress.size();
            for (auto i = 0u; i < vectorSize; ++i) {
                *outFile << stress.readElement(i) << " ";
            }
            *outFile << endl;
        }
    }

    if (computeEnergy) {
        *outFile << "ENERGY " << title << endl;
        for (auto& energy : energies) {
            auto vectorSize = energy->size();
            for (auto i = 0u; i < vectorSize; ++i) {
                *outFile << energy->readElement(i) << " ";
            }
            *outFile << endl;
        }
    }
}

void FlexBody::writeBodyInfo(std::ofstream * outFile)
{
// DONE: set type to MESH (alt:GALERKIN), export the connectivity and read it in mknixPost
//   Implemented 2D; 3D shows strange connectivity when read in mknixpost.
//   Previous code is use for 3D meshes till fixed.
    if (Simulation::getDim() == 2) {
        *outFile << "\t" << "MESH" << " "
        << this->title << endl;
        if (bodyPoints.size() == 0) { //regular one domain meshes
            *outFile << "\t" << "\t" << "NODES "
            << this->getNode(0)->getNumber() << " "
            << this->getLastNode()->getNumber() << endl;
        }
        else { // multiple domains
            *outFile << "\t" << "\t" << "NODES "
            << this->getBodyPoint(0)->getNumber() << " "
            << this->getLastBodyPoint()->getNumber() << endl;
        }
        *outFile << "\t" << "\t" << "CELLS "
        << this->cells.size() << endl;
        for (auto& cell : cells) {
            cell.second->outputConnectivityToFile(outFile);
        }
        *outFile << "\n\t" << "END" << "MESH" << endl;
    }
    else {
        //    if(itFlexBodies->second->getType() == "MESHFREE"){
        //    *outFile << "\t" << itFlexBodies->second->getType() << " "
        *outFile << "\t" << "MESHFREE" << " "
        << this->title << endl;
        if (bodyPoints.size() == 0) { //regular one domain meshes
            *outFile << "\t" << "\t" << "NODES "
            << this->getNode(0)->getNumber() << " "
            << this->getLastNode()->getNumber() << endl;
        }
        else { // multiple domains
            *outFile << "\t" << "\t" << "NODES "
            << this->getBodyPoint(0)->getNumber() << " "
            << this->getLastBodyPoint()->getNumber() << endl;
        }

        //TODO: implement for 3D. For the moment
        if (Simulation::getDim() == 2) {
            *outFile << "\t" << "\t" << "BOUNDARY "
            << getBoundarySize() << endl;
            if (getBoundarySize() > 0) {
                int firstNode = getBoundaryFirstNode()->getNumber();
                //int actualNode = firstNode;
                //Point * p_actualNode = getBoundaryFirstNode();
                *outFile << "\t" << "\t" << firstNode << " ";
                //    for( int i=1; i<boundary.size(); ++i )
                //    {
                //      actualNode = p_actualNode->getNumber();
                //      *outFile << actualNode << " " ;
                //      cout <<actualNode << endl;
                //      p_actualNode = getBoundaryNextNode( p_actualNode );
                //    }
                // Following code works only with either MESH or one CELLS file. Testing an alternative below...
                //    for(int i=0; i<getBoundarySize(); ++i){
                //      p_actualNode = getBoundaryNextNode( p_actualNode );
                //      actualNode = p_actualNode->getNumber();
                for (auto& boundary : linearBoundary) {
                    auto p_actualNode = boundary.second;
                    auto actualNode = p_actualNode->getNumber();
                    *outFile << actualNode << " ";
                    cout << "BoundaryNextNode: " << actualNode << endl;
                }
                //    while( actualNode != firstNode );
                *outFile << endl;// << "\t" << "\t" << "ENDBOUNDARY" << endl;
            }
        }
        *outFile << "\n\t" << "END" << "MESHFREE" << endl;
    }

//    *outFile << "\t" << "END" << itFlexBodies->second->getType() << endl;
}


void FlexBody::writeBoundaryNodes(std::vector<Point *>& boundary_nodes)
{
    for (auto boundary : linearBoundary) {
        boundary_nodes.push_back(boundary.first);
    }
}

void FlexBody::writeBoundaryConnectivity(std::vector<std::vector<Point *> >& connectivity_nodes)
{
    connectivity_nodes.push_back(std::vector<Point *>());
    Point * first_node = linearBoundary.begin()->first;
    connectivity_nodes.back().push_back(first_node);
    Point * next_node = linearBoundary[first_node];
    while (next_node != first_node) {
        connectivity_nodes.back().push_back(next_node);
        next_node = linearBoundary[next_node];
    }
}

}
