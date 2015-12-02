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
#include "systemchain.h"
#include "bodyrigid1D.h"
#include "constraintdistance.h"
#include "constraintfixedcoordinates.h"

#include <simulation/simulation.h>

namespace mknix {

SystemChain::SystemChain()
{
}


SystemChain::SystemChain(const char * title_in) : System(title_in)
{
}


SystemChain::~SystemChain()
{
}

Node * SystemChain::getNode(int node_i)
{
    if (node_i == 0) return this->rigidBodies.begin()->second->getNode(-1);
    if (node_i == 1) return this->rigidBodies.rbegin()->second->getNode(-2);
    return nullptr;
}


void SystemChain::populate(Simulation * theSimulation, std::string& energyKeyword)
{
    std::stringstream barName;
    std::string rigidTitle;
    int lastNode;
    double xA, yA, zA, xB, yB, zB;
    currentLength = std::pow(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) + std::pow(z1 - z0, 2), 0.5);
    // Define or update the timeLenght[0] value
    timeLenghts[0] = currentLength;

    for (int i = 0; i < segments; ++i) {
        cout << i * (length / segments) << " < " << currentLength << " ?" << endl;
        if (i * (length / segments) < currentLength) {
            xA = x0 + i * (length / segments) / currentLength * (x1 - x0);
            yA = y0 + i * (length / segments) / currentLength * (y1 - y0);
            zA = z0 + i * (length / segments) / currentLength * (z1 - z0);
        }
        else {
            xA = x1;
            yA = y1;
            zA = z1;
        }
        if ((i + 1) * (length / segments) < currentLength) {
            xB = x0 + (i + 1) * (length / segments) / currentLength * (x1 - x0);
            yB = y0 + (i + 1) * (length / segments) / currentLength * (y1 - y0);
            zB = z0 + (i + 1) * (length / segments) / currentLength * (z1 - z0);
        }
        else {
            xB = x1;
            yB = y1;
            zB = z1;
        }

        lastNode = theSimulation->nodes.end()->first;
        theSimulation->nodes[lastNode] = new Node(lastNode, xA, yA, zA);
        theSimulation->nodes[lastNode + 1] = new Node(lastNode + 1, xB, yB, zB);

        barName.str(std::string());
        barName << this->title << i;
        rigidTitle = barName.str();

        this->rigidBodies[rigidTitle] = new RigidBar(rigidTitle,
                                                     theSimulation->nodes[lastNode + 0],
                                                     theSimulation->nodes[lastNode + 1],
                                                     mass / segments
        );

        localConstraintDistance[rigidTitle] = new ConstraintDistance(
                theSimulation->nodes[lastNode + 0],
                theSimulation->nodes[lastNode + 1],
                Simulation::alpha,
                Simulation::constraintMethod
        );
        // Storing local pointers for classifing specific ConstraintDistance
        // Can be optimized using castings below
        this->constraints[rigidTitle] = localConstraintDistance[rigidTitle];

        this->rigidBodies[rigidTitle]->setOutput(energyKeyword); //chapuza
        // Spherical joints between the bars:
        if (i > 0) { // No joint with only one bar
            barName << this->title << i + 1;
            rigidTitle = barName.str();
            this->constraints[rigidTitle] = new ConstraintFixedCoordinates(
                    theSimulation->nodes[lastNode + 0],
                    theSimulation->nodes[lastNode - 1],
                    Simulation::alpha,
                    Simulation::constraintMethod
            );
        }
    }

}


void SystemChain::update(double theTime)
{
    std::stringstream barName;
    // compute new length. Temp store
    double updateLength = interpolate1D(theTime, timeLenghts);
    cout << "currentLength = " << currentLength
    << ", updateLength = " << updateLength << endl;

    if (updateLength != currentLength) {
        if (updateLength > length || updateLength < 0) {
            cerr << "ERROR: Trying to extend CHAIN beyond its lenght!!!" << endl;
        }
        else { // Commented are possible optimizations (opt)
//       if( updateLength > currentLength ){ // EXTENSION for opt
            // Locate bar. Compute index with currentLength
            //int firstConstraint = floor(currentLength * segments / length);
            int lastConstraint = ceil(updateLength * segments / length);
            for (int i = 0; i < lastConstraint; ++i) { // for opt: Should say i=firstConstraint
                barName.str(std::string());
                barName << this->title << i;
                if (i < lastConstraint - 1) { // Extend constraint to maximum
                    this->localConstraintDistance[barName.str()]
                            ->setLenght(length / segments);
                }
                else { // Must be extended precisely
                    this->localConstraintDistance[barName.str()]
                            ->setLenght(updateLength - i * length / segments - 1E-8);
                }
            }
//       }
//       if( updateLength < currentLength ){ // COMPRESION
//         // Locate bar. Compute index with currentLength
//         int firstConstraint = ceil( currentLength*segments/length );
//         int lastConstraint = floor( updateLength*segments/length );
//         for (int i=firstConstraint; i > lastConstraint; --i){
//           barName.str(std::string());
//           barName << this->title << i;
//           if(i>lastConstraint+1){ // Compress constraint to minimum
//             this->localConstraintDistance[barName.str()]->setLenght(0.);
//           }
//           else{ // Must be extended precisely
//             this->localConstraintDistance[barName.str()]
//             ->setLenght(updateLength-lastConstraint*length/segments);
//           }
//         }
//       }
            // Update currentLength
            currentLength = updateLength;
        }
    }
}


}
