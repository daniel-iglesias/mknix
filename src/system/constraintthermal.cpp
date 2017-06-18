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
#include "constraintthermal.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix {

ConstraintThermal::ConstraintThermal()
    : Constraint()
{
}


ConstraintThermal::ConstraintThermal(double& alpha_in, std::string& method_in)
    : Constraint(alpha_in, method_in, 1)
{
}


ConstraintThermal::~ConstraintThermal()
{
}

void ConstraintThermal::assembleInternalForces
(VectorX< data_type > & globalInternalForces)
{
    int nodesSize = nodes.size();
    int i, m;
    size_t k, counter(0);
    for (i=0; i<nodesSize; ++i) {
        for (k=0; k < nodes[i]->getSupportSize(0); ++k) {
            if (nodes[i]->getThermalNumber() >= 0 ) {
                for (m=0; m<dim; ++m) {
//           cout << endl <<
//                   "dim*nodes[i]->getSupportNodeNumber(0,k) + m = " <<
//                   dim*nodes[i]->getSupportNodeNumber(0,k) + m << endl <<
//                   "dim*counter + m = " <<
//                   dim*counter + m << endl;
                    globalInternalForces( dim*nodes[i]->getSupportNodeNumber(0,k) + m)
                    += internalForces.readElement(dim*counter + m);
                }
            }
            ++counter;
        }
//        for (m=0; m<dim; ++m){
//          globalInternalForces( dim*nodes[i]->getThermalNumber() + m)
//            += internalForces.readElement(dim*i + m);
    }
}

void ConstraintThermal::assembleTangentMatrix(SparseMatrix< data_type > & globalTangent)
{
    int nodesSize = nodes.size();
    size_t supportNodesSize_i, supportNodesSize_j;
    int i, j, m, n;
    size_t k, l, counter_i(0), counter_j(0);
    for ( i=0; i<nodesSize; ++i ) {
        supportNodesSize_i = nodes[i]->getSupportSize(0);
        for (k=0; k < supportNodesSize_i; ++k) {
            for ( j=0; j<nodesSize; ++j ) {
                supportNodesSize_j = nodes[j]->getSupportSize(0);
                for (l=0; l < supportNodesSize_j; ++l) {
                    if (nodes[i]->getThermalNumber() >= 0 ) {
                        if (nodes[j]->getThermalNumber() >= 0 ) {
                            for ( m=0; m<dim; ++m ) {
                                for ( n=0; n<dim; ++n ) {
                                    globalTangent(dim*nodes[i]->getSupportNodeNumber(0,k) + m,
                                                  dim*nodes[j]->getSupportNodeNumber(0,l) + n)
                                    += stiffnessMatrix.readElement(dim*counter_i + m,
                                                                   dim*counter_j + n);
                                }
                            }
                        }
                    }
                    ++counter_j;
                }
//              for ( m=0; m<dim; ++m ){
//                for ( n=0; n<dim; ++n ){
//                  globalTangent(dim*nodes[i]->getThermalNumber() + m,
//                                dim*nodes[j]->getThermalNumber() + n)
//                  += stiffnessMatrix.readElement(dim*i + m,
//                                                 dim*j + n);
//                }
//              }
            }
            counter_j=0;
            ++counter_i;
        }
    }
}

}
