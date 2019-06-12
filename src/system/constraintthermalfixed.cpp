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
#include "constraintthermalfixed.h"

#include <core/node.h>
#include <simulation/simulation.h>

namespace mknix
{

ConstraintThermalFixed::ConstraintThermalFixed()
    : ConstraintThermal()
{
}


ConstraintThermalFixed::ConstraintThermalFixed( Node* a_in, Node* b_in, double& alpha_in, std::string& method_in )
    : ConstraintThermal(alpha_in, method_in)
{
    nodes.push_back( a_in );
    nodes.push_back( b_in );
    size_type total_support_nodes = a_in->getSupportSize(0) + b_in->getSupportSize(0);

//   cout << "node A: " << nodes[0]->getNumber() << endl;
//   cout << "node A support size: " << a_in->getSupportSize(0) << endl;
//   cout << "node B: " << nodes[1]->getNumber() << endl;
//   cout << "node B support size: " << b_in->getSupportSize(0) << endl;

    Tt = nodes[1]->getTemp() - nodes[0]->getTemp() ;

    this->stiffnessMatrix.resize(total_support_nodes*dim,total_support_nodes*dim);
    this->internalForces.resize(total_support_nodes*dim);

    this->lambda.resize(1); // dimension is DoF, so only temperature
    this->lambda[0]=0.0;

    this->phi.resize(1);

    this->phi_q.resize(1);
    this->phi_q[0].resize(total_support_nodes);

    this->phi_qq.resize(1);
    this->phi_qq[0].resize(total_support_nodes,total_support_nodes);
}


ConstraintThermalFixed::~ConstraintThermalFixed()
{
}

void ConstraintThermalFixed::calcPhi()
{
    Tt = nodes[1]->getTemp() ;

    this->phi[0] = Tt - nodes[0]->getTemp() ;

//  cout << endl << "Tt = " << nodes[1]->getTemp() << " - " << nodes[0]->getTemp()<< " = " << Tt << endl;
//  cout << endl << "Phi in fixedcoordinates = " << phi[0] << endl;
}

void ConstraintThermalFixed::calcPhiq()
{
    size_type i,j;
    for(i=0; i<nodes[0]->getSupportSize(0); ++i)
    {
        this->phi_q[0](i+0) = -nodes[0]->getShapeFunValue(0,i);
    }
    // Now i=supportNodes[0].size()
    for(j=0; j<nodes[1]->getSupportSize(0); ++j)
    {
        this->phi_q[0](i+j) = nodes[1]->getShapeFunValue(0,j);
    }
//  this->phi_q[0](0) = -1.0;
//  this->phi_q[0](1) =  1.0;
//
}

void ConstraintThermalFixed::calcPhiqq()
{
}

bool ConstraintThermalFixed::checkAugmented()
{
    // TODO: Define as parameters in the input file: Aug-stage T limit (Tt), max_iter_augmented and delta_T for stage1, and delta_T for stage2.
    if (method == "AUGMENTED")
    {
//         if(Tt>300){
        //         double energy(0);
        double delta(0);
        for (auto i = 0u; i < phi.size(); ++i)
        {
            delta += std::abs(phi[i]);
            //             energy += phi[i] * phi[i];
            //             cout << endl << "Phi: " << phi[0] << "\t\t\t\t";
        }
        //         energy *= 0.5 * alpha;
        //         if (energy <= 2E5) {
        if ( (delta <= 10) && (Tt>=300) )
        {
            iter_augmented=0;
            return 1;
        }
        if ( ( (iter_augmented == 2) || (delta <= 20) ) && (Tt<300) )
        {
            iter_augmented=0;
            return 1;
        }
        else
        {
            //             cout << endl << "energy: " << energy << "\t\t\t\t";
            for (auto i = 0u; i < phi.size(); ++i)
            {
                lambda[i] += alpha * phi[i];
            }
            ++iter_augmented;
            return 0;
        }
    }
    else
    {
        iter_augmented=0;
        return 1;
    }
}

}

