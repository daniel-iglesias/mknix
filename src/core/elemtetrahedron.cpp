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
#include "elemtetrahedron.h"
#include "material.h"
#include "node.h"
#include "gausspoint3D.h"

namespace mknix {

ElemTetrahedron::ElemTetrahedron()
        : CellTetrahedron()
{
}


ElemTetrahedron::ElemTetrahedron(Material& material_in,
                                 double alpha_in,
                                 int nGPoints_in,
                                 Node * n1_in,
                                 Node * n2_in,
                                 Node * n3_in,
                                 Node * n4_in
)
        : CellTetrahedron(material_in,
                          std::string(),
                          alpha_in,
                          nGPoints_in,
                          n1_in,
                          n2_in,
                          n3_in,
                          n4_in
) { }


ElemTetrahedron::~ElemTetrahedron()
{
}


void ElemTetrahedron::initialize(std::vector<Node *>& nodes_in)
{
//   cout << "CellTriang points " << this->points << endl;
    gPoints_MC.clear();
    this->createGaussPoints_MC();

    // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
    int i = 0;
    for (auto& point : gPoints) {
        for (i = 0; i < 4; ++i) {
            point->addSupportNode(dynamic_cast<Node *>(bodyPoints[i]));
        }
    }
    for (auto& point : gPoints_MC) {
        for (i = 0; i < 4; ++i) {
            point->addSupportNode(dynamic_cast<Node *>(this->bodyPoints[i]));
        }
    }
}


void ElemTetrahedron::computeShapeFunctions()
{
    for (auto& point : gPoints) {
        point->fillFEmatrices();
    }
    for (auto& point : gPoints_MC) {
        point->fillFEmatrices();
    }
}

void ElemTetrahedron::createGaussPoints_MC()
{
    int nGPoints_MC = 4;
    if (nGPoints == 1) {
        nGPoints_MC = 4;
    } else if (nGPoints == 4) nGPoints_MC = 4; // TODO: More needed

    lmx::DenseMatrix<double> gCoef(size_type(nGPoints_MC), 5);

    // reference: http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
    if (nGPoints_MC == 1) {
        gCoef(0, 0) = .25;
        gCoef(0, 1) = .25;
        gCoef(0, 2) = .25;
        gCoef(0, 3) = .25;
        gCoef(0, 4) = 1.;
    }
    else if (nGPoints_MC == 4) {
        gCoef(0, 0) = .585410196624969;
        gCoef(0, 1) = gCoef(0, 2) = gCoef(0, 3) = .138196601125011;
        gCoef(0, 4) = .25;
        gCoef(1, 1) = .585410196624969;
        gCoef(1, 0) = gCoef(1, 2) = gCoef(1, 3) = .138196601125011;
        gCoef(1, 4) = .25;
        gCoef(2, 2) = .585410196624969;
        gCoef(2, 1) = gCoef(2, 0) = gCoef(2, 3) = .138196601125011;
        gCoef(2, 4) = .25;
        gCoef(3, 3) = .585410196624969;
        gCoef(3, 1) = gCoef(3, 2) = gCoef(3, 0) = .138196601125011;
        gCoef(3, 4) = .25;
    }
// more on ... http://electromagnetics.biz/2D%20Gauss.txt
//   else if(nGPoints == 12){
//   }
//   else if(nGPoints == 20){
//   }

    for (int i = 0; i < nGPoints_MC; ++i) {
        gPoints_MC.push_back
                (new GaussPoint3D
                         (this->alpha,
                          gCoef(i, 4),
                          jacobian,
                          mat,
                          i,
                          bodyPoints[0]->getX() * gCoef(i, 0)
                          + bodyPoints[1]->getX() * gCoef(i, 1)
                          + bodyPoints[2]->getX() * gCoef(i, 2)
                          + bodyPoints[3]->getX() * gCoef(i, 3),
                          bodyPoints[0]->getY() * gCoef(i, 0)
                          + bodyPoints[1]->getY() * gCoef(i, 1)
                          + bodyPoints[2]->getY() * gCoef(i, 2)
                          + bodyPoints[3]->getY() * gCoef(i, 3),
                          bodyPoints[0]->getZ() * gCoef(i, 0)
                          + bodyPoints[1]->getZ() * gCoef(i, 1)
                          + bodyPoints[2]->getZ() * gCoef(i, 2)
                          + bodyPoints[3]->getZ() * gCoef(i, 3),
                          dc,
                          false
                         )
                );
    }
}


}
