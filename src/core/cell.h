/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of Nemesis.                                             *
 *                                                                            *
 *  Nemesis is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  Nemesis is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with Nemesis.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#ifndef CELL_H
#define CELL_H

#include <vector>
#include <string>
#include <common.h>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file cell.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {
template <typename T> class Vector;
template <typename T> class Matrix;
}

namespace mknix {

class LoadThermalBody;
class Material;
class GaussPoint;
class Node;
class Point;

/**
@author Daniel Iglesias
*/
class Cell {

public:
    Cell();

    Cell( Material&, std::string, double, int );

    virtual ~Cell();

    int getMaterialId();

    bool setMaterialIfLayer( Material&, double );

    virtual void initialize( std::vector<Node*> & );

    virtual void computeShapeFunctions(  );

    void computeCapacityGaussPoints(  );

    void assembleCapacityGaussPoints( lmx::Matrix<data_type> & );

    void assembleCapacityGaussPointsWithMap(data_type *globalCapacity,
                                            int *matrixMap,
                                            int number_nodes);

    void presenceCapacityGaussPoints(int* presence_matrix, int number_nodes);

    void computeConductivityGaussPoints(  );

    void assembleConductivityGaussPoints( lmx::Matrix<data_type> & );

    void assembleConductivityGaussPointsWithMap(data_type *globalConductivity,
                                                int *matrixMap,
                                                int number_nodes);

    void computeQextGaussPoints( LoadThermalBody* );

    void assembleQextGaussPoints( lmx::Vector<data_type> & );

    void computeMGaussPoints(  );

    void assembleMGaussPoints( lmx::Matrix<data_type> & );

    void computeFintGaussPoints(  );

    void computeNLFintGaussPoints(  );

    void assembleFintGaussPoints( lmx::Vector<data_type> & );

    void computeFextGaussPoints(  );

    void assembleFextGaussPoints( lmx::Vector<data_type> & );

    void computeKGaussPoints(  );

    void computeNLKGaussPoints(  );

    void assembleKGaussPoints( lmx::Matrix<data_type> & );

    void assembleRGaussPoints( lmx::Vector<data_type> &, int );

    void assembleNLRGaussPoints( lmx::Vector<data_type> &, int );

    double calcPotentialEGaussPoints( const lmx::Vector<data_type> & );

    double calcKineticEGaussPoints( const lmx::Vector<data_type> & );

    double calcElasticEGaussPoints(  );

    void outputConnectivityToFile(std::ofstream*);

    virtual void gnuplotOut( std::ofstream&, std::ofstream& ) = 0;

    void gnuplotOutStress( std::ofstream& );

protected:
    Material* mat;
    std::string formulation;
    double alpha;
    int nGPoints; /**< number of Gauss Points */
    std::vector<GaussPoint*> gPoints;
    std::vector<GaussPoint*> gPoints_MC;
    double jacobian;
    std::vector< Point* > bodyPoints;
    double dc;

};

} //Namespace mknix

#endif
