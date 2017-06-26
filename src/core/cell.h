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
#include <gpu/cpu_run_type.h>

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

    double getCondPhi(int gp, int deriv, int node);

    double getCapPhi(int gp, int deriv, int node);

    double getWeight(int gp);

    double getWeightMC(int gp);

    double getJacobian();

    double getJacobianP(int gp);

    double getJacobianMC(int gp);

    int getSupportSizeMC();
    int getSupportSize();

    bool setMaterialIfLayer( Material&, double );

    virtual void initialize( std::vector<Node*> & );

    virtual void computeShapeFunctions(  );

    void computeCapacityGaussPoints(  );

    std::vector<double> getShapeCij();
    std::vector<double> getCij();
    std::vector<double> getTempsCij();//FOR DEBUG ONLY
    double getCFactor();

    int getNumPoints_MC();
    int getNumPoints();

    void assembleCapacityGaussPoints( lmx::Matrix<data_type> & );

    void assembleCapacityGaussPointsWithMap(data_type *globalCapacity,
                                            uint *matrixMap,
                                            int number_nodes);

    void presenceCapacityGaussPoints(int* presence_matrix, int number_nodes);

    void computeConductivityGaussPoints(  );

    std::vector<double> getShapeHij();
    std::vector<double> getHij();
    std::vector<double> getTempsHij();//FOR DEBUG ONLY
    double getHFactor();

    void assembleConductivityGaussPoints( lmx::Matrix<data_type> & );

    void assembleConductivityGaussPointsWithMap(data_type *globalConductivity,
                                                uint *matrixMap,
                                                int number_nodes);

    void presenceConductivityGaussPoints(int* presence_matrix, int number_nodes);

    void computeQextGaussPoints( LoadThermalBody* );

    void assembleQextGaussPoints( lmx::Vector<data_type> & );
    void assembleQextGaussPoints( VectorX<data_type> & );

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

    double calcPotentialEGaussPoints( const VectorX<data_type> & );

    double calcKineticEGaussPoints( const lmx::Vector<data_type> & );

    double calcElasticEGaussPoints(  );

    void outputConnectivityToFile(std::ofstream*);

    virtual void gnuplotOut( std::ofstream&, std::ofstream& ) = 0;

    void gnuplotOutStress( std::ofstream& );

    void mapThermalNodesMC(int* thermalMap,
                           int supportNodeSize,
                           int cell_index);

    void mapThermalNodes(int* thermalMap,
                         int supportNodeSize,
                         int cell_index);
    //
    void mapNodesMC(uint* Map,
                    int supportNodeSize,
                    int cell_index,
                    int total_nodes);
    //
    void mapNodes(uint* Map,
                  int supportNodeSize,
                  int cell_index,
                  int total_nodes);
    //

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
