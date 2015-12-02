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

#ifndef POINT_H
#define POINT_H

#include <vector>
#include <map>
#include <string>

#include "common.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file point.h

  \brief Point of interest.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

class Node;

class ShapeFunction;

class ShapeFunctionRBF;

class Simulation;


/**
@author Daniel Iglesias
*/
class Point
{

public:
    // TODO: This frind class stuff should be removed. It is dangerous!
    friend class GaussPoint;

    friend class ShapeFunction;

    friend class ShapeFunctionRBF;

    friend class ShapeFunctionMLS;

    friend class ShapeFunctionTetrahedron;

    friend class ShapeFunctionTriangle;

    friend class ShapeFunctionLinear;

    friend class ShapeFunctionLinearX;

public:
    Point();

    Point(const Point& point_in);

    Point(const Point * point_in);

    Point(int i, double coor_x, double coor_y, double coor_z);

    Point(int dim, int i, double coor_x, double coor_y, double coor_z,
          double alpha_in, double dc_in);

    Point(int dim, int i, double coor_x, double coor_y, double coor_z,
          double jacobian_in, double alpha_in, double dc_in);

    virtual ~Point();

    inline const int& getDim() const
    {
        return this->dim;
    }

    inline const double& getX() const
    {
        return this->X;
    }

    inline const double& getY() const
    {
        return this->Y;
    }

    inline const double& getZ() const
    {
        return this->Z;
    }

    inline const int& getNumber() const
    {
        return this->num;
    }

    virtual double getTemp() const;

    virtual double getConf(int) const;

    double distance(Point&) const;

    void setShapeFunType(std::string type_in)
    {
        shapeFunType = type_in;
    }

    std::string getShapeFunType()
    {
        return shapeFunType;
    }

    virtual size_t getSupportSize(int deriv = 0)
    {
        return supportNodesSize;
    }

    virtual int getSupportNodeNumber(int deriv, int s_node);

    void setDc(double dc_in)
    {   // set dc only if is a greater (conservative) value:
        if (dc_in > dc) dc = dc_in;
    }

    void setAlphai(double alphai_in)
    {
        alphai = alphai_in;
    }

    void addSupportNode(Node *);

    void findSupportNodes(std::vector<Node *>&);

    // for rectangular domains
    void findSupportNodes(std::vector<Node *>&,
                          double, double, double, double);

    virtual void shapeFunSolve(std::string, double);

    void shapeFunSolve(double q_in) { shapeFunSolve(shapeFunType, q_in); }

    void setJacobian(double jac_in)
    {
        jacobian = jac_in;
    }

    void gnuplotOut(std::ofstream&);

    const std::vector<Node *>& getSupportNodes() const
    {
        return supportNodes;
    }

protected:
    int dim;
    int num;
    double X, Y, Z;
    double alphai;
    double dc;
    std::string shapeFunType;
    ShapeFunction * shapeFun;
    std::vector<Node *> supportNodes;
    size_t supportNodesSize;
    double jacobian;
};

} //Namespace mknix

#endif
