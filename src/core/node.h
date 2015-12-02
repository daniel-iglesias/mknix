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

#ifndef MKNIXNODE_H
#define MKNIXNODE_H

#include "point.h"
#include "LMX/lmx.h"

namespace mknix {

/**
  @author AUTHORS <MAILS>
*/
class Node : public Point {
public:
    Node();

    Node(const Node&);

    Node(const Node*);

    Node(int i_in, double x_in, double y_in, double z_in);

    ~Node();

    void addWeight( double w_in ) {
        this->weight += w_in;
    }

    inline const double& getWeight() const {
        return this->weight;
    }

    inline const int& getThermalNumber() const {
        return this->thermalNum;
    }

    inline const void setNumber( int num_in ) {
        num = num_in;
    }

    inline const void setThermalNumber( int num_in ) {
        thermalNum = num_in;
    }

//     inline const double& getx() const {
//         return this->x;
//     }
// 
//     inline const double& gety() const {
//         return this->y;
//     }
// 
//     inline const double& getz() const {
//         return this->z;
//     }

    inline const double& getqt() const {
        return this->qt;
    }

    inline double getqx(int i) const
    {
        if(i==0) return qx;
        if(i==1) return qy;
        if(i==2) return qz;
        else return 1E10;
    }

    inline double getUx() const {
        return getConf(0)-X;
    }

    inline double getUy() const {
        return getConf(1)-Y;
    }

    inline double getUz() const {
        return getConf(2)-Z;
    }

    inline double getU(int i) const
    {
        if(i==0) return getConf(0)-X;
        if(i==1) return getConf(1)-Y;
        if(i==2) return getConf(2)-Z;
        else return 1E10;
    }

    double getConf(int dof) const override;

    double getTemp() const override;

    size_t getSupportSize( int deriv );

    int getSupportNodeNumber( int deriv, int s_node );

    double getShapeFunValue( int deriv, int s_node );

//  TODO: Check the following functions:
    inline void setUx(double ux_in) {
        qx = X + ux_in;
    }

    inline void setUy(double uy_in) {
        qy = Y + uy_in;
    }

    inline void setUz(double uz_in) {
        qz = Z + uz_in;
    }

    void setqx(const lmx::Vector<data_type>& globalConf, int dim);

    void setqt( const lmx::Vector<data_type>& globalConf );

    void setqt( double& temp_in ){
      qt = temp_in;
    }
//     inline void setx(double x_in) {
//         qx = x_in;
//     }
// 
//     inline void sety(double y_in) {
//         y = y_in;
//     }
// 
//     inline void setz(double z_in) {
//         z = z_in;
//     }
// 
    inline void setX(double X_in) {
        X = X_in;
        qx = X_in;
    }

    inline void setY(double Y_in) {
        Y = Y_in;
        qy = Y_in;
    }

    inline void setZ(double Z_in) {
        Z = Z_in;
        qz = Z_in;
    }

private:
    int thermalNum;
    double qx, qy, qz, qt;
//    double X, Y, Z;
    double weight;

////// Members possibly needed for additional
//     std::string shapeFunTypeThermal;
//     std::vector<Node*> supportNodesThermal;
//     ShapeFunction* shapeFunThermal;
//     int supportNodesSizeThermal;
};

}

#endif
