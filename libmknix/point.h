//-- Licencia --
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
class Point {

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
    
    Point( const Point& point_in);

    Point( const Point* point_in);

    Point( int i, double coor_x, double coor_y, double coor_z );

    Point( int dim, int i, double coor_x, double coor_y, double coor_z,
           double alpha_in, double dc_in );

    Point( int dim, int i, double coor_x, double coor_y, double coor_z,
           double jacobian_in, double alpha_in, double dc_in );

    ~Point();

    inline const int& getDim() const {
        return this->dim;
    }

    inline const double& getX() const {
        return this->X;
    }

    inline const double& getY() const {
        return this->Y;
    }

    inline const double& getZ() const {
        return this->Z;
    }

    inline const int& getNumber() const {
        return this->num;
    }

    virtual double getTemp();

    virtual double getConf(int);
    
    double distance( Point& );

    void setShapeFunType( std::string type_in )
    {
        shapeFunType = type_in;
    }

    std::string getShapeFunType( )
    {
        return shapeFunType;
    }

    virtual size_t getSupportSize( int deriv=0 )
    {
        return supportNodesSize;
    }

    virtual int getSupportNodeNumber( int deriv, int s_node );

    void setDc( double dc_in )
    {   // set dc only if is a greater (conservative) value:
        if (dc_in > dc) dc = dc_in;
    }

    void setAlphai( double alphai_in )
    {
        alphai = alphai_in;
    }

    void addSupportNode( Node* );

    void findSupportNodes ( std::vector< Node* > & );

    // for rectangular domains
    void findSupportNodes ( std::vector< Node* > &,
                            double, double, double, double );

    void shapeFunSolve( std::string, double );

    void shapeFunSolve( double q_in )
    { shapeFunSolve( shapeFunType, q_in ); }

    void setJacobian(double jac_in)
    {
        jacobian = jac_in;
    }

    void gnuplotOut( std::ofstream& );

protected:
    int dim, num;
    double X, Y, Z;
    double alphai;
    double dc;
    std::string shapeFunType;
    std::vector<Node*> supportNodes;
    ShapeFunction* shapeFun;
    int supportNodesSize;
    double jacobian;


};

} //Namespace mknix

#endif
