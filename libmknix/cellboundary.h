//-- Licencia --
#ifndef CELLBOUNDARY_H
#define CELLBOUNDARY_H

#include <vector>
#include <string>
#include "common.h"

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

class LoadThermalBoundary1D;
class GaussPointBoundary;
class Node;
class Point;

/**
@author Daniel Iglesias
*/
class CellBoundary {

public:
    CellBoundary();

    CellBoundary( std::string, double, int );

    virtual ~CellBoundary();

    virtual void initialize( std::vector<Node*> & );

    virtual void computeShapeFunctions(  );

    void computeQextGaussPoints( LoadThermalBoundary1D* );

    void assembleQextGaussPoints( lmx::Vector<data_type> & );

//     void outputConnectivityToFile(std::ofstream*);
    
//     virtual void gnuplotOut( std::ofstream&, std::ofstream& ) = 0;

//     void gnuplotOutStress( std::ofstream& );

protected:
    std::string formulation;
    double alpha;
    int nGPoints; /**< number of Gauss Points */
    std::vector<GaussPointBoundary*> gPoints;
    double jacobian;
    std::vector< Point* > bodyPoints;
    double dc;

};

} //Namespace mknix

#endif
