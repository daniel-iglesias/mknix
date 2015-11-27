//-- Licencia --
#ifndef CELLRECT_H
#define CELLRECT_H

#include "LMX/lmx_mat_dense_matrix.h"
#include "cell.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file cellrect.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace lmx {
template <typename T> class DenseMatrix;
}

namespace mknix {

/**
@author Daniel Iglesias
*/
class CellRect : public Cell {
private:
    double Ax,Ay,Az;
    double minX, minY, minZ, maxX, maxY, maxZ;
    lmx::DenseMatrix<double> points; /**< position of vertex points */

public:
    CellRect();

    CellRect( Material&,
              std::string,
              double, int,
              double, double,
              double, double,
              double, double,
              double, double,
              double, double,
              double, double,
              double, double );

    CellRect( Material&,
              std::string,
              double, int,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double,
              double, double, double
            );

    ~CellRect();

//    void initialize( std::map<int,Node*> & );

    void initialize( std::vector<Node*> & );

    void gnuplotOut( std::ofstream&, std::ofstream& );

private:
    void createGaussPoints( double, double, double = 0 );


};

} //Namespace mknix

#endif
