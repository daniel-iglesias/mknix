//-- Licencia --
#ifndef CELLBOUNDARYLINEAR_H
#define CELLBOUNDARYLINEAR_H

#include "LMX/cofe_TensorRank2.h"
#include "cellboundary.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file cellboundarylinear.h

  \brief Background cells for boundary integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

/**
@author Daniel Iglesias
*/
class CellBoundaryLinear : public CellBoundary {
protected:
    cofe::TensorRank2<3,double> points; /**< position of vertex points */

public:
    CellBoundaryLinear();

    CellBoundaryLinear( 
                std::string,
                double,
		int,
                Point*,
                Point*
              );

    CellBoundaryLinear( 
                std::string,
                double,
		int,
                Point*,
                Point*,
		double
              );

    ~CellBoundaryLinear();

    virtual void initialize( std::vector<Node*> & );

    virtual void computeShapeFunctions(  );

//     void gnuplotOut( std::ofstream&, std::ofstream& );

protected:
    void createGaussPoints( );


};

} //Namespace mknix

#endif
