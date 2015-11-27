//-- Licencia --
#ifndef CELLTRIANG_H
#define CELLTRIANG_H

#include "LMX/cofe_TensorRank2.h"
#include "cell.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file celltriang.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

/**
@author Daniel Iglesias
*/
class CellTriang : public Cell {
protected:
    cofe::TensorRank2<3,double> points; /**< position of vertex points */

public:
    CellTriang();

    CellTriang( Material&,
                std::string,
                double,
		int,
                Point*,
                Point*,
                Point*
              );

    CellTriang( Material&,
                std::string,
                double,
		int,
                Point*,
                Point*,
                Point*,
		double
              );

    virtual ~CellTriang();

    void gnuplotOut( std::ofstream&, std::ofstream& );

protected:
    void createGaussPoints( );


};

} //Namespace mknix

#endif
