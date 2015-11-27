//-- Licencia --
#ifndef CELLTETRAHEDRON_H
#define CELLTETRAHEDRON_H

#include "LMX/cofe_TensorRank2.h"
#include "cell.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file celltetrahedron.h

  \brief Background cells for integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)


namespace mknix {

/**
@author Daniel Iglesias
*/
class CellTetrahedron : public Cell {
protected:
    cofe::TensorRank2<3,double> points; /**< position of vertex points */

public:
    CellTetrahedron();

    CellTetrahedron( Material&,
                     std::string,
                     double, int,
                     Point*,
                     Point*,
                     Point*,
                     Point*
                   );

    virtual ~CellTetrahedron();

    void gnuplotOut( std::ofstream&, std::ofstream& );

protected:
    void createGaussPoints( );


};

} //Namespace mknix

#endif
