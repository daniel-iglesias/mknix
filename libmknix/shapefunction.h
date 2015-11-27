//-- Licencia --
#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H

#include "LMX/lmx.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file shapefunction.h

  \brief Function for interpolation of basic variables.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

class Point;

/**
@author Daniel Iglesias
*/
class ShapeFunction {

protected:
    int dim;
    lmx::DenseMatrix<double> phi;
    Point* gp;

public:
    ShapeFunction();

    ShapeFunction( const ShapeFunction* );

    ShapeFunction( Point* );

    virtual ~ShapeFunction();

    virtual void calc() = 0;

    virtual double getPhi(int i, int j) // i = derivative order, j = node
    {
        return phi.readElement(i, j);
    }

    virtual void setPhi(double value_in, int i, int j) // i = derivative order, j = node
    {
        phi.writeElement(value_in, i, j);
    }

    virtual void outputValues();

    virtual void gnuplotOut();


};

} //Namespace mknix

#endif
