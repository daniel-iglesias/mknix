//-- Licencia --
#ifndef MKNIXGAUSSPOINT_H
#define MKNIXGAUSSPOINT_H

#include "LMX/lmx.h"
#include "point.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file gausspoint.h

  \brief Point for numerical integration.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

class Material;
class LoadThermalBody;

/**
@author Daniel Iglesias
*/
class GaussPoint : public Point {

public:
    GaussPoint();

    GaussPoint( int dim_in, double alpha_in, double weight_in, double jacobian_in,
                Material* mat_in, int num, double coor_x, double coor_y,
                double dc_in, bool );

    GaussPoint( int dim_in, double alpha_in, double weight_in, double jacobian_in,
                Material* mat_in, int num, double coor_x, double coor_y,
                double coor_z, double dc_in, bool );

    virtual ~GaussPoint();

    virtual void shapeFunSolve( std::string, double ) override;

    virtual void fillFEmatrices( )=0;

    void computeCij( );

    void computeHij( );

    void computeQext( LoadThermalBody* );

    virtual void computeFint( )=0;

    virtual void computeFext( )=0;

    virtual void computeMij( )=0;

    virtual void computeKij( )=0;

    virtual void computeStress( )=0;

    virtual void computeNLStress( )=0;

    virtual void computeNLFint( )=0;

    virtual void computeNLKij( )=0;

    void assembleCij( lmx::Matrix<data_type> & );

    void assembleHij( lmx::Matrix<data_type> & );

    void assembleQext( lmx::Vector<data_type> & );

    virtual void assembleMij( lmx::Matrix<data_type> & )=0;

    virtual void assembleKij( lmx::Matrix<data_type> & )=0;

    virtual void assembleRi( lmx::Vector<data_type> &, int )=0;

    virtual void assembleFint( lmx::Vector<data_type> & )=0;

    virtual void assembleFext( lmx::Vector<data_type> & )=0;

    virtual double calcPotentialE( const lmx::Vector<data_type> & )=0;

    virtual double calcKineticE( const lmx::Vector<data_type> & )=0;

    virtual double calcElasticE( )=0;

    void gnuplotOutStress( std::ofstream& );

protected:
    int num;
    double weight;
    Material* mat;
    bool stressPoint;

    lmx::DenseMatrix<data_type> B;
    lmx::DenseMatrix<data_type> C;
    lmx::DenseMatrix<data_type> H;
    lmx::DenseMatrix<data_type> M;
    lmx::DenseMatrix<data_type> K;

    lmx::Vector<data_type> tension;
    lmx::Vector<data_type> r; // = integral( Phi^T * tension )dA

    lmx::Vector<data_type> Qext;
    lmx::Vector<data_type> fint;
    lmx::Vector<data_type> fext;
};

} //Namespace mknix

#endif
