//-- Licencia --
#ifndef MKNIXGAUSSPOINT3D_H
#define MKNIXGAUSSPOINT3D_H

#include "LMX/lmx.h"
#include "gausspoint.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file gausspoint.3Dh

  \brief Point for numerical integration in 3D cells.

  \author Daniel Iglesias

 */
//////////////////////////////////////////// Doxygen file documentation (end)

namespace mknix {

class Material;
class LoadThermal;

/**
@author Daniel Iglesias
*/
class GaussPoint3D : public GaussPoint {

public:
    GaussPoint3D();

    GaussPoint3D( double alpha_in, double weight_in, double jacobian_in,
                  Material* mat_in, int num, double coor_x, double coor_y,
                  double coor_z, double dc_in, bool stressPoint_in );

    ~GaussPoint3D();

    void shapeFunSolve( std::string, double );

    void fillFEmatrices( );

    void computeMij( );

    void computeKij( );

    void computeStress( );

    void computeNLStress( );

    void computeFint( );

    void computeFext( );

    void computeNLFint( );

    void computeNLKij( );

    void assembleMij( lmx::Matrix<data_type> & );

    void assembleKij( lmx::Matrix<data_type> & );

    void assembleRi( lmx::Vector<data_type> &, int );

    void assembleFint( lmx::Vector<data_type> & );

    void assembleFext( lmx::Vector<data_type> & );

    double calcPotentialE( const lmx::Vector<data_type> & );

    double calcKineticE( const lmx::Vector<data_type> & );

    double calcElasticE( );

private:
    void initializeMatVecs();

private:
    cofe::TensorRank2<3,double> F3;
    cofe::TensorRank2Sym<3,double> S3;
    cofe::TensorRank2<3,double> P3;
    cofe::TensorRank2Sym<3,double> sigma3;
};

} //Namespace mknix

#endif
