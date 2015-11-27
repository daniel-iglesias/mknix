//-- Licencia --
#include "celltriang.h"
#include "material.h"
#include "node.h"
#include "gausspoint2D.h"

namespace mknix {

CellTriang::CellTriang()
{
}


CellTriang::CellTriang(Material& material_in,
                       std::string formulation_in,
                       double alpha_in,
                       int nGPoints_in,
                       Point * n1_in,
                       Point * n2_in,
                       Point * n3_in
)
        : Cell(material_in, formulation_in, alpha_in, nGPoints_in
)
{
    this->bodyPoints.push_back(n1_in);
    this->bodyPoints.push_back(n2_in);
    this->bodyPoints.push_back(n3_in);

    points(0, 0) = n1_in->getX();
    points(0, 1) = n1_in->getY();
    points(0, 2) = 1.;
    points(1, 0) = n2_in->getX();
    points(1, 1) = n2_in->getY();
    points(1, 2) = 1.;
    points(2, 0) = n3_in->getX();
    points(2, 1) = n3_in->getY();
    points(2, 2) = 1.;

    this->jacobian = 0.5 * points.determinant();

    dc = (n1_in->distance(*n2_in)
          + n1_in->distance(*n3_in)
          + n2_in->distance(*n3_in)) / 3.;

    this->createGaussPoints();
}

CellTriang::CellTriang(Material& material_in,
                       std::string formulation_in,
                       double alpha_in,
                       int nGPoints_in,
                       Point * n1_in,
                       Point * n2_in,
                       Point * n3_in,
                       double dc_in
)
        : Cell(material_in, formulation_in, alpha_in, nGPoints_in
)
{
    this->bodyPoints.push_back(n1_in);
    this->bodyPoints.push_back(n2_in);
    this->bodyPoints.push_back(n3_in);

    points(0, 0) = n1_in->getX();
    points(0, 1) = n1_in->getY();
    points(0, 2) = 1.;
    points(1, 0) = n2_in->getX();
    points(1, 1) = n2_in->getY();
    points(1, 2) = 1.;
    points(2, 0) = n3_in->getX();
    points(2, 1) = n3_in->getY();
    points(2, 2) = 1.;

    this->jacobian = 0.5 * points.determinant();

    dc = dc_in;

    this->createGaussPoints();
}

CellTriang::~CellTriang()
{
}


void CellTriang::createGaussPoints()
{
    lmx::DenseMatrix<double> gCoef(size_type(nGPoints), 4);

    // reference: http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
    if (nGPoints == 1) {
        gCoef(0, 0) = 1. / 3.;
        gCoef(0, 1) = 1. / 3.;
        gCoef(0, 2) = 1. / 3.;
        gCoef(0, 3) = 1.;
    }
    else if (nGPoints == 3) {
        gCoef(0, 0) = 2. / 3.;
        gCoef(0, 1) = 1. / 6.;
        gCoef(0, 2) = 1. / 6.;
        gCoef(0, 3) = 1. / 3.;
        gCoef(1, 0) = 1. / 6.;
        gCoef(1, 1) = 2. / 3.;
        gCoef(1, 2) = 1. / 6.;
        gCoef(1, 3) = 1. / 3.;
        gCoef(2, 0) = 1. / 6.;
        gCoef(2, 1) = 1. / 6.;
        gCoef(2, 2) = 2. / 3.;
        gCoef(2, 3) = 1. / 3.;
// From Zienkiewicz, Taylor:
//         gCoef(0,0) = 0.5;
//         gCoef(0,1) = 0.5;
//         gCoef(0,2) = 0. ;
//         gCoef(0,3) = 1./3.;
//         gCoef(1,0) = 0. ;
//         gCoef(1,1) = 0.5;
//         gCoef(1,2) = 0.5;
//         gCoef(1,3) = 1./3.;
//         gCoef(2,0) = 0.5;
//         gCoef(2,1) = 0. ;
//         gCoef(2,2) = 0.5;
//         gCoef(2,3) = 1./3.;
    }
    else if (nGPoints == 4) {
        gCoef(0, 0) = 1. / 3.;
        gCoef(0, 1) = 1. / 3.;
        gCoef(0, 2) = 1. / 3.;
        gCoef(0, 3) = -27. / 48.;
        gCoef(1, 0) = 0.6;
        gCoef(1, 1) = 0.2;
        gCoef(1, 2) = 0.2;
        gCoef(1, 3) = 25. / 48.;
        gCoef(2, 0) = 0.2;
        gCoef(2, 1) = 0.6;
        gCoef(2, 2) = 0.2;
        gCoef(2, 3) = 25. / 48.;
        gCoef(3, 0) = 0.2;
        gCoef(3, 1) = 0.2;
        gCoef(3, 2) = 0.6;
        gCoef(3, 3) = 25. / 48.;
    }
    else if (nGPoints == 6) { //must be checked
        double alpha1 = 0.8168475730;
        double alpha2 = 0.0915762135;
        double beta1 = 0.1081030182;
        double beta2 = 0.4459484909;
        gCoef(3, 0) = beta1;
        gCoef(3, 1) = beta1;
        gCoef(3, 2) = alpha1;
        gCoef(3, 3) = 0.0549758718;
        gCoef(1, 0) = alpha1;
        gCoef(1, 1) = beta1;
        gCoef(1, 2) = beta1;
        gCoef(1, 3) = 0.0549758718;
        gCoef(2, 0) = beta1;
        gCoef(2, 1) = alpha1;
        gCoef(2, 2) = beta1;
        gCoef(2, 3) = 0.0549758718;
        gCoef(3, 0) = beta2;
        gCoef(3, 1) = beta2;
        gCoef(3, 2) = alpha2;
        gCoef(3, 3) = 0.1116907998;
        gCoef(4, 0) = alpha2;
        gCoef(4, 1) = beta2;
        gCoef(4, 2) = beta2;
        gCoef(4, 3) = 0.1116907998;
        gCoef(5, 0) = beta2;
        gCoef(5, 1) = alpha2;
        gCoef(5, 2) = beta2;
        gCoef(5, 3) = 0.1116907998;
    }
    else if (nGPoints == 7) {
        double alpha1 = 0.0597158717;
        double alpha2 = 0.7974269853;
        double beta1 = 0.4701420641;
        double beta2 = 0.1012865073;
        gCoef(0, 0) = 1. / 3.;
        gCoef(0, 1) = 1. / 3.;
        gCoef(0, 2) = 1. / 3.;
        gCoef(0, 3) = 0.225;
        gCoef(1, 0) = alpha1;
        gCoef(1, 1) = beta1;
        gCoef(1, 2) = beta1;
        gCoef(1, 3) = 0.1323941527;
        gCoef(2, 0) = beta1;
        gCoef(2, 1) = alpha1;
        gCoef(2, 2) = beta1;
        gCoef(2, 3) = 0.1323941527;
        gCoef(3, 0) = beta1;
        gCoef(3, 1) = beta1;
        gCoef(3, 2) = alpha1;
        gCoef(3, 3) = 0.1323941527;
        gCoef(4, 0) = alpha2;
        gCoef(4, 1) = beta2;
        gCoef(4, 2) = beta2;
        gCoef(4, 3) = 0.1259391805;
        gCoef(5, 0) = beta2;
        gCoef(5, 1) = alpha2;
        gCoef(5, 2) = beta2;
        gCoef(5, 3) = 0.1259391805;
        gCoef(6, 0) = beta2;
        gCoef(6, 1) = beta2;
        gCoef(6, 2) = alpha2;
        gCoef(6, 3) = 0.1259391805;
    }
    else if (nGPoints == 12) {
        double alpha1 = 0.873821971016996;
        double beta1 = 0.063089014491502;
        double weight1 = 0.05084490637027;
        double alpha2 = 0.501426509658179;
        double beta2 = 0.249286745170910;
        double weight2 = 0.116786275726379;
        //double alpha3 = 0.6365024969121399;
        //double beta3 = 0.310352451033785;
        double gamma3 = 0.053145049844816;
        double weight3 = 0.082851075618374;
        gCoef(0, 0) = alpha1;
        gCoef(0, 1) = beta1;
        gCoef(0, 2) = beta1;
        gCoef(0, 3) = weight1;
        gCoef(1, 0) = beta1;
        gCoef(1, 1) = alpha1;
        gCoef(1, 2) = beta1;
        gCoef(1, 3) = weight1;
        gCoef(2, 0) = beta1;
        gCoef(2, 1) = beta1;
        gCoef(2, 2) = alpha1;
        gCoef(2, 3) = weight1;
        gCoef(3, 0) = alpha2;
        gCoef(3, 1) = beta2;
        gCoef(3, 2) = beta2;
        gCoef(3, 3) = weight2;
        gCoef(4, 0) = beta2;
        gCoef(4, 1) = alpha2;
        gCoef(4, 2) = beta2;
        gCoef(4, 3) = weight2;
        gCoef(5, 0) = beta2;
        gCoef(5, 1) = beta2;
        gCoef(5, 2) = alpha2;
        gCoef(5, 3) = weight2;
        gCoef(6, 0) = alpha1;
        gCoef(6, 1) = beta1;
        gCoef(6, 2) = gamma3;
        gCoef(6, 3) = weight3;
        gCoef(7, 0) = gamma3;
        gCoef(7, 1) = alpha1;
        gCoef(7, 2) = beta1;
        gCoef(7, 3) = weight3;
        gCoef(8, 0) = beta1;
        gCoef(8, 1) = gamma3;
        gCoef(8, 2) = alpha1;
        gCoef(8, 3) = weight3;
        gCoef(9, 0) = alpha1;
        gCoef(9, 1) = gamma3;
        gCoef(9, 2) = beta1;
        gCoef(9, 3) = weight3;
        gCoef(10, 0) = beta1;
        gCoef(10, 1) = alpha1;
        gCoef(10, 2) = gamma3;
        gCoef(10, 3) = weight3;
        gCoef(11, 0) = gamma3;
        gCoef(11, 1) = beta1;
        gCoef(11, 2) = alpha1;
        gCoef(11, 3) = weight3;

    }
//   else if(nGPoints == 20){
//   }

    for (int i = 0; i < nGPoints; ++i) {
        gPoints.push_back(new GaussPoint2D(this->alpha,
                                           gCoef(i, 3),
                                           jacobian,
                                           mat,
                                           i,
                                           points(0, 0) * gCoef(i, 0) + points(1, 0) * gCoef(i, 1) +
                                           points(2, 0) * gCoef(i, 2),
                                           points(0, 1) * gCoef(i, 0) + points(1, 1) * gCoef(i, 1) +
                                           points(2, 1) * gCoef(i, 2),
                                           dc,
                                           true
        )
        );
    }
}


void CellTriang::gnuplotOut(std::ofstream& data, std::ofstream& gpdata)
{
    for (int i = 0; i < 3; ++i) {
        data << points(i, 0) << " " << points(i, 1) << " 0\n";
    }
    data << points(0, 0) << " " << points(0, 1) << " 0\n";
    data << std::endl;

    for (auto& point : gPoints) {
        point->gnuplotOut(gpdata);
    }
}

} //Namespace mknix
