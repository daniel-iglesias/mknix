//-- Licencia --
#include "cellrect.h"
#include "material.h"
#include "node.h"
#include "gausspoint2D.h"
#include <string>

namespace mknix {

CellRect::CellRect()
{
}


CellRect::CellRect( Material& material_in,
                    std::string formulation_in,
                    double alpha_in,
                    int nGPoints_in,
                    double x1_in, double y1_in,
                    double x2_in, double y2_in,
                    double x3_in, double y3_in,
                    double x4_in, double y4_in,
                    double dcx_in,double dcy_in,
                    double minX_in,double minY_in,
                    double maxX_in,double maxY_in
                  )
    : Cell( material_in
            , formulation_in
            , alpha_in
            , nGPoints_in
          )
    , Az(0)
    , minX(minX_in)
    , minY(minY_in)
    , minZ(0)
    , maxX(maxX_in)
    , maxY(maxY_in)
    , maxZ(0)
{
    points.resize(4,2);
    points.writeElement(x1_in,0,0);
    points(0,1) = y1_in;
    points(1,0) = x2_in;
    points(1,1) = y2_in;
    points(2,0) = x3_in;
    points(2,1) = y3_in;
    points(3,0) = x4_in;
    points(3,1) = y4_in;

    Ax = x2_in - x1_in;
    Ay = y3_in - y2_in;

    this->jacobian = Ax/2. * Ay/2.;

    this->createGaussPoints( dcx_in, dcy_in );

}

CellRect::CellRect( Material& material_in,
                    std::string formulation_in,
                    double alpha_in,
                    int nGPoints_in,
                    double x1_in, double y1_in, double z1_in,
                    double x2_in, double y2_in, double z2_in,
                    double x3_in, double y3_in, double z3_in,
                    double x4_in, double y4_in, double z4_in,
                    double x5_in, double y5_in, double z5_in,
                    double x6_in, double y6_in, double z6_in,
                    double x7_in, double y7_in, double z7_in,
                    double x8_in, double y8_in, double z8_in,
                    double dcx_in,double dcy_in,double dcz_in,
                    double minX_in,double minY_in,double minZ_in,
                    double maxX_in,double maxY_in,double maxZ_in
                  )
    : Cell( material_in
            , formulation_in
            , alpha_in
            , nGPoints_in
          )
    , minX(minX_in)
    , minY(minY_in)
    , minZ(minZ_in)
    , maxX(maxX_in)
    , maxY(maxY_in)
    , maxZ(maxZ_in)
{
    points.resize(8,3);
    points.writeElement(x1_in,0,0);
    points(0,1) = y1_in;
    points(0,2) = z1_in;
    points(1,0) = x2_in;
    points(1,1) = y2_in;
    points(1,2) = z2_in;
    points(2,0) = x3_in;
    points(2,1) = y3_in;
    points(2,2) = z3_in;
    points(3,0) = x4_in;
    points(3,1) = y4_in;
    points(3,2) = z4_in;
    points.writeElement(x5_in,4,0);
    points(4,1) = y5_in;
    points(4,2) = z5_in;
    points(5,0) = x6_in;
    points(5,1) = y6_in;
    points(5,2) = z6_in;
    points(6,0) = x7_in;
    points(6,1) = y7_in;
    points(6,2) = z7_in;
    points(7,0) = x8_in;
    points(7,1) = y8_in;
    points(7,2) = z8_in;

    // only valid for rectangular blocks
    Ax = x2_in - x1_in;
    Ay = y3_in - y2_in;
    Az = z5_in - z1_in;

    this->jacobian = Ax/2. * Ay/2. * Az/2.;

    this->createGaussPoints( dcx_in, dcy_in, dcz_in );
}

CellRect::~CellRect()
{
}


void CellRect::createGaussPoints( double dcx_in, double dcy_in, double dcz_in )
{
    lmx::DenseMatrix<double> gCoef( (size_type)nGPoints,2); // col 1: offset from center point,  col 2: weight factor.

    if (nGPoints == 1) {
        gCoef(0,0) = 0.;
        gCoef(0,1) = 1.;
    }
    else if(nGPoints == 2) {
        gCoef(0,0) = -std::sqrt(1./3.);
        gCoef(0,1) = 1.;
        gCoef(1,0) =  std::sqrt(1./3.);
        gCoef(1,1) = 1.;
    }
    else if(nGPoints == 3) {
        gCoef(0,0) = -0.77459666924148337703;
        gCoef(0,1) = 0.88888888888888888888;
        gCoef(1,0) =  0.00000000000000000000;
        gCoef(1,1) = 0.55555555555555555555;
        gCoef(2,0) =  0.77459666924148337703;
        gCoef(2,1) = 0.88888888888888888888;
    }
    else if(nGPoints == 4) {
        gCoef(0,0) = -0.86113631159405257522;
        gCoef(0,1) = 0.34785484513745385737;
        gCoef(1,0) = -0.33998104358485626480;
        gCoef(1,1) = 0.65214515486254614262;
        gCoef(2,0) =  0.33998104358485626480;
        gCoef(2,1) = 0.65214515486254614262;
        gCoef(3,0) =  0.86113631159405257522;
        gCoef(3,1) = 0.34785484513745385737;
    }
    else if(nGPoints == 5) {
        gCoef(0,0) = -0.90617984593866399279;
        gCoef(0,1) = 0.23692688505618908751;
        gCoef(1,0) = -0.53846931010568309103;
        gCoef(1,1) = 0.47862867049936646804;
        gCoef(2,0) =  0.00000000000000000000;
        gCoef(2,1) = 0.56888888888888888889;
        gCoef(3,0) =  0.53846931010568309103;
        gCoef(2,1) = 0.47862867049936646804;
        gCoef(4,0) =  0.90617984593866399279;
        gCoef(3,1) = 0.23692688505618908751;
    }
    else if(nGPoints == 6) {
        gCoef(0,0) = -0.932469514203152;
        gCoef(0,1) = 0.23692688505618908751;
        gCoef(1,0) = -0.661209386466265;
        gCoef(1,1) = 0.47862867049936646804;
        gCoef(2,0) = -0.238619186083197;
        gCoef(2,1) = 0.56888888888888888889;
        gCoef(3,0) = -gCoef(2,0)       ;
        gCoef(3,1) = gCoef(2,1);
        gCoef(4,0) = -gCoef(1,0)       ;
        gCoef(4,1) = gCoef(1,1);
        gCoef(5,0) = -gCoef(0,0)       ;
        gCoef(5,1) = gCoef(0,1);
    }
    else if(nGPoints == 8) {
        gCoef(0,0) = -0.96028985649753623168;
        gCoef(0,1) = 0.10122853629037625915;
        gCoef(1,0) = -0.79666647741362673959;
        gCoef(1,1) = 0.22238103445337447054;
        gCoef(2,0) = -0.52553240991632898581;
        gCoef(2,1) = 0.31370664587788728733;
        gCoef(3,0) = -0.18343464249564980493;
        gCoef(3,1) = 0.36268378337836198296;
        gCoef(4,0) = -gCoef(3,0)            ;
        gCoef(4,1) = gCoef(3,1);
        gCoef(5,0) = -gCoef(2,0)            ;
        gCoef(5,1) = gCoef(2,1);
        gCoef(6,0) = -gCoef(1,0)            ;
        gCoef(6,1) = gCoef(1,1);
        gCoef(7,0) = -gCoef(0,0)            ;
        gCoef(7,1) = gCoef(0,1);
    }
    else if(nGPoints == 12) {
        gCoef( 0,0) = -0.981560634246719;
        gCoef( 0,1) = 0.047175336386511;
        gCoef( 1,0) = -0.904117256370475;
        gCoef( 1,1) = 0.106939325995318;
        gCoef( 2,0) = -0.769902674194305;
        gCoef( 2,1) = 0.160078328543346;
        gCoef( 3,0) = -0.587317954286617;
        gCoef( 3,1) = 0.203167426723066;
        gCoef( 4,0) = -0.367831498998180;
        gCoef( 4,1) = 0.233492536538355;
        gCoef( 5,0) = -0.125233408511469;
        gCoef( 5,1) = 0.249147045813403;
        gCoef( 6,0) = -gCoef(5,0)       ;
        gCoef( 6,1) = gCoef(5,1);
        gCoef( 7,0) = -gCoef(4,0)       ;
        gCoef( 7,1) = gCoef(4,1);
        gCoef( 8,0) = -gCoef(3,0)       ;
        gCoef( 8,1) = gCoef(3,1);
        gCoef( 9,0) = -gCoef(2,0)       ;
        gCoef( 9,1) = gCoef(2,1);
        gCoef(10,0) = -gCoef(1,0)       ;
        gCoef(10,1) = gCoef(1,1);
        gCoef(11,0) = -gCoef(0,0)       ;
        gCoef(11,1) = gCoef(0,1);
    }
    else if(nGPoints == 20) {
        gCoef( 0,0) = -0.993128599185094924786;
        gCoef( 0,1) = 0.017614007139152118312;
        gCoef( 1,0) = -0.963971927277913791268;
        gCoef( 1,1) = 0.040601429800386941331;
        gCoef( 2,0) = -0.912234428251325905868;
        gCoef( 2,1) = 0.062672048334109063570;
        gCoef( 3,0) = -0.839116971822218823395;
        gCoef( 3,1) = 0.083276741576704748725;
        gCoef( 4,0) = -0.746331906460150792614;
        gCoef( 4,1) = 0.101930119817240435037;
        gCoef( 5,0) = -0.636053680726515025453;
        gCoef( 5,1) = 0.118194531961518417312;
        gCoef( 6,0) = -0.510867001950827098004;
        gCoef( 6,1) = 0.131688638449176626898;
        gCoef( 7,0) = -0.373706088715419560673;
        gCoef( 7,1) = 0.142096109318382051329;
        gCoef( 8,0) = -0.227785851141645078080;
        gCoef( 8,1) = 0.149172986472603746788;
        gCoef( 9,0) = -0.076526521133497333755;
        gCoef( 9,1) = 0.152753387130725850698;
        gCoef(10,0) = -gCoef(9,0)             ;
        gCoef(10,1) = gCoef(9,1);
        gCoef(11,0) = -gCoef(8,0)             ;
        gCoef(11,1) = gCoef(8,1);
        gCoef(12,0) = -gCoef(7,0)             ;
        gCoef(12,1) = gCoef(7,1);
        gCoef(13,0) = -gCoef(6,0)             ;
        gCoef(13,1) = gCoef(6,1);
        gCoef(14,0) = -gCoef(5,0)             ;
        gCoef(14,1) = gCoef(5,1);
        gCoef(15,0) = -gCoef(4,0)             ;
        gCoef(15,1) = gCoef(4,1);
        gCoef(16,0) = -gCoef(3,0)             ;
        gCoef(16,1) = gCoef(3,1);
        gCoef(17,0) = -gCoef(2,0)             ;
        gCoef(17,1) = gCoef(2,1);
        gCoef(18,0) = -gCoef(1,0)             ;
        gCoef(18,1) = gCoef(1,1);
        gCoef(19,0) = -gCoef(0,0)             ;
        gCoef(19,1) = gCoef(0,1);
    }

    for (unsigned int i=0; i<gCoef.rows(); ++i) {
        for (unsigned int j=0; j<gCoef.rows(); ++j) {
            gPoints.push_back( new GaussPoint2D( this->alpha,
                                                 gCoef(i,1)*gCoef(j,1),
                                                 jacobian,
                                                 mat,
                                                 i * nGPoints + j,
                                                 points.readElement(0,0) + Ax/2 * ( 1 + gCoef(i,0) ) ,
                                                 points.readElement(0,1) + Ay/2 * ( 1 + gCoef(j,0) ),
                                                 dcx_in,
						 true
                                               )
                             );

        }
    }
}


//void CellRect::initialize( std::map<int,Node*> & nodes_in )
//{
////   cout << "CellRect points " << this->points << endl;
//
//  // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
//  for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
//        it != gPoints.end();
//        ++it)
//  {
//    (*it)->findSupportNodes( nodes_in, minX, maxX, minY, maxY );
//  }
//}

void CellRect::initialize( std::vector<Node*> & nodes_in )
{
//   cout << "CellRect points " << this->points << endl;

    // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->findSupportNodes( nodes_in, minX, maxX, minY, maxY );
//    (*it)->findSupportNodes( nodes_in );
    }
}

void CellRect::gnuplotOut( std::ofstream& data, std::ofstream& gpdata )
{
//   cout << "Cell points" << points.rows() << ", " << points.cols() << endl;
//   cout << "Cell points " << this->points << endl;

    for (unsigned int i=0; i<points.rows(); ++i) {
        data << points(i,0) << " " << points(i,1) << " 0\n";
    }
    data << points(0,0) << " " << points(0,1) << " 0\n";
    data << std::endl;

    for(std::vector<GaussPoint*>::iterator it=gPoints.begin();
            it!=gPoints.end();
            ++it) {
        (*it)->gnuplotOut( gpdata );
    }
}


} //Namespace mknix
