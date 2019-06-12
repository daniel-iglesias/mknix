//-- Licencia --
#include "celltetrahedron.h"
#include "material.h"
#include "node.h"
#include "gausspoint3D.h"
#include <string>

namespace mknix
{

CellTetrahedron::CellTetrahedron()
{
}


CellTetrahedron::CellTetrahedron( Material& material_in,
                                  std::string formulation_in,
                                  double alpha_in,
                                  int nGPoints_in,
                                  Point* n1_in,
                                  Point* n2_in,
                                  Point* n3_in,
                                  Point* n4_in
                                )
    : Cell( material_in
            , formulation_in
            , alpha_in
            , nGPoints_in
          )
{
    bodyPoints.push_back( n1_in );
    bodyPoints.push_back( n2_in );
    bodyPoints.push_back( n3_in );
    bodyPoints.push_back( n4_in );

    points(0,0) = n1_in->getX() - n4_in->getX();
    points(0,1) = n1_in->getY() - n4_in->getY();
    points(0,2) = n1_in->getZ() - n4_in->getZ();

    points(1,0) = n2_in->getX() - n4_in->getX();
    points(1,1) = n2_in->getY() - n4_in->getY();
    points(1,2) = n2_in->getZ() - n4_in->getZ();

    points(2,0) = n3_in->getX() - n4_in->getX();
    points(2,1) = n3_in->getY() - n4_in->getY();
    points(2,2) = n3_in->getZ() - n4_in->getZ();

// From Cramer's rule (source http://en.wikipedia.org/wiki/Parallelepiped)
    this->jacobian = points.determinant() / 6.;

    dc = ( n1_in->distance( *n2_in )
           +n1_in->distance( *n3_in )
           +n1_in->distance( *n4_in )
           +n2_in->distance( *n3_in )
           +n2_in->distance( *n4_in )
           +n4_in->distance( *n4_in ) ) / 6.;
//   dc = 0.9;
//  cout << points(0,2) << ", "
//  << points(1,2) << ", "
//  << points(2,2) << "\n ";
//  cout << "dc = "<< dc << "\n";

    this->createGaussPoints( );
}

CellTetrahedron::~CellTetrahedron()
{
}


void CellTetrahedron::createGaussPoints( )
{
    lmx::DenseMatrix<double> gCoef( size_type(nGPoints), 5);

    // reference: http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
    if (nGPoints == 1)
    {
        gCoef(0,0) = .25;
        gCoef(0,1) = .25;
        gCoef(0,2) = .25;
        gCoef(0,3) = .25;
        gCoef(0,4) = 1.;
    }
    else if(nGPoints == 4)
    {
        gCoef(0,0) = .585410196624969;
        gCoef(0,1) = gCoef(0,2) = gCoef(0,3) = .138196601125011;
        gCoef(0,4) = .25;
        gCoef(1,1) = .585410196624969;
        gCoef(1,0) = gCoef(1,2) = gCoef(1,3) = .138196601125011;
        gCoef(1,4) = .25;
        gCoef(2,2) = .585410196624969;
        gCoef(2,1) = gCoef(2,0) = gCoef(2,3) = .138196601125011;
        gCoef(2,4) = .25;
        gCoef(3,3) = .585410196624969;
        gCoef(3,1) = gCoef(3,2) = gCoef(3,0) = .138196601125011;
        gCoef(3,4) = .25;
    }
// more on ... http://electromagnetics.biz/2D%20Gauss.txt
//   else if(nGPoints == 12){
//   }
//   else if(nGPoints == 20){
//   }

    for (int i=0; i<nGPoints; ++i)
    {
        gPoints.push_back
        ( new GaussPoint3D
          ( this->alpha,
            gCoef(i,4),
            jacobian,
            mat,
            i,
            bodyPoints[0]->getX() * gCoef(i,0)
            + bodyPoints[1]->getX() * gCoef(i,1)
            + bodyPoints[2]->getX() * gCoef(i,2)
            + bodyPoints[3]->getX() * gCoef(i,3),
            bodyPoints[0]->getY() * gCoef(i,0)
            + bodyPoints[1]->getY() * gCoef(i,1)
            + bodyPoints[2]->getY() * gCoef(i,2)
            + bodyPoints[3]->getY() * gCoef(i,3),
            bodyPoints[0]->getZ() * gCoef(i,0)
            + bodyPoints[1]->getZ() * gCoef(i,1)
            + bodyPoints[2]->getZ() * gCoef(i,2)
            + bodyPoints[3]->getZ() * gCoef(i,3),
            dc,
            true
          )
        );
    }
}


// void CellTetrahedron::initialize( std::map<int,Point*> & nodes_in )
// {
// //   cout << "CellTetrahedron points " << this->points << endl;
//
//   // This function can be joined with assembleGaussPoints so the
//   // Gpoints are iterated only once...
//   for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
//         it != gPoints.end();
//         ++it)
//   {
//     (*it)->findSupportPoints( nodes_in );
//   }
// }



void CellTetrahedron::gnuplotOut( std::ofstream& data, std::ofstream& gpdata )
{
//  for (int i=0; i<3; ++i){
//    data << points(i,0) << " " << points(i,1) << " 0\n";
//  }
//  data << points(0,0) << " " << points(0,1) << " 0\n";
//  data << std::endl;
//
//  for(std::vector<GaussPoint*>::iterator it=gPoints.begin();
//      it!=gPoints.end();
//      ++it){
//        (*it)->gnuplotOut( gpdata );
//      }
}

} //Namespace mknix
