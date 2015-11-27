//-- Licencia --
#include "gausspoint2D.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctiontriangle.h"
#include "simulation.h"
#include "system.h"
#include "loadthermalbody.h"

namespace mknix {

GaussPoint2D::GaussPoint2D()
{
}


GaussPoint2D::GaussPoint2D(double alpha_in, double weight_in, double jacobian_in, Material * mat_in, int num_in,
                           double coor_x, double coor_y, double dc_in, bool stressPoint_in
)
        : GaussPoint(2, alpha_in, weight_in, jacobian_in, mat_in, num_in, coor_x, coor_y, dc_in, stressPoint_in
)
{
//   this->jacobian = jacobian_in;
//  cout << "jacobian = " << jacobian << endl;
    F2.beUnityTensor();
}

GaussPoint2D::GaussPoint2D(double alpha_in, double weight_in, double jacobian_in, Material * mat_in, int num_in,
                           double coor_x, double coor_y, double coor_z, double dc_in, bool stressPoint_in
)
        : GaussPoint(2, alpha_in, weight_in, jacobian_in, mat_in, num_in, coor_x, coor_y, coor_z, dc_in, stressPoint_in
)
{
//   this->jacobian = jacobian_in;
//  cout << "jacobian = " << jacobian << endl;
    F2.beUnityTensor();
}

GaussPoint2D::~GaussPoint2D()
{
}


void GaussPoint2D::shapeFunSolve(std::string type_in, double q_in)
{
    initializeMatVecs();
    GaussPoint::shapeFunSolve(type_in, q_in);
    //////////////// Filling matrix B = L * Phi
    for (auto i = 0u; i < supportNodesSize; ++i) {
        B.writeElement(shapeFun->getPhi(1, i), 0, 0 + 2 * i);
        B.writeElement(shapeFun->getPhi(2, i), 1, 1 + 2 * i);
        B.writeElement(shapeFun->getPhi(2, i), 2, 0 + 2 * i);
        B.writeElement(shapeFun->getPhi(1, i), 2, 1 + 2 * i);
//	  cout << supportNodes[i]->getx() << "\t"
//			  << supportNodes[i]->gety() << "\t"
//			  << shapeFun->getPhi(0,i) << "\t"
//			  << shapeFun->getPhi(1,i) << "\t"
//			  << shapeFun->getPhi(2,i) << "\n";
//      cout << B << endl;
    }
    if (stressPoint) {
        if (Simulation::getSmoothingType() == "CONSTANT") {
            for (auto i = 0u; i < supportNodesSize; ++i) {
                supportNodes[i]->addWeight(weight * jacobian);
            }
        }
        if (Simulation::getSmoothingType() == "LOCAL") {
            for (auto i = 0u; i < supportNodesSize; ++i) {
                supportNodes[i]->addWeight(weight * shapeFun->getPhi(0, i) * jacobian);
            }
        }
    }
}


void GaussPoint2D::fillFEmatrices()
{
    initializeMatVecs();
    cout << "GP in (" << this->X << ", " << this->Y << ")" << endl;

    if (!shapeFun) {
        shapeFun = new ShapeFunctionTriangle(this);
    }
    shapeFun->calc();
    //////////////// Filling matrix B = L * Phi
    /////////////// Dim(B) = 3x6
    for (auto i = 0u; i < supportNodesSize; ++i) {
// 	  cout << supportNodes[i]->getX() << "\t"
// 			  << supportNodes[i]->getY() << "\t"
// 			  << shapeFun->getPhi(0,i) << "\t"
// 			  << shapeFun->getPhi(1,i) << "\t"
// 			  << shapeFun->getPhi(2,i) << "\n";
        B.writeElement(shapeFun->getPhi(1, i), 0, 0 + 2 * i);
        B.writeElement(shapeFun->getPhi(2, i), 1, 1 + 2 * i);
        B.writeElement(shapeFun->getPhi(2, i), 2, 0 + 2 * i);
        B.writeElement(shapeFun->getPhi(1, i), 2, 1 + 2 * i);
//    cout << B << endl;
        if (stressPoint) {
            if (Simulation::getSmoothingType() == "CONSTANT") {
                supportNodes[i]->addWeight(weight * jacobian);
// 	      cout << "i=" << i << ", weight = " << weight << ", jacobian = "<< jacobian << endl;
// 	      cout << "B matrix = " << B << endl;
            }
            if (Simulation::getSmoothingType() == "LOCAL") {
                supportNodes[i]->addWeight(weight * shapeFun->getPhi(0, i) * jacobian);
            }
        }
    }
}


void GaussPoint2D::computeMij()
{
    //////////////// Calculation of Mass matrix:
    //////////////// M = rho * wg * N^T * N * |Jc|
    //////////////// Mij = rho * wg * Ni * Nj * |Jc| = M(2*i + m, 2*j + n), m,n = 0,1.
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            M.writeElement(mat->getDensity() * weight * shapeFun->getPhi(0, i) * shapeFun->getPhi(0, j) * jacobian,
                           2 * i, 2 * j);
            M.writeElement(mat->getDensity() * weight * shapeFun->getPhi(0, i) * shapeFun->getPhi(0, j) * jacobian,
                           2 * i + 1, 2 * j + 1);
        }
    }
//   cout << "Jacobian: " << jacobian << endl;
//   cout << "weight: " << weight << endl;
//   cout << "matrix M: " << M << endl;

}


void GaussPoint2D::computeKij()
{
//////////////// Calculation of Tangent matrix:
// Kij = wg * Bi^T * D * Bj * |Jc| = K(2*i + m, 2*j + n)
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 0u; m < 2; ++m) {
                for (auto n = 0u; n < 2; ++n) {
                    //           K(2*i + m, 2*j + n) = Kij(m,n);
                    K.writeElement(
                            (B.readElement(m, m + 2 * i) * mat->getD().readElement(m, n) * B.readElement(n, n + 2 * j)
                             + B.readElement(2, m + 2 * i) * mat->getD().readElement(2, 2) * B.readElement(2, n + 2 * j)
                            ) * weight * jacobian,
                            2 * i + m,
                            2 * j + n);
                }
            }
            //       cout<< endl << "Ks, node=" << this->i << ", i=" << i << ", j=" << j
            //           << endl << Kij << endl << KijTest << endl;
        }
    }
//   cout << "Jacobian: " << jacobian << endl;
//   cout << "weight: " << weight << endl;
//   cout << "matrix K: " << K << endl;
}


void GaussPoint2D::computeStress()
{
    tension.clear();

    for (auto i = 0u; i < supportNodesSize; ++i) {
        tension.addElement(
                mat->getD().readElement(0, 0) * B.readElement(0, 0 + 2 * i) * supportNodes[i]->getUx()
                + mat->getD().readElement(0, 1) * B.readElement(1, 1 + 2 * i) * supportNodes[i]->getUy(),
                0);
        tension.addElement(
                mat->getD().readElement(1, 0) * B.readElement(0, 0 + 2 * i) * supportNodes[i]->getUx()
                + mat->getD().readElement(1, 1) * B.readElement(1, 1 + 2 * i) * supportNodes[i]->getUy(),
                1);
        tension.addElement(
                mat->getD().readElement(2, 2) * B.readElement(2, 0 + 2 * i) * supportNodes[i]->getUx()
                + mat->getD().readElement(2, 2) * B.readElement(2, 1 + 2 * i) * supportNodes[i]->getUy(),
                2);
// 	cout << "computeStress() " << i << " " << B.readElement(1,1+2*i) << " " << supportNodes[i]->getUy() << endl ;
    }
//     cout << " B = " << B << endl;
//     cout << "TENSION = " << tension << endl;

    if (Simulation::getSmoothingType() == "CONSTANT") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
// 	    cout << "i=" << i << ", weight = " << weight << ", jacobian = "<< jacobian 
// 		 << "supportNodes[i]->getWeight() = " << supportNodes[i]->getWeight() << endl;
            r.writeElement((weight * jacobian) / supportNodes[i]->getWeight()
                           * tension.readElement(0), 3 * i + 0);
            r.writeElement((weight * jacobian) / supportNodes[i]->getWeight()
                           * tension.readElement(1), 3 * i + 1);
            r.writeElement((weight * jacobian) / supportNodes[i]->getWeight()
                           * tension.readElement(2), 3 * i + 2);
        }
    }
    else if (Simulation::getSmoothingType() == "LOCAL") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * tension.readElement(0) * jacobian, 3 * i + 0);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * tension.readElement(1) * jacobian, 3 * i + 1);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * tension.readElement(2) * jacobian, 3 * i + 2);
        }
    }
    else if (Simulation::getSmoothingType() == "GLOBAL") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * tension.readElement(0) * jacobian, 3 * i + 0);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * tension.readElement(1) * jacobian, 3 * i + 1);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * tension.readElement(2) * jacobian, 3 * i + 2);
        }
    }
}

void GaussPoint2D::computeNLStress()
{
    sigma2.beProductTraOf(P2, F2);
//     sigma2.beProductTraOf( F2*S2, F2 );
    sigma2 *= 1. / F2.determinant();

//  cout << "\nVector tension: \n" << sigma2 << endl;
    if (Simulation::getSmoothingType() == "CONSTANT") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma2(0, 0) * jacobian, 3 * i + 0);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma2(1, 1) * jacobian, 3 * i + 1);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma2(0, 1) * jacobian, 3 * i + 2);
        }
    }
    else if (Simulation::getSmoothingType() == "LOCAL") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma2(0, 0) * jacobian, 3 * i + 0);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma2(1, 1) * jacobian, 3 * i + 1);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma2(0, 1) * jacobian, 3 * i + 2);
        }
    }
    else if (Simulation::getSmoothingType() == "GLOBAL") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma2(0, 0) * jacobian, 3 * i + 0);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma2(1, 1) * jacobian, 3 * i + 1);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma2(0, 1) * jacobian, 3 * i + 2);
        }
    }
}

void GaussPoint2D::computeFint()
{
    fint.fillIdentity(0.);
    computeKij();
    lmx::Vector<data_type> disp(2 * supportNodesSize);
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto j = 0u; j < 2; ++j) {
            disp.writeElement(supportNodes[a]->getU(j), 2 * a + j);
        }
    }
//  cout << "disp = " << disp;

    fint.mult(K, disp);
}


void GaussPoint2D::computeFext()
{
    // Mass matrix must be computed previously
    fext.fillIdentity(0.);
    lmx::Vector<data_type> gravity(2 * supportNodesSize);
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto j = 0u; j < 2; ++j) {
            gravity.writeElement(-Simulation::getGravity(j), 2 * a + j);
        }
    }
    fext.mult(M, gravity);

//   cout << "M:" << M << endl;
//   cout << "gravity:" << gravity << endl;
}


void GaussPoint2D::computeNLFint()
{
    F2.zero();
    for (auto i = 0u; i < 2; ++i) {
        for (auto j = 0u; j < 2; ++j) {
            for (auto k = 0u; k < supportNodesSize; ++k) {
                F2(i, j) += supportNodes[k]->getqx(i) * shapeFun->getPhi(j + 1, k);
            }
        }
    }
    mat->computeS(S2, F2, this->getTemp());
    P2.beProductOf(F2, S2);
    fint.fillIdentity(0.);
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto j = 0u; j < 2; ++j) {
            for (auto k = 0u; k < 2; ++k) { //dot product
                fint.addElement(weight * jacobian * P2(j, k) * shapeFun->getPhi(k + 1, a),
                                2 * a + j
                ); //maybe P(j,i)?
            }
        }
    }
}

void GaussPoint2D::computeNLKij()
{
    // It fills the shared K matrix...
    // Constitutive component (similar to linear case):
    data_type temp1, temp2;
    K.reset();

    ////////////// Calculation of Tangent matrix:
    ////////////// Kij = wg * |Jc| * phidot_i * csym * phidot_j
    //////////////     = K(2*i + m, 2*j + n)
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto b = 0u; b < supportNodesSize; ++b) {
            for (auto m = 0u; m < 2; ++m) {
                for (auto n = 0u; n < 2; ++n) {
                    for (auto i = 0u; i < 2; ++i) {
                        for (auto j = 0u; j < 2; ++j) {
//                K(2*a + m, 2*b + n) = Kab(m,n);
//                K(2*a + m, 2*b + n) +=
//                    weight * F(m,i)
//                    * (  B.readElement(i,i+2*a)
//                         * mat->getC().readElement(i,j)
//                         * B.readElement(j,j+2*b)
//                         +  B.readElement(2,i+2*a)
//                         * mat->getC().readElement(2,2)
//                         * B.readElement(2,j+2*b) )
//                    * F(n,j)
//                    * jacobian;
                            //optimized:
                            temp1 = B.readElement(i, i + 2 * a);
                            temp1 *= mat->getC().readElement(i, j);
                            temp1 *= B.readElement(j, j + 2 * b);
                            temp2 = B.readElement(2, i + 2 * a);
                            temp2 *= mat->getC().readElement(2, 2);
                            temp2 *= B.readElement(2, j + 2 * b);
                            temp1 += temp2;
                            temp1 *= weight;
                            temp1 *= F2(m, i);
                            temp1 *= F2(n, j);
                            temp1 *= jacobian;
                            K.addElement(temp1, 2 * a + m, 2 * b + n);
                        }
                    }
                }
            }
        }
    }
    // Initial Stress component:
    data_type phiSphi = 0.;
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto b = 0u; b < supportNodesSize; ++b) {
            for (auto i = 0u; i < 2; ++i) {
                for (auto j = 0u; j < 2; ++j) {
//          phiSphi += shapeFun->getPhi(i+1,a) * S2(i,j) *
//                     shapeFun->getPhi(j+1,b);
                    phiSphi += shapeFun->getPhi(i + 1, a) * S2(i, j)
                               * shapeFun->getPhi(j + 1, b);
                }
            }
            phiSphi *= weight;
            phiSphi *= jacobian;
            K.addElement(phiSphi, 2 * a + 0, 2 * b + 0);
            K.addElement(phiSphi, 2 * a + 1, 2 * b + 1);
            phiSphi = 0.;
        }
    }

    // NUMERICAL APROXIMATION...
//   double Ax = 1E-7;
//   cofe::TensorRank2<2,double> F_back, F_forw;
//   cofe::TensorRank2Sym<2,double> S_back, S_forw;
//   cofe::TensorRank2<2,double> P_back, P_forw;
//   lmx::Vector<double> fint_back, fint_forw;
//   fint_back.resize( 2*supportNodesSize );
//   fint_forw.resize( 2*supportNodesSize );
//
//   // The hard way...
//   for(int z=0; z<supportNodesSize; ++z){
//     for(int y=0; y<2; ++y){
//       if(y==0) supportNodes[z]->setUx( supportNodes[z]->getUx() + Ax );
//       else     supportNodes[z]->setUy( supportNodes[z]->getUy() + Ax );
//       F_back.zero();
//       F_forw.zero();
//       for(int i=0; i<2; ++i){
//         for(int j=0; j<2; ++j){
//           for(int k=0; k<supportNodesSize; ++k){
//               F_forw(i,j) += supportNodes[k]->getx(i) * shapeFun->getPhi(j+1,k);
//           }
//         }
//       }
//       if(y==0) supportNodes[z]->setUx( supportNodes[z]->getUx() - 2*Ax );
//       else     supportNodes[z]->setUy( supportNodes[z]->getUy() - 2*Ax );
//       for(int i=0; i<2; ++i){
//         for(int j=0; j<2; ++j){
//           for(int k=0; k<supportNodesSize; ++k){
//             F_back(i,j) += supportNodes[k]->getx(i) * shapeFun->getPhi(j+1,k);
//           }
//         }
//       }
//       if(y==0) supportNodes[z]->setUx( supportNodes[z]->getUx() + Ax );
//       else     supportNodes[z]->setUy( supportNodes[z]->getUy() + Ax );
//
//       mat->computeS( S_back, F_back );
//       mat->computeS( S_forw, F_forw );
//       P_back.beProductOf( F_back, S_back );
//       P_forw.beProductOf( F_forw, S_forw );
//
//       fint_back.fillIdentity( 0. );
//       fint_forw.fillIdentity( 0. );
//       for(int a=0; a<supportNodesSize; ++a){
//         for(int j=0; j<2; ++j){
//           for(int k=0; k<2; ++k){ //dot product
//             fint_back(2*a+j) += weight * jacobian * P_back(j,k) * shapeFun->getPhi(k+1,a); //maybe P(j,i)?
//             fint_forw(2*a+j) += weight * jacobian * P_forw(j,k) * shapeFun->getPhi(k+1,a); //maybe P(j,i)?
//           }
//         }
//       }
//       for ( int a=0; a<2*supportNodesSize; ++a ){
//         K(a,2*z+y) = ( fint_forw(a) - fint_back(a) ) / (2*Ax);
//         K(a,2*z+y) = ( fint_forw(a) - fint(a) ) / (1*Ax);
//       }
//     }
//   }

    // Maybe faster (if it works :-)...
//   for(int z=0; z<supportNodesSize; ++z){
//     for(int y=0; y<2; ++y){
//
//       F_back.zero();
//       F_forw.zero();
//       for(int j=0; j<2; ++j){
//         F_back(y,j) += (supportNodes[z]->getU(y) - Ax) * shapeFun->getPhi(j+1,z);
//         F_forw(y,j) += (supportNodes[z]->getU(y) + Ax) * shapeFun->getPhi(j+1,z);
//       }
//
//       mat->computeS( S_back, F_back );
//       mat->computeS( S_forw, F_forw );
//       P_back.beProductOf( F_back, S_back );
//       P_forw.beProductOf( F_forw, S_forw );
//
//       fint_back = fint;
//       fint_forw = fint;
//       for(int a=0; a<supportNodesSize; ++a){
//         for(int j=0; j<2; ++j){
//           for(int k=0; k<2; ++k){ //dot product
//             fint_back(2*a+j) += weight * jacobian * P_back(j,k) * shapeFun->getPhi(k+1,a);
//             fint_forw(2*a+j) += weight * jacobian * P_forw(j,k) * shapeFun->getPhi(k+1,a);
//           }
//         }
//       }
//       for ( int a=0; a<2*supportNodesSize; ++a ){
//         K(a,2*z+y) = ( fint_forw(a) - fint_back(a) ) / (2*Ax);
// //         K(a,2*z+y) = ( fint_forw(a) - fint(a) ) / (1*Ax);
//       }
//     }
//   }

}


void GaussPoint2D::assembleMij(lmx::Matrix<data_type>& globalMass)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
//     cout<<"node: "<< supportNodes[i]->nodeNumber() << " "
//         <<supportNodes[i]->getx() << " "
//         <<supportNodes[i]->gety() << " " << endl ;
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 0u; m < 2; ++m) {
                for (auto n = 0u; n < 2; ++n) {
                    globalMass.addElement(M.readElement(2 * i + m, 2 * j + n),
                                          2 * supportNodes[i]->getNumber() + m,
                                          2 * supportNodes[j]->getNumber() + n
                    );
                }
            }
        }
    }
//   cout << endl;
}


void GaussPoint2D::assembleKij(lmx::Matrix<data_type>& globalTangent)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
//     cout<<"node: "<< supportNodes[i]->nodeNumber() << " "
//         <<supportNodes[i]->getx() << " "
//         <<supportNodes[i]->gety() << " " << endl ;
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 0u; m < 2; ++m) {
                for (auto n = 0u; n < 2; ++n) {
                    globalTangent.addElement(K.readElement(2 * i + m, 2 * j + n),
                                             2 * supportNodes[i]->getNumber() + m,
                                             2 * supportNodes[j]->getNumber() + n
                    );
                }
            }
        }
    }
//   cout << globalTangent << endl;
}

void GaussPoint2D::assembleRi(lmx::Vector<data_type>& bodyR, int firstNode)
{
//   cout << "Size = " << bodyR.size() << "solving tensions..." << endl;
//   cout << "Size = " << r.size() << "solving tensions..." << endl;

//     cout << "bodyR before: " << bodyR << endl;
    if (Simulation::getSmoothingType() == "OFF") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            for (auto m = 0u; m < 3; ++m) {
                bodyR.writeElement(r.readElement(3 * i + m),
                                   3 * (supportNodes[i]->getNumber() - firstNode) + m
                );
            }
        }
    }
    else {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            for (auto m = 0u; m < 3; ++m) {
                bodyR.addElement(r.readElement(3 * i + m),
                                 3 * (supportNodes[i]->getNumber() - firstNode) + m
                );
            }
        }
    }
//     cout << "bodyR after: " << bodyR << endl;
}


void GaussPoint2D::assembleFint(lmx::Vector<data_type>& globalFint)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 2; ++m) {
            globalFint.addElement(fint.readElement(2 * i + m),
                                  2 * supportNodes[i]->getNumber() + m
            );
        }
    }
}


void GaussPoint2D::assembleFext(lmx::Vector<data_type>& globalFext)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 2; ++m) {
            globalFext.addElement(fext.readElement(2 * i + m),
                                  2 * supportNodes[i]->getNumber() + m
            );
        }
    }
}


double GaussPoint2D::calcPotentialE(const lmx::Vector<data_type>& q)
{
    double potential = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 2; ++m) {
            potential += q.readElement(2 * supportNodes[i]->getNumber() + m)
                         * fext.readElement(2 * i + m);
        }
    }
    return potential;
}


double GaussPoint2D::calcKineticE(const lmx::Vector<data_type>& qdot)
{
    double kinetic = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto n = 0u; n < 2; ++n) {
            for (auto j = 0u; j < supportNodesSize; ++j) {
                for (auto m = 0u; m < 2; ++m) {
                    kinetic -= 0.5
                               * qdot.readElement(2 * supportNodes[i]->getNumber() + n)
                               * M.readElement(2 * i + n, 2 * j + m)
                               * qdot.readElement(2 * supportNodes[j]->getNumber() + m);
                }
            }
        }
    }

    return kinetic;
}

double GaussPoint2D::calcElasticE()
{
    double elastic;
    elastic = mat->computeEnergy(F2);

    elastic *= -1. * weight * jacobian;

//  cout << "GP " << i << ": E_elastic = " << elastic << endl
//      << "F = " << F;

    return elastic;

}


void GaussPoint2D::initializeMatVecs() // TODO: Must be optimized for mech and thermal independency
{
    B.resize(3, 2 * supportNodesSize);
    tension.resize(3);
    r.resize(3 * supportNodesSize);

    C.resize(supportNodesSize, supportNodesSize);
    H.resize(supportNodesSize, supportNodesSize);
    Qext.resize(supportNodesSize);
    M.resize(2 * supportNodesSize, 2 * supportNodesSize);
    K.resize(2 * supportNodesSize, 2 * supportNodesSize);
    fext.resize(2 * supportNodesSize);
    fint.resize(2 * supportNodesSize);
}


} //Namespace mknix
