//-- Licencia --
#include "gausspoint3D.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctiontetrahedron.h"

#include <simulation/simulation.h>
#include <system/system.h>
#include <system/loadthermalbody.h>

namespace mknix {

GaussPoint3D::GaussPoint3D()
{
}


GaussPoint3D::GaussPoint3D(double alpha_in,
                           double weight_in,
                           double jacobian_in,
                           Material * mat_in,
                           int num_in,
                           double coor_x,
                           double coor_y,
                           double coor_z,
                           double dc_in, bool stressPoint_in
)
        : GaussPoint(3,
                     alpha_in,
                     weight_in,
                     jacobian_in,
                     mat_in,
                     num_in,
                     coor_x,
                     coor_y,
                     coor_z,
                     dc_in,
                     stressPoint_in
)
{
//   this->jacobian = jacobian_in;
//  cout << "jacobian = " << jacobian << endl;
    F3.beUnityTensor();
}

GaussPoint3D::~GaussPoint3D()
{
}


void GaussPoint3D::shapeFunSolve(std::string type_in, double q_in)
{
    initializeMatVecs();
    GaussPoint::shapeFunSolve(type_in, q_in);
    //////////////// Filling matrix B = L * Phi
    for (auto i = 0u; i < supportNodesSize; ++i) {
        B.writeElement(shapeFun->getPhi(1, i), 0, 0 + 3 * i);
        B.writeElement(shapeFun->getPhi(2, i), 1, 1 + 3 * i);
        B.writeElement(shapeFun->getPhi(3, i), 2, 2 + 3 * i);
        B.writeElement(shapeFun->getPhi(2, i), 3, 0 + 3 * i);
        B.writeElement(shapeFun->getPhi(1, i), 3, 1 + 3 * i);
        B.writeElement(shapeFun->getPhi(3, i), 4, 1 + 3 * i);
        B.writeElement(shapeFun->getPhi(2, i), 4, 2 + 3 * i);
        B.writeElement(shapeFun->getPhi(3, i), 5, 0 + 3 * i);
        B.writeElement(shapeFun->getPhi(1, i), 5, 2 + 3 * i);
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


void GaussPoint3D::fillFEmatrices()
{
    initializeMatVecs();

    if (!shapeFun) {
        shapeFun = new ShapeFunctionTetrahedron(this);
    }
    shapeFun->calc();
    //////////////// Filling matrix B = L * Phi
    /////////////// Dim(B) = 6x12
    for (auto i = 0u; i < supportNodesSize; ++i) {
        B.writeElement(shapeFun->getPhi(1, i), 0, 0 + 3 * i);
        B.writeElement(shapeFun->getPhi(2, i), 1, 1 + 3 * i);
        B.writeElement(shapeFun->getPhi(3, i), 2, 2 + 3 * i);
        B.writeElement(shapeFun->getPhi(2, i), 3, 0 + 3 * i);
        B.writeElement(shapeFun->getPhi(1, i), 3, 1 + 3 * i);
        B.writeElement(shapeFun->getPhi(3, i), 4, 1 + 3 * i);
        B.writeElement(shapeFun->getPhi(2, i), 4, 2 + 3 * i);
        B.writeElement(shapeFun->getPhi(3, i), 5, 0 + 3 * i);
        B.writeElement(shapeFun->getPhi(1, i), 5, 2 + 3 * i);
        if (stressPoint) {
            if (Simulation::getSmoothingType() == "CONSTANT") {
                supportNodes[i]->addWeight(weight * jacobian);
            }
            if (Simulation::getSmoothingType() == "LOCAL") {
                supportNodes[i]->addWeight(weight * shapeFun->getPhi(0, i) * jacobian);
            }
        }
    }
}


void GaussPoint3D::computeMij()
{
    //////////////// Calculation of Mass matrix:
    //////////////// M = rho * wg * N^T * N * |Jc|
    //////////////// Mij = rho * wg * Ni * Nj * |Jc| = M(3*i + m, 3*j + n), m,n = 0,1.
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            M.writeElement(mat->getDensity() * weight * shapeFun->getPhi(0, i) * shapeFun->getPhi(0, j) * jacobian,
                           3 * i, 3 * j);
            M.writeElement(mat->getDensity() * weight * shapeFun->getPhi(0, i) * shapeFun->getPhi(0, j) * jacobian,
                           3 * i + 1, 3 * j + 1);
            M.writeElement(mat->getDensity() * weight * shapeFun->getPhi(0, i) * shapeFun->getPhi(0, j) * jacobian,
                           3 * i + 2, 3 * j + 2);
        }
    }
//   cout << "X "<< this->X << ", Y " << this->Y << endl;
//   cout << "Jacobian: " << jacobian << endl;
//   cout << "weight: " << weight << endl;
//   cout << "matrix M: " << M << endl;

}


void GaussPoint3D::computeKij()
{
    K.reset();
    // Kij = wg * Bi^T * D * Bj * |Jc| = K(3*i + m, 3*j + n)
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 0u; m < 3; ++m) {
                for (auto n = 0u; n < 3; ++n) {
                    for (auto a = 0u; a < 6; ++a) {
                        for (auto b = 0u; b < 6; ++b) {
                            // K(3*i + m, 3*j + n) = Kij(m,n)
                            // = w * J * Bi(a,m) * D(a,b) * Bj(b,n);
                            K.addElement(weight *
                                         B.readElement(a, m + 3 * i) *
                                         mat->getD().readElement(a, b) *
                                         B.readElement(b, n + 3 * j) *
                                         jacobian,
                                         3 * i + m,
                                         3 * j + n);
                        }
                    }
                }
            }
            //       cout<< endl << "Ks, node=" << this->i << ", i=" << i << ", j=" << j
            //           << endl << Kij << endl << KijTest << endl;
        }
    }
}


void GaussPoint3D::computeStress()
{
    tension.clear();

    for (auto i = 0u; i < supportNodesSize; ++i) {
        tension.addElement(
                mat->getD().readElement(0, 0) * B.readElement(0, 0 + 2 * i) * supportNodes[i]->getUx() +
                mat->getD().readElement(0, 1) * B.readElement(1, 1 + 2 * i) * supportNodes[i]->getUy() +
                mat->getD().readElement(0, 2) * B.readElement(2, 2 + 2 * i) * supportNodes[i]->getUz(),
                0);
        tension.addElement(
                mat->getD().readElement(1, 0) * B.readElement(0, 0 + 2 * i) * supportNodes[i]->getUx() +
                mat->getD().readElement(1, 1) * B.readElement(1, 1 + 2 * i) * supportNodes[i]->getUy() +
                mat->getD().readElement(1, 2) * B.readElement(2, 2 + 2 * i) * supportNodes[i]->getUz(),
                1);
        tension.addElement(
                mat->getD().readElement(2, 0) * B.readElement(0, 0 + 2 * i) * supportNodes[i]->getUx() +
                mat->getD().readElement(2, 1) * B.readElement(1, 1 + 2 * i) * supportNodes[i]->getUy() +
                mat->getD().readElement(2, 2) * B.readElement(2, 2 + 2 * i) * supportNodes[i]->getUz(),
                2);
        tension.addElement(
                mat->getD().readElement(3, 3) * B.readElement(3, 0 + 2 * i) * supportNodes[i]->getUx() +
                mat->getD().readElement(3, 3) * B.readElement(3, 1 + 2 * i) * supportNodes[i]->getUy(),
                3);
        tension.addElement(
                mat->getD().readElement(4, 4) * B.readElement(4, 1 + 2 * i) * supportNodes[i]->getUy() +
                mat->getD().readElement(4, 4) * B.readElement(4, 2 + 2 * i) * supportNodes[i]->getUz(),
                4);
        tension.addElement(
                mat->getD().readElement(5, 5) * B.readElement(5, 0 + 2 * i) * supportNodes[i]->getUx() +
                mat->getD().readElement(5, 5) * B.readElement(5, 2 + 2 * i) * supportNodes[i]->getUz(),
                5);
    }
//  cout << "\nVector tension: \n" << tension << endl;

    if (Simulation::getSmoothingType() == "CONSTANT") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            for (auto j = 0u; j < 6; ++j) {
                r.writeElement(weight / supportNodes[i]->getWeight()
                               * tension.readElement(j) * jacobian, 6 * i + j);
            }
        }
    }
}

void GaussPoint3D::computeNLStress()
{
    //   sigma2.beProductTraOf( P3, F3 );
    sigma3.beProductTraOf(F3 * S3, F3);
    sigma3 *= 1. / F3.determinant();

//  cout << "\nVector tension: \n" << sigma2 << endl;
    if (Simulation::getSmoothingType() == "CONSTANT") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma3(0, 0) * jacobian, 6 * i + 0);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma3(1, 1) * jacobian, 6 * i + 1);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma3(2, 2) * jacobian, 6 * i + 2);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma3(0, 1) * jacobian, 6 * i + 3);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma3(1, 2) * jacobian, 6 * i + 4);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * sigma3(2, 2) * jacobian, 6 * i + 5);
        }
    }
    else if (Simulation::getSmoothingType() == "LOCAL") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma3(0, 0) * jacobian, 6 * i + 0);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma3(1, 1) * jacobian, 6 * i + 1);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma3(2, 2) * jacobian, 6 * i + 2);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma3(0, 1) * jacobian, 6 * i + 3);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma3(1, 2) * jacobian, 6 * i + 4);
            r.writeElement(weight / supportNodes[i]->getWeight()
                           * shapeFun->getPhi(0, i) * sigma3(2, 2) * jacobian, 6 * i + 5);
        }
    }
    else if (Simulation::getSmoothingType() == "GLOBAL") {
        for (auto i = 0u; i < supportNodesSize; ++i) {
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma3(0, 0) * jacobian, 6 * i + 0);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma3(1, 1) * jacobian, 6 * i + 1);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma3(2, 2) * jacobian, 6 * i + 2);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma3(0, 1) * jacobian, 6 * i + 3);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma3(1, 2) * jacobian, 6 * i + 4);
            r.writeElement(weight * mat->getDensity()
                           * shapeFun->getPhi(0, i) * sigma3(2, 2) * jacobian, 6 * i + 5);
        }
    }
}

void GaussPoint3D::computeFint()
{
    fint.reset();
    computeKij();
    lmx::Vector<data_type> disp(3 * supportNodesSize);
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto j = 0u; j < 3; ++j) {
            disp.writeElement(supportNodes[a]->getU(j), 3 * a + j);
        }
    }
//  cout << "disp = " << disp;

    fint.mult(K, disp);
}


void GaussPoint3D::computeFext()
{
    // Mass matrix must be computed previously
    fext.reset();
    lmx::Vector<data_type> gravity(3 * supportNodesSize);
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto j = 0u; j < 3; ++j) {
            gravity.writeElement(-Simulation::getGravity(j), 3 * a + j);
        }
    }
    fext.mult(M, gravity);

//   cout << "M:" << M << endl;
//   cout << "gravity:" << gravity << endl;
}


void GaussPoint3D::computeNLFint()
{
    F3.zero();
    for (auto i = 0u; i < 3; ++i) {
        for (auto j = 0u; j < 3; ++j) {
            for (auto k = 0u; k < supportNodesSize; ++k) {
                F3(i, j) += supportNodes[k]->getqx(i) * shapeFun->getPhi(j + 1, k);
            }
        }
    }
    mat->computeS(S3, F3);
    P3.beProductOf(F3, S3);
    fint.reset();
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto j = 0u; j < 3; ++j) {
            for (auto k = 0u; k < 3; ++k) { //dot product
                fint.addElement(weight * jacobian * P3(j, k) * shapeFun->getPhi(k + 1, a),
                                3 * a + j
                ); //maybe P(j,i)?
            }
        }
    }
}

void GaussPoint3D::computeNLKij()
{
    // It fills the shared K matrix...
    // Constitutive component (similar to linear case):
    data_type temp1, temp2;
    K.reset();


    ////////////// Calculation of Tangent matrix:
    ////////////// Kab = wg * |Jc| * phidot_a * csym * phidot_b
    //////////////     = K(3*a + m, 3*b + n)
    //Exactly the same as in 2D but less efficient...
    for (auto a = 0u; a < supportNodesSize; ++a) {
        for (auto b = 0u; b < supportNodesSize; ++b) {
            for (auto m = 0u; m < 3; ++m) {
                for (auto n = 0u; n < 3; ++n) {
                    for (auto i = 0u; i < 3; ++i) {
                        for (auto j = 0u; j < 3; ++j) {
//                             for (  k=0; k<3; ++k ) {
//                                 for ( l=0; l<3; ++l ) {
// //                   K(3*a + m, 3*b + n) = Kab(m,n);
//                                     K.addElement( weight * jacobian * ( F3(m,i)
//                                                        * shapeFun->getPhi(j+1,a) * mat->getCsym(i,j,k,l)
//                                                        * shapeFun->getPhi(k+1,b) * F3(n,l) ),
// 						  3*a+m, 3*b+n);
//                                 }
//                             }
                            //optimized
                            // temp2 (index 3,4,5) was sometimes zero.
                            // mat->getD() was called instead of mat->getC()
                            // (Bug probably fixed)
                            temp1 = B.readElement(i, i + 3 * a);
                            temp1 *= mat->getC().readElement(i, j);
                            temp1 *= B.readElement(j, j + 3 * b);
                            temp2 = B.readElement(3, i + 3 * a);
                            temp2 *= mat->getC().readElement(3, 3);
                            temp2 *= B.readElement(3, j + 3 * b);
                            temp1 += temp2;
                            temp2 = B.readElement(4, i + 3 * a);
                            temp2 *= mat->getC().readElement(4, 4);
                            temp2 *= B.readElement(4, j + 3 * b);
                            temp1 += temp2;
                            temp2 = B.readElement(5, i + 3 * a);
                            temp2 *= mat->getC().readElement(5, 5);
                            temp2 *= B.readElement(5, j + 3 * b);
                            temp1 += temp2;
                            temp1 *= weight;
                            temp1 *= F3(m, i);
                            temp1 *= F3(n, j);
                            temp1 *= jacobian;
                            K.addElement(temp1, 3 * a + m, 3 * b + n);
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
            for (auto i = 0u; i < 3; ++i) {
                for (auto j = 0u; j < 3; ++j) {
//          phiSphi += shapeFun->getPhi(i+1,a) * S2(i,j) *
//                     shapeFun->getPhi(j+1,b);
                    phiSphi += shapeFun->getPhi(i + 1, a) * S3(i, j)
                               * shapeFun->getPhi(j + 1, b);
                }
            }
            phiSphi *= weight;
            phiSphi *= jacobian;
            K.addElement(phiSphi, 3 * a + 0, 3 * b + 0);
            K.addElement(phiSphi, 3 * a + 1, 3 * b + 1);
            K.addElement(phiSphi, 3 * a + 2, 3 * b + 2);
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


void GaussPoint3D::assembleMij(lmx::Matrix<data_type>& globalMass)
{
// 	cout << "MATRIX M in assembly:" << M << endl;
    for (auto i = 0u; i < supportNodesSize; ++i) {
//     cout<<"node: "<< supportNodes[i]->getNumber() << " "
//         <<supportNodes[i]->getX() << " "
//         <<supportNodes[i]->getY() << " " << endl ;
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 0u; m < 3; ++m) {
                for (auto n = 0u; n < 3; ++n) {
                    globalMass.addElement(M.readElement(3 * i + m, 3 * j + n),
                                          3 * supportNodes[i]->getNumber() + m,
                                          3 * supportNodes[j]->getNumber() + n
                    );
                }
            }
        }
    }
//   cout << "globalMass " << globalMass << endl;
}


void GaussPoint3D::assembleKij(lmx::Matrix<data_type>& globalTangent)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
//     cout<<"node: "<< supportNodes[i]->nodeNumber() << " "
//         <<supportNodes[i]->getx() << " "
//         <<supportNodes[i]->gety() << " " << endl ;
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 0u; m < 3; ++m) {
                for (auto n = 0u; n < 3; ++n) {
                    globalTangent.addElement(K.readElement(3 * i + m, 3 * j + n),
                                             3 * supportNodes[i]->getNumber() + m,
                                             3 * supportNodes[j]->getNumber() + n
                    );
                }
            }
        }
    }
//   cout << globalTangent << endl;
}

void GaussPoint3D::assembleRi(lmx::Vector<data_type>& bodyR, int firstNode)
{
//   cout << "Size = " << bodyR.size() << "solving tensions..." << endl;
//   cout << "Size = " << r.size() << "solving tensions..." << endl;

    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 6; ++m) {
            bodyR.addElement(r.readElement(6 * i + m),
                             6 * (supportNodes[i]->getNumber() - firstNode) + m
            );
        }
    }
}


void GaussPoint3D::assembleFint(lmx::Vector<data_type>& globalFint)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 3; ++m) {
            globalFint.addElement(fint.readElement(3 * i + m),
                                  3 * supportNodes[i]->getNumber() + m
            );
        }
    }
}


void GaussPoint3D::assembleFext(lmx::Vector<data_type>& globalFext)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 3; ++m) {
            globalFext.addElement(fext.readElement(3 * i + m),
                                  3 * supportNodes[i]->getNumber() + m
            );
        }
    }
}


double GaussPoint3D::calcPotentialE(const lmx::Vector<data_type>& q)
{
    double potential = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 3; ++m) {
            potential += q.readElement(3 * supportNodes[i]->getNumber() + m)
                         * fext.readElement(3 * i + m);
        }
    }
    return potential;
}

double GaussPoint3D::calcPotentialE(const VectorX<data_type>& q)
{
    double potential = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto m = 0u; m < 3; ++m) {
            potential += q[3 * supportNodes[i]->getNumber() + m)]
                         * fext.readElement(3 * i + m);
        }
    }
    return potential;
}


double GaussPoint3D::calcKineticE(const lmx::Vector<data_type>& qdot)
{
    double kinetic = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto n = 0u; n < 3; ++n) {
            for (auto j = 0u; j < supportNodesSize; ++j) {
                for (auto m = 0u; m < 3; ++m) {
                    kinetic -= 0.5
                               * qdot.readElement(3 * supportNodes[i]->getNumber() + n)
                               * M.readElement(3 * i + n, 3 * j + m)
                               * qdot.readElement(3 * supportNodes[j]->getNumber() + m);
                }
            }
        }
    }
    return kinetic;
}

double GaussPoint3D::calcElasticE()
{
    double elastic;
    elastic = mat->computeEnergy(F3);

    elastic *= -1. * weight * jacobian;

//  cout << "GP " << i << ": E_elastic = " << elastic << endl
//      << "F = " << F;

    return elastic;

}


void GaussPoint3D::initializeMatVecs() // Must be optimized for mech and thermal independency
{
    B.resize(6, 3 * supportNodesSize);
    tension.resize(6);
    r.resize(6 * supportNodesSize);

    C.resize(supportNodesSize, supportNodesSize);
    H.resize(supportNodesSize, supportNodesSize);
    Qext.resize(supportNodesSize);
    M.resize(3 * supportNodesSize, 3 * supportNodesSize);
    K.resize(3 * supportNodesSize, 3 * supportNodesSize);
    fext.resize(3 * supportNodesSize);
    fint.resize(3 * supportNodesSize);
}


} //Namespace mknix
