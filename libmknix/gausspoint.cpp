//-- Licencia --
#include "gausspoint.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctionMLS.h"
#include "shapefunctiontriangle.h"
#include "shapefunctiontetrahedron.h"
#include "simulation.h"
#include "system.h"
#include "loadthermalbody.h"

namespace mknix {

GaussPoint::GaussPoint()
{
}


GaussPoint::GaussPoint( int dim_in,
                        double alpha_in,
                        double weight_in,
                        double jacobian_in,
                        Material* mat_in,
                        int num_in,
                        double coor_x,
                        double coor_y,
                        double dc_in,
			bool stressPoint_in
                      )
    : Point(dim_in, num_in, coor_x, coor_y, 0., jacobian_in, alpha_in, dc_in)
    , weight(weight_in)
    , mat(mat_in)
    , stressPoint(stressPoint_in)
{
}

GaussPoint::GaussPoint( int dim_in,
                        double alpha_in,
                        double weight_in,
                        double jacobian_in,
                        Material* mat_in,
                        int num_in,
                        double coor_x,
                        double coor_y,
                        double coor_z,
                        double dc_in,
			bool stressPoint_in
                      )
    : Point(dim_in, num_in, coor_x, coor_y, coor_z, jacobian_in, alpha_in, dc_in)
    , weight(weight_in)
    , mat(mat_in)
    , stressPoint(stressPoint_in)
{
}

GaussPoint::~GaussPoint()
{
}


void GaussPoint::shapeFunSolve( std::string type_in, double q_in )
{
    q_in = 0.5; // TODO: Check WTF I did this!!
    if (!shapeFun) {
        if(type_in == "RBF"){
	    alphai=3.5; // For the validation triangle, this works better: phi closer to 1.
            this->shapeFun = new ShapeFunctionRBF(supportNodesSize,
                                                  0,
                                                  0, // RBF type
                                                  alphai,
                                                  dc,
                                                  q_in,
                                                  this );
	}
        else if(type_in == "MLS"){
// 	    alphai=3.5; // For the validation triangle, this works better: phi closer to 1.
            this->shapeFun = new ShapeFunctionMLS(supportNodesSize,
                                                  1,
                                                  1, // weight type
                                                  alphai,
                                                  dc,
                                                  this );
	}
	cout << "INFO AT shapeFunSolve IN GaussPoint: (x, y) = "
	    << this->X << ", " << this->Y << endl;
	cout << "\t alphai = " << alphai << ", "
	    << "dc = " << dc << ", "
	    << "q_in = "<< q_in
	    << endl;
	cout << "\t Number of Support Nodes = " << supportNodesSize  << endl;

        shapeFun->calc();
    }
}



void GaussPoint::computeCij( )
{
    //////////////// Calculation of Capacity matrix:
    //////////////// M = rho * Cp * wg * N^T * N * |Jc|
    //////////////// Mij = rho * Cp * wg * Ni * Nj * |Jc| = M(i ,j)
//   cout << mat->getDensity() << " " << mat->getCapacity() << " = density, capacity" << endl;
    int i, j;
//     C.fillIdentity(0.);
    for (i=0; i<supportNodesSize; ++i) {
//         for (j=0; j<supportNodesSize; ++j) {
//             C.writeElement( mat->getDensity() * mat->getCapacity() * weight * shapeFun->getPhi(0,i)
//                             * shapeFun->getPhi(0,j) * std::abs(jacobian), i, j );	
// Lumped matrix:
//          C.addElement( mat->getDensity() * mat->getCapacity() * weight * shapeFun->getPhi(0,i)
//                             * shapeFun->getPhi(0,j) * std::abs(jacobian), i, i );	
//         }
// Faster lumped matrix:
         C.writeElement( mat->getDensity() * mat->getCapacity(supportNodes[i]->getTemp()) * weight * shapeFun->getPhi(0,i)
                            * std::abs(jacobian), i, i );	
// 	  cout << i << "," << j << " = "
// 	       << mat->getDensity() << "*" 
// 	       << mat->getCapacity() << "*" 
// 	       << weight  << "*" 
// 	       << shapeFun->getPhi(0,i) << "*" 
// 	       << shapeFun->getPhi(0,j)  << "*" 
// 	       << std::abs(jacobian)  << " = "
// 	       << C.readElement(i,j) << endl;
    }
//     HRZ lumped matrix [Hinton et al. (1976)] (if implemented, should be in Cell class)

}

void GaussPoint::computeHij( )
{
    int i, j, m, n;
    double avgTemp=0;
    int max_deriv_index = dim+1;
    H.fillIdentity(0.);
    for ( i=0; i<supportNodesSize; ++i ) {
      avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0,i);
    }
    // Hij = wg * grad(N_j) * kappa * grad(N_I) * |Jc|
    for ( i=0; i<supportNodesSize; ++i ) {
        for ( j=0; j<supportNodesSize; ++j ) {
// 	  if(supportNodes[i]->getThermalNumber() == 39){
// 	    cout << "KAPPA in " << i << "," << j << " = " ;
// 	    cout << 0.5*( mat->getKappa(supportNodes[i]->getTemp()) + mat->getKappa(supportNodes[j]->getTemp()) );
// 	    cout << " with nodes temperatures:  ";
// 	    cout << supportNodes[i]->getTemp() << " and ";
// 	    cout << supportNodes[j]->getTemp() << endl;
// 	  }
            for (  m=1; m<max_deriv_index; ++m ) {
//         for ( n=1; n<max_deriv_index; ++n ){
                //           K(2*i + m, 2*j + n) = Kij(m,n);
                H.addElement(weight * (  shapeFun->getPhi(m,i) 
// 					 * mat->getKappa(supportNodes[i]->getTemp()) 
// 					 + .5*mat->getKappa(supportNodes[j]->getTemp())/*.readElement(m,n)*/
					 * mat->getKappa(avgTemp) /*.readElement(m,n)*/
// 					 * 0.5*( mat->getKappa(supportNodes[i]->getTemp()) + mat->getKappa(supportNodes[j]->getTemp()) )/*.readElement(m,n)*/
                                         * shapeFun->getPhi(m,j) ) * std::abs(jacobian),
                             i,
                             j);
// 	  cout << i << "," << j << " = "
// 	       << mat->getDensity() << "*" 
// 	       << mat->getKappa() << "*" 
// 	       << weight  << "*" 
// 	       << shapeFun->getPhi(m,i) << "*" 
// 	       << shapeFun->getPhi(m,j)  << "*" 
// 	       << jacobian  << " = "
// 	       << H.readElement(i,j) << endl;
//         }
            }
//         H.writeElement( weight * mat->getKappa() * jacobian *
//                         (  shapeFun->getPhi(1,i) * shapeFun->getPhi(1,j) +
//                          + shapeFun->getPhi(2,i) * shapeFun->getPhi(2,j)
// 			),
//                                   i,
//                                   j);
        }
    }
//     cout << "GP (" << this->getX() << ", " << this->getY() << ")" << endl;
//     cout << "H = " << H << endl;
}

void GaussPoint::computeQext( LoadThermalBody* loadThermalBody_in )
{
    int i;
    for (i=0; i<supportNodesSize; ++i) {
        // Qi = wg * r * N_I * |Jc|
        Qext.writeElement( weight * shapeFun->getPhi(0,i) * loadThermalBody_in->getLoadThermalBody( this )
                           * std::abs(jacobian), i);
//       M.writeElement( weight * shapeFun->getPhi(0,i) * shapeFun->getPhi(0,j) * jacobian, dim*i+2, dim*j+2 );
    }
}


void GaussPoint::assembleCij( lmx::Matrix<data_type> & globalCapacity )
{
    int i, j;
    for (i=0; i<supportNodesSize; ++i) {
        for (j=0; j<supportNodesSize; ++j) {
            globalCapacity.addElement( C.readElement(i, j),
                                       supportNodes[i]->getThermalNumber(),
                                       supportNodes[j]->getThermalNumber()
                                     );
        }
    }
}

void GaussPoint::assembleHij( lmx::Matrix<data_type> & globalConductivity )
{
    int i, j;
    for (i=0; i<supportNodesSize; ++i) {
        for (j=0; j<supportNodesSize; ++j) {
            globalConductivity.addElement( H.readElement(i, j),
                                           supportNodes[i]->getThermalNumber(),
                                           supportNodes[j]->getThermalNumber()
                                         );
        }
    }
}

void GaussPoint::assembleQext( lmx::Vector<data_type> & globalHeat)
{
    int i;
    for (i=0; i<supportNodesSize; ++i) {
        globalHeat.addElement( Qext.readElement(i),
                               supportNodes[i]->getThermalNumber()
                             );
    }
}


void GaussPoint::gnuplotOutStress( std::ofstream & gptension )
{
    gptension << X << " " << Y << " " << tension(0) << endl;
}


} //Namespace mknix
