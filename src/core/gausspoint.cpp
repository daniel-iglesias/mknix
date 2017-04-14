//-- Licencia --
#include "gausspoint.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctionMLS.h"

#include <simulation/simulation.h>
#include <system/system.h>
#include <system/loadthermalbody.h>

namespace mknix {

GaussPoint::GaussPoint()
{
}


GaussPoint::GaussPoint(int dim_in,
                       double alpha_in,
                       double weight_in,
                       double jacobian_in,
                       Material * mat_in,
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

GaussPoint::GaussPoint(int dim_in,
                       double alpha_in,
                       double weight_in,
                       double jacobian_in,
                       Material * mat_in,
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


void GaussPoint::shapeFunSolve(std::string type_in, double q_in)
{
    q_in = 0.5; // Original RBF
    if (!shapeFun) {
        if (type_in == "RBF") {
//             alphai = 3.5; // For the validation triangle, this works better: phi closer to 1.
            this->shapeFun = new ShapeFunctionRBF(supportNodesSize,
                                                  0,
                                                  0, // RBF type
                                                  alphai,
                                                  dc,
                                                  q_in,
                                                  this);
        }
        else if (type_in == "MLS") {
// 	    alphai=3.5; // For the validation triangle, this works better: phi closer to 1.
            this->shapeFun = new ShapeFunctionMLS(supportNodesSize,
                                                  1,
                                                  1, // weight type
                                                  alphai,
                                                  dc,
                                                  this);
        }
        /*cout << "INFO AT shapeFunSolve IN GaussPoint: (x, y) = "
        << this->X << ", " << this->Y << endl;
        cout << "\t alphai = " << alphai << ", "
        << "dc = " << dc << ", "
        << "q_in = " << q_in
        << endl;
        cout << "\t Number of Support Nodes = " << supportNodesSize << endl;*/

        shapeFun->calc();
    }
}


void GaussPoint::computeCij()
{
    // TODO: not sure why it's needed to be done here too, but for the moment is required for positive validation
    avgTemp = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
    }
    double avgFactor = mat->getDensity() * mat->getCapacity(avgTemp) * weight * std::abs(jacobian);
    // BUG: This is a test for adding lower conductivity layer for y > -1E-5 (MAST-U CFC tiles)
//     if (this->Y > -20E-5) {avgFactor*=0.001;}
    //////////////// Calculation of Capacity matrix:
    //////////////// M = rho * Cp * wg * N^T * N * |Jc|
    //////////////// Mij = rho * Cp * wg * Ni * Nj * |Jc| = M(i ,j)
    int j;
//   cout << mat->getDensity() << " " << mat->getCapacity() << " = density, capacity" << endl;
//     C.reset();
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (j=0; j<supportNodesSize; ++j) {
            C.writeElement( avgFactor * shapeFun->getPhi(0,i) * shapeFun->getPhi(0,j), i, j );
        }
/////////////////////////////////
// Do not use Lumped in ALICIA //
/////////////////////////////////
// Lumped matrix:
//          C.addElement( mat->getDensity() * mat->getCapacity() * weight * shapeFun->getPhi(0,i)
//                             * shapeFun->getPhi(0,j) * std::abs(jacobian), i, i );
//         }
// Faster lumped matrix:
//         C.writeElement(
//                 mat->getDensity() * mat->getCapacity(supportNodes[i]->getTemp()) * weight * shapeFun->getPhi(0, i)
//                 * std::abs(jacobian), i, i);
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

void GaussPoint::computeHij()
{
    int max_deriv_index = dim + 1;
    H.reset();
    avgTemp = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
    }
    double avgFactor = mat->getKappa(avgTemp) * weight * std::abs(jacobian);
      // BUG:
      // This is a test for adding lower conductivity layer for y > -1E-5 (MAST-U CFC tiles)
      // if (this->Y > -142E-6) {avgFactor*=0.08;}
      // if (this->Y > -0.7E-5) {avgFactor*=0.005;}
      // This is a test for adding W higher conductivity layer for y > -1E-5 (Tile 6)
      //
    // Hij = wg * grad(N_j) * kappa * grad(N_I) * |Jc|
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
// 	  if(supportNodes[i]->getThermalNumber() == 39){
// 	    cout << "KAPPA in " << i << "," << j << " = " ;
// 	    cout << 0.5*( mat->getKappa(supportNodes[i]->getTemp()) + mat->getKappa(supportNodes[j]->getTemp()) );
// 	    cout << " with nodes temperatures:  ";
// 	    cout << supportNodes[i]->getTemp() << " and ";
// 	    cout << supportNodes[j]->getTemp() << endl;
// 	  }
            for (auto m = 1; m < max_deriv_index; ++m) {
//         for ( n=1; n<max_deriv_index; ++n ){
                //           K(2*i + m, 2*j + n) = Kij(m,n);
                H.addElement( (shapeFun->getPhi(m, i)
                                       // 					 * mat->getKappa(supportNodes[i]->getTemp())
                                       // 					 + .5*mat->getKappa(supportNodes[j]->getTemp())/*.readElement(m,n)*/
//                                        * mat->getKappa(avgTemp) /*.readElement(m,n)*/
                                       // 					 * 0.5*( mat->getKappa(supportNodes[i]->getTemp()) + mat->getKappa(supportNodes[j]->getTemp()) )/*.readElement(m,n)*/
                                * shapeFun->getPhi(m, j)) * avgFactor,
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

void GaussPoint::computeQext(LoadThermalBody * loadThermalBody_in)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        // Qi = wg * r * N_I * |Jc|
        Qext.writeElement(weight * shapeFun->getPhi(0, i) * loadThermalBody_in->getLoadThermalBody(this)
                          * std::abs(jacobian), i);
//       M.writeElement( weight * shapeFun->getPhi(0,i) * shapeFun->getPhi(0,j) * jacobian, dim*i+2, dim*j+2 );
    }
}


void GaussPoint::assembleCij(lmx::Matrix<data_type>& globalCapacity)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            globalCapacity.addElement(C.readElement(i, j),
                                      supportNodes[i]->getThermalNumber(),
                                      supportNodes[j]->getThermalNumber()
            );
        }
    }
}

void GaussPoint::assembleCijWithMap(data_type *globalCapacity,
                                    int *matrixMap,
                                    int num_nodes)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
           //we recover the point position in the sparse matrix
           int rowNode = supportNodes[i]->getThermalNumber();
           int colNode = supportNodes[j]->getThermalNumber();
           int positionGlobal = rowNode + colNode * num_nodes;
          // int debugPos = matrixMap[positionGlobal];
           int mypos = matrixMap[positionGlobal];
           //we add the value directly into position
            globalCapacity[mypos] += C.readElement(i, j);
        }
    }
}

void GaussPoint::presenceCij(int* presenceMatrix, int num_nodes)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            int rowNode = supportNodes[i]->getThermalNumber();
            int colNode = supportNodes[j]->getThermalNumber();
            int positionGlobal = rowNode + colNode * num_nodes;//col storage format
            presenceMatrix[positionGlobal] = 1;
        }
    }
}

void GaussPoint::assembleHij(lmx::Matrix<data_type>& globalConductivity)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            globalConductivity.addElement(H.readElement(i, j),
                                          supportNodes[i]->getThermalNumber(),
                                          supportNodes[j]->getThermalNumber()
            );
        }
    }
}

void GaussPoint::assembleHijWithMap(data_type *globalConductivity,
                                    int *matrixMap,
                                    int num_nodes)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
          //we recover the point position in the sparse matrix
          int rowNode = supportNodes[i]->getThermalNumber();
          int colNode = supportNodes[j]->getThermalNumber();
          int positionGlobal = rowNode + colNode * num_nodes;
          int mypos = matrixMap[positionGlobal];
          //we add the value directly into position
          globalConductivity[mypos] += H.readElement(i, j);
        }
    }
}

void GaussPoint::assembleQext(lmx::Vector<data_type>& globalHeat)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        globalHeat.addElement(Qext.readElement(i),
                              supportNodes[i]->getThermalNumber()
        );
    }
}

double  GaussPoint::getNodePhi(int deriv, int node)
{//for SOA purpouses
  return shapeFun->getPhi(deriv, node);
}

double  GaussPoint::getWeight()
{//for SOA purpouses
  return weight;
}

/*void GaussPoint::mapThermalNumber(std::vector<int> &_locaThermalNumbers m)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        int positionGlobal = supportNodes[i]->getThermalNumber;
        _locaThermalNumbers[positionGlobal] = 1;
        );
    }
}*/


void GaussPoint::gnuplotOutStress(std::ofstream& gptension)
{
    gptension << X << " " << Y << " " << tension(0) << endl;
}


} //Namespace mknix
