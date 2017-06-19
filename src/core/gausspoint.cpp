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
  //std::cout << "REALLY?";
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
   //C.reset();//aparently this doesnt matter
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

std::vector<double> GaussPoint::getShapeCij(){
  std::vector<double> myShapeCij(supportNodesSize *supportNodesSize);
  for (auto i = 0; i < supportNodesSize; ++i)
      for (auto j = 0; j<supportNodesSize; ++j)
        myShapeCij[i*supportNodesSize+j] =  shapeFun->getPhi(0,i) * shapeFun->getPhi(0,j);
   return myShapeCij;
}

std::vector<double> GaussPoint::getCij(){
  std::vector<double> myCij(supportNodesSize *supportNodesSize);
  for (auto i = 0; i < supportNodesSize; ++i)
      for (auto j = 0; j<supportNodesSize; ++j)
        myCij[i*supportNodesSize+j] =  C.readElement(i, j);
   return myCij;
}
std::vector<double> GaussPoint::getTemps(){
  std::vector<double> myTemps(supportNodesSize);
  for (auto i = 0; i < supportNodesSize; ++i)
        myTemps[i] =  supportNodes[i]->getTemp();
   return myTemps;
}

double GaussPoint::getCFactor(){
  avgTemp = 0;
  for (auto i = 0u; i < supportNodesSize; ++i) {
      avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
  }
  double avgFactor = mat->getDensity() * mat->getCapacity(avgTemp) * weight * std::abs(jacobian);
  //double avgFactor =  mat->getCapacity(avgTemp);OK!
  //double avgFactor =  weight ;OK!
  //double avgFactor =  std::abs(jacobian);
  return avgFactor;
}

void GaussPoint::computeHij()
{
  //std::cout << "WTF ";
    int max_deriv_index = dim + 1;
    H.reset();
    avgTemp = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
    }
    double avgFactor = mat->getKappa(avgTemp) * weight * std::abs(jacobian);

    for (auto i = 0u; i < supportNodesSize; ++i) {
        for (auto j = 0u; j < supportNodesSize; ++j) {
            for (auto m = 1; m < max_deriv_index; ++m) {
                H.addElement( (shapeFun->getPhi(m, i) * shapeFun->getPhi(m, j)) * avgFactor,i,j);
            }
        }
    }

}

double GaussPoint::getHFactor()
{
    int max_deriv_index = dim + 1;
    avgTemp = 0;
    for (auto i = 0u; i < supportNodesSize; ++i) {
        avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
    }
    double avgFactor = mat->getKappa(avgTemp) * weight * std::abs(jacobian);
    return avgFactor;
}

std::vector<double> GaussPoint::getShapeHij(){
  int max_deriv_index = dim + 1;
  std::vector<double> myShapeHij(supportNodesSize *supportNodesSize);
  for (auto i = 0; i < supportNodesSize; ++i){
      for (auto j = 0; j<supportNodesSize; ++j){
         for (auto m = 1; m < max_deriv_index; ++m){ myShapeHij[i*supportNodesSize+j] +=  shapeFun->getPhi(m, i)*shapeFun->getPhi(m, j);}
      }
  }
  return myShapeHij;
}

std::vector<double> GaussPoint::getHij(){
  std::vector<double> myHij(supportNodesSize *supportNodesSize);
  /*int max_deriv_index = dim + 1;
  avgTemp = 0;
  for (auto i = 0u; i < supportNodesSize; ++i) {
      avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
  }
  double avgFactor = mat->getKappa(avgTemp) * weight * std::abs(jacobian);

  for (auto i = 0u; i < supportNodesSize; ++i) {
      for (auto j = 0u; j < supportNodesSize; ++j) {
          for (auto m = 1; m < max_deriv_index; ++m) {
              myHij[i *supportNodesSize + j] = shapeFun->getPhi(m, i)* shapeFun->getPhi(m, j) * avgFactor;
          }
      }
  }*/
  for (auto i = 0; i < supportNodesSize; ++i){
      for (auto j = 0; j<supportNodesSize; ++j){
        myHij[i*supportNodesSize+j] =  H.readElement(i,j);
      }
    }
   return myHij;
}

void GaussPoint::computeQext(LoadThermalBody * loadThermalBody_in)
{
  //std::cout << "void GaussPoint::computeQext" <<std::endl;
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
                                    uint *matrixMap,
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
                                    uint *matrixMap,
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

void GaussPoint::presenceHij(int* presenceMatrix, int num_nodes)
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

void GaussPoint::assembleQext(lmx::Vector<data_type>& globalHeat)
{
  //std::cout << " GaussPoint::assembleQext" <<std::endl; This is used!
    for (auto i = 0u; i < supportNodesSize; ++i) {
        globalHeat.addElement(Qext.readElement(i),
                              supportNodes[i]->getThermalNumber()
        );
    }
}
void GaussPoint::assembleQext(VectorX<data_type>& globalHeat)
{
  //std::cout << " GaussPoint::assembleQext" <<std::endl; This is used!
    for (auto i = 0u; i < supportNodesSize; ++i) {
        globalHeat[supportNodes[i]->getThermalNumber()] += Qext.readElement(i);
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

double  GaussPoint::getJacobian()
{//for SOA purpouses
  return jacobian;
}

int  GaussPoint::getSupportSize()
{//for SOA purpouses
  return supportNodesSize;
}

void GaussPoint::mapThermalNumbers(int* thermalMap, int my_position)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        thermalMap[my_position + i] = supportNodes[i]->getThermalNumber();
    }
}
void GaussPoint::mapNodeNumbers(uint* thermalMap,
                                int my_position,
                                int total_nodes)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
       int row = supportNodes[i]->getThermalNumber();
       for (auto j = 0u; j < supportNodesSize; ++j) {
         int col = supportNodes[j]->getThermalNumber();
         int out_id = my_position + i * supportNodesSize + j;
         thermalMap[out_id] = row * total_nodes + col;
      }
    }
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

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// templating part /////////////////////////////////

/*template std::vector<double> GaussPoint::getpCij<double>();
template std::vector<float> GaussPoint::getpCij<float>();
template std::vector<double> GaussPoint::getHij<double>();
template std::vector<float> GaussPoint::getHij<float>();*/

} //Namespace mknix
