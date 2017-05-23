//-- Licencia --
#include "LMX/lmx.h"
#include "cell.h"
#include "material.h"
#include "gausspoint.h"
#include "node.h"

#include <system/loadthermalbody.h>

namespace mknix {

Cell::Cell()
{
}


Cell::Cell(Material& material_in,
           std::string formulation_in,
           double alpha_in,
           int nGPoints_in)
        : mat(&material_in)
        , formulation(formulation_in)
        , alpha(alpha_in)
        , nGPoints(nGPoints_in)
        , dc(0)
{
}

Cell::~Cell()
{
    for (auto& point : gPoints) {
        delete point;
    }
    /*
    for (auto& point : gPoints_MC) {
        delete point;
    }
    */
}

int Cell::getMaterialId(){return mat->getMaterialId();}

double Cell::getNodePhi(int gp, int deriv, int node)
{
  return gPoints[gp]->getNodePhi(deriv, node);
}

double Cell::getCellPhi(int gp, int deriv, int node)
{
  return gPoints_MC[gp]->getNodePhi(deriv, node);
}

double Cell::getWeightMC(int gp)
{
  return gPoints_MC[gp]->getWeight();
}

double Cell::getWeight(int gp)
{
  return gPoints[gp]->getWeight();
}

double Cell::getJacobian()
{
  return jacobian;
}

double Cell::getJacobianMC(int gp)
{
  return gPoints_MC[gp]->getJacobian();
}
double Cell::getJacobianP(int gp)
{
  return gPoints[gp]->getJacobian();
}

int Cell::getSupportSizeMC(){gPoints_MC[0]->getSupportSize();}
int Cell::getSupportSize(){gPoints[0]->getSupportSize();}

bool Cell::setMaterialIfLayer(Material& newMat, double thickness)
{
    bool changed(0);
    // First we check if the minimum distance between nodes is less than the thickness
    //   This assumes that the layer is composed by the smallest elements of the mesh.
    for (auto& point1 : bodyPoints){
        for (auto& point2 : bodyPoints){
            if(point1 != point2){ // Avoid comparing a point with itself
//                 cout << point1->distance(*point2) << endl;
                if (point1->distance(*point2) < thickness){

                    // Given the case, we iterate in the gPoints to change the Material
                    for (auto& gPoint : gPoints) {
                        gPoint->setMaterial( newMat );
                    }
                    for (auto& gPoint : gPoints_MC) {
                        gPoint->setMaterial( newMat );
                    }
                    changed=1;
                }
            }
        }
    }
    return changed;
}


// Only for Meshfree Cells, function is specialized for FEM elements
void Cell::initialize(std::vector<Node *>& nodes_in)
{
    // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
    for (auto& point : gPoints) {
        gPoints_MC.push_back(point); // use same GP for all matrices
        point->findSupportNodes(nodes_in);
    }
    // Set the dc and alpha parameteres for cell nodes. This way, the values will
    // be greater than zero only for meshfree nodes which need shapefunctions to be
    // calculated. The type of shape function is also set here.
//   std::vector<Node*>::iterator it_p;
//   for ( it_p = nodes.begin();
//         it_p != nodes.end();
//         ++it_p)
//   {
//     (*it_p)->setAlphai( alpha );
//     (*it_p)->setDc( dc );
//     if( formulation == "RPIM" )
//       (*it_p)->setShapeFunType( "RBF" );
//     else if( formulation == "EFG" )
//       (*it_p)->setShapeFunType( "MLS" );
//   }

}


void Cell::computeShapeFunctions()
{
    for (auto& point : gPoints) {
        if (formulation == "RPIM") {
            point->shapeFunSolve("RBF", 1.03);
        } else if (formulation == "EFG") {
            point->shapeFunSolve("MLS", 1.03);
        }
    }
}


void Cell::computeCapacityGaussPoints()
{
    for (auto& point : gPoints_MC) {
        point->computeCij();
    }
}
std::vector<double> Cell::getShapeCij(){
  std::cout << "Cell::getShapeCij has " << gPoints_MC.size()<< " gausspoints"<<std::endl;
  return gPoints_MC[0]->getShapeCij();}
  std::vector<double> Cell::getCij(){return gPoints_MC[0]->getCij();}//FOR DEBUG ONLY
  std::vector<double> Cell::getTempsCij(){return gPoints_MC[0]->getTemps();}//FOR DEBUG ONLY
  double Cell::getCFactor(){return gPoints_MC[0]->getCFactor();}//FOR DEBUG ONLY

int Cell::getNumPoints_MC(){return gPoints_MC.size();}
int Cell::getNumPoints(){return gPoints.size();}

void Cell::presenceCapacityGaussPoints(int* presence_matrix, int number_nodes)
{
  for (auto& point : gPoints_MC) {
      point->presenceCij(presence_matrix, number_nodes);
  }
}

void Cell::assembleCapacityGaussPoints(lmx::Matrix<data_type>& globalCapacity)
{
    for (auto& point : gPoints_MC) {
        point->assembleCij(globalCapacity);
    }
}

void Cell::assembleCapacityGaussPointsWithMap(data_type *globalCapacity,
                                              uint *matrixMap,
                                              int number_nodes)
{
    for (auto& point : gPoints_MC) {
        point->assembleCijWithMap(globalCapacity,
                                  matrixMap,
                                  number_nodes);
    }
}


void Cell::computeConductivityGaussPoints()
{
    for (auto& point : gPoints) {
        point->computeHij();
    }
}
std::vector<double> Cell::getShapeHij(){return gPoints[0]->getShapeHij();}
std::vector<double> Cell::getHij(){return gPoints[0]->getHij();}
std::vector<double> Cell::getTempsHij(){return gPoints[0]->getTemps();}//FOR DEBUG ONLY
double Cell::getHFactor(){return gPoints[0]->getHFactor();}//FOR DEBUG ONLY

void Cell::presenceConductivityGaussPoints(int* presence_matrix, int number_nodes)
{
  for (auto& point : gPoints) {
      point->presenceHij(presence_matrix, number_nodes);
  }
}

void Cell::assembleConductivityGaussPoints(lmx::Matrix<data_type>& globalConductivity)
{
    for (auto& point : gPoints) {
        point->assembleHij(globalConductivity);
    }
}

void Cell::assembleConductivityGaussPointsWithMap(data_type *globalConductivity,
                                                  uint *matrixMap,
                                                  int number_nodes)
{
    for (auto& point : gPoints) {
        point->assembleHijWithMap(globalConductivity,
                                  matrixMap,
                                  number_nodes);
    }
}


void Cell::computeQextGaussPoints(LoadThermalBody * loadThermalBody_in)
{
    for (auto& point : gPoints) {
        point->computeQext(loadThermalBody_in);
    }
}

void Cell::assembleQextGaussPoints(lmx::Vector<data_type>& globalQext)
{
    for (auto& point : gPoints) {
        point->assembleQext(globalQext);
    }
}


void Cell::computeMGaussPoints()
{
    for (auto& point : gPoints_MC) {
        point->computeMij();
    }
}


void Cell::assembleMGaussPoints(lmx::Matrix<data_type>& globalMass)
{
    for (auto& point : gPoints_MC) {
        point->assembleMij(globalMass);
    }
}


void Cell::computeFintGaussPoints()
{
    for (auto& point : gPoints) {
        point->computeFint();
    }
}


void Cell::computeNLFintGaussPoints()
{
    for (auto& point : gPoints) {
        point->computeNLFint();
    }
}


void Cell::assembleFintGaussPoints(lmx::Vector<data_type>& globalFint)
{
    for (auto& point : gPoints) {
        point->assembleFint(globalFint);
    }
}


void Cell::computeFextGaussPoints()
{
    for (auto& point : gPoints_MC) {
        point->computeFext();
    }
}


void Cell::assembleFextGaussPoints(lmx::Vector<data_type>& globalFext)
{
    for (auto& point : gPoints_MC) {
        point->assembleFext(globalFext);
    }
}


void Cell::computeKGaussPoints()
{
    for (auto& point : gPoints) {
        point->computeKij();
    }
}


void Cell::computeNLKGaussPoints()
{
    for (auto& point : gPoints) {
        point->computeNLKij();
    }
}


void Cell::assembleKGaussPoints(lmx::Matrix<data_type>& globalTangent)
{
    for (auto& point : gPoints) {
        point->assembleKij(globalTangent);
    }
}


void Cell::assembleRGaussPoints(lmx::Vector<data_type>& globalStress,
                                int firstNode
)
{
    for (auto& point : gPoints) {
        point->computeStress();
        point->assembleRi(globalStress, firstNode);
    }
}


void Cell::assembleNLRGaussPoints(lmx::Vector<data_type>& globalStress,
                                  int firstNode
)
{
    for (auto& point : gPoints) {
        point->computeNLStress();
        point->assembleRi(globalStress, firstNode);
    }
}


double Cell::calcPotentialEGaussPoints(const lmx::Vector<data_type>& q)
{
    double potentialEnergy = 0;

    for (auto& point : gPoints) {
        potentialEnergy += point->calcPotentialE(q);
    }
    return potentialEnergy;
}


double Cell::calcKineticEGaussPoints(const lmx::Vector<data_type>& qdot)
{
    double kineticEnergy = 0;

    for (auto& point : gPoints) {
        kineticEnergy += point->calcKineticE(qdot);
    }
    return kineticEnergy;
}


double Cell::calcElasticEGaussPoints()
{
    double elasticEnergy = 0;

    for (auto& point : gPoints) {
        elasticEnergy += point->calcElasticE();
    }
    return elasticEnergy;
}


void Cell::outputConnectivityToFile(std::ofstream * outfile)
{
    *outfile << "\t\t\t";
    for (auto& point : bodyPoints) {
        *outfile << point->getNumber() << " ";
    }
    *outfile << std::endl;
}


void Cell::gnuplotOutStress(std::ofstream& gptension)
{
    int counter;
    for (auto& point : gPoints) {
        ++counter;
        point->gnuplotOutStress(gptension);
        if (counter % 4 == 0) gptension << endl;
    }
}

void Cell::mapThermalNodesMC(int* thermalMap,
                            int supportNodeSize,
                            int cell_index )
{
  int np = gPoints_MC.size();
  for (auto& point : gPoints_MC) {
    int position = cell_index * np * supportNodeSize;
    point->mapThermalNumbers(thermalMap, position);
  }
}

void Cell::mapThermalNodes(int* thermalMap,
                           int supportNodeSize,
                           int cell_index)
{
  int np = gPoints.size();
  for (auto& point : gPoints) {
    int position = cell_index * np * supportNodeSize;
    point->mapThermalNumbers(thermalMap,position);
  }
}

void Cell::mapNodesMC(uint* Map,
                      int supportNodeSize,
                      int cell_index,
                      int total_nodes)
{
  int np = gPoints_MC.size();
  for (auto& point : gPoints_MC) {
    int position = cell_index * np * supportNodeSize * supportNodeSize;
    point->mapNodeNumbers(Map,
                          position,
                          total_nodes);
  }
}

void Cell::mapNodes(uint* Map,
                    int supportNodeSize,
                    int cell_index,
                    int total_nodes)
{
  int np = gPoints.size();
  for (auto& point : gPoints) {
    int position = cell_index * np * supportNodeSize * supportNodeSize;
    point->mapNodeNumbers(Map,
                          position,
                          total_nodes);
  }
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// templating part /////////////////////////////////

/*template std::vector<double> Cell::getCij<double>();
template std::vector<float> Cell::getCij<float>();
template std::vector<double> Cell::getHij<double>();
template std::vector<float> Cell::getHij<float>();*/

} //Namespace mknix
