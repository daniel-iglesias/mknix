//-- Licencia --
#include "LMX/lmx.h"
#include "cell.h"
#include "material.h"
#include "gausspoint.h"
#include "node.h"

#include <system/loadthermalbody.h>

namespace mknix
{

/**
 * @brief Default constructor for Cell.
 */
Cell::Cell()
{
}


/**
 * @brief Constructs a Cell with a material, formulation type, influence radius factor, and number of Gauss points.
 * @param material_in Reference to the material assigned to this cell.
 * @param formulation_in String identifying the meshfree formulation (e.g. "EFG", "RPIM").
 * @param alpha_in Influence radius scaling factor.
 * @param nGPoints_in Number of integration Gauss points per direction.
 */
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

/**
 * @brief Destructor. Releases all dynamically allocated Gauss points.
 */
Cell::~Cell()
{
    for (auto& point : gPoints)
    {
        delete point;
    }
    /*
    for (auto& point : gPoints_MC) {
        delete point;
    }
    */
}

/**
 * @brief Reassigns the material for all Gauss points if the cell belongs to a thin layer.
 *
 * Detects whether the cell is within a layer of the given thickness by checking if any
 * pair of body-points is closer than the thickness. If so, the material is replaced.
 *
 * @param newMat The replacement material.
 * @param thickness Maximum inter-node distance that qualifies the cell as part of the layer.
 * @return True if the material was changed, false otherwise.
 */
bool Cell::setMaterialIfLayer(Material& newMat, double thickness)
{
    bool changed(0);
    // First we check if the minimum distance between nodes is less than the thickness
    //   This assumes that the layer is composed by the smallest elements of the mesh.
    for (auto& point1 : bodyPoints)
    {
        for (auto& point2 : bodyPoints)
        {
            if(point1 != point2)  // Avoid comparing a point with itself
            {
//                 cout << point1->distance(*point2) << endl;
                if (point1->distance(*point2) < thickness)
                {

                    // Given the case, we iterate in the gPoints to change the Material
                    for (auto& gPoint : gPoints)
                    {
                        gPoint->setMaterial( newMat );
                    }
                    for (auto& gPoint : gPoints_MC)
                    {
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
/**
 * @brief Initializes the cell by finding support nodes for each Gauss point (meshfree formulations).
 * @param nodes_in Vector of domain nodes used to build the support neighbourhood.
 */
void Cell::initialize(std::vector<Node *>& nodes_in)
{
    // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
    for (auto& point : gPoints)
    {
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


/**
 * @brief Computes and stores shape function values at each Gauss point using the selected formulation.
 */
void Cell::computeShapeFunctions()
{
    for (auto& point : gPoints)
    {
        if (formulation == "RPIM")
        {
            point->shapeFunSolve("RBF", 1.03);
        }
        else if (formulation == "EFG")
        {
            point->shapeFunSolve("MLS", 1.03);
        }
    }
}


/**
 * @brief Computes the local thermal capacity (heat capacity) contribution at each Gauss point.
 */
void Cell::computeCapacityGaussPoints()
{
    for (auto& point : gPoints_MC)
    {
        point->computeCij();
    }
}

/**
 * @brief Assembles Gauss-point capacity contributions into the global thermal capacity matrix.
 * @param globalCapacity Global thermal capacity matrix to be updated.
 */
void Cell::assembleCapacityGaussPoints(lmx::Matrix<data_type>& globalCapacity)
{
    for (auto& point : gPoints_MC)
    {
        point->assembleCij(globalCapacity);
    }
}


/**
 * @brief Computes the local thermal conductivity contribution at each Gauss point.
 */
void Cell::computeConductivityGaussPoints()
{
    for (auto& point : gPoints)
    {
        point->computeHij();
    }
}

/**
 * @brief Assembles Gauss-point conductivity contributions into the global conductivity matrix.
 * @param globalConductivity Global thermal conductivity matrix to be updated.
 */
void Cell::assembleConductivityGaussPoints(lmx::Matrix<data_type>& globalConductivity)
{
    for (auto& point : gPoints)
    {
        point->assembleHij(globalConductivity);
    }
}


/**
 * @brief Computes the external thermal load vector contribution at each Gauss point.
 * @param loadThermalBody_in Vector of body thermal loads applied to this cell.
 */
void Cell::computeQextGaussPoints(std::vector<LoadThermalBody*> loadThermalBody_in)
{
    for (auto& point : gPoints)
    {
        point->computeQext(loadThermalBody_in);
    }
}

/**
 * @brief Assembles external thermal load contributions from Gauss points into the global heat vector.
 * @param globalQext Global external heat load vector to be updated.
 */
void Cell::assembleQextGaussPoints(lmx::Vector<data_type>& globalQext)
{
    for (auto& point : gPoints)
    {
        point->assembleQext(globalQext);
    }
}


/**
 * @brief Computes the local mass matrix contribution at each Gauss point.
 */
void Cell::computeMGaussPoints()
{
    for (auto& point : gPoints_MC)
    {
        point->computeMij();
    }
}


/**
 * @brief Assembles Gauss-point mass contributions into the global mass matrix.
 * @param globalMass Global mass matrix to be updated.
 */
void Cell::assembleMGaussPoints(lmx::Matrix<data_type>& globalMass)
{
    for (auto& point : gPoints_MC)
    {
        point->assembleMij(globalMass);
    }
}


/**
 * @brief Computes the linear internal force vector contribution at each Gauss point.
 */
void Cell::computeFintGaussPoints()
{
    for (auto& point : gPoints)
    {
        point->computeFint();
    }
}


/**
 * @brief Computes the nonlinear internal force vector contribution at each Gauss point.
 */
void Cell::computeNLFintGaussPoints()
{
    for (auto& point : gPoints)
    {
        point->computeNLFint();
    }
}


/**
 * @brief Assembles internal force contributions from Gauss points into the global internal force vector.
 * @param globalFint Global internal force vector to be updated.
 */
void Cell::assembleFintGaussPoints(lmx::Vector<data_type>& globalFint)
{
    for (auto& point : gPoints)
    {
        point->assembleFint(globalFint);
    }
}


/**
 * @brief Computes the external force vector contribution at each Gauss point.
 */
void Cell::computeFextGaussPoints()
{
    for (auto& point : gPoints_MC)
    {
        point->computeFext();
    }
}


/**
 * @brief Assembles external force contributions from Gauss points into the global external force vector.
 * @param globalFext Global external force vector to be updated.
 */
void Cell::assembleFextGaussPoints(lmx::Vector<data_type>& globalFext)
{
    for (auto& point : gPoints_MC)
    {
        point->assembleFext(globalFext);
    }
}


/**
 * @brief Computes the linear tangent (stiffness) matrix contribution at each Gauss point.
 */
void Cell::computeKGaussPoints()
{
    for (auto& point : gPoints)
    {
        point->computeKij();
    }
}


/**
 * @brief Computes the nonlinear tangent matrix contribution at each Gauss point.
 */
void Cell::computeNLKGaussPoints()
{
    for (auto& point : gPoints)
    {
        point->computeNLKij();
    }
}


/**
 * @brief Assembles tangent matrix contributions from Gauss points into the global tangent matrix.
 * @param globalTangent Global tangent (stiffness) matrix to be updated.
 */
void Cell::assembleKGaussPoints(lmx::Matrix<data_type>& globalTangent)
{
    for (auto& point : gPoints)
    {
        point->assembleKij(globalTangent);
    }
}


/**
 * @brief Computes linear stress at each Gauss point and assembles the result into a body stress vector.
 * @param globalStress Body-level stress resultant vector to be updated.
 * @param firstNode Global index of the first node in this body.
 */
void Cell::assembleRGaussPoints(lmx::Vector<data_type>& globalStress,
                                int firstNode
                               )
{
    for (auto& point : gPoints)
    {
        point->computeStress();
        point->assembleRi(globalStress, firstNode);
    }
}


/**
 * @brief Computes nonlinear (large-deformation) stress at each Gauss point and assembles the result.
 * @param globalStress Body-level stress resultant vector to be updated.
 * @param firstNode Global index of the first node in this body.
 */
void Cell::assembleNLRGaussPoints(lmx::Vector<data_type>& globalStress,
                                  int firstNode
                                 )
{
    for (auto& point : gPoints)
    {
        point->computeNLStress();
        point->assembleRi(globalStress, firstNode);
    }
}


/**
 * @brief Accumulates the potential (gravitational/body-force) energy across all Gauss points.
 * @param q Current global displacement state vector.
 * @return Total potential energy contribution from this cell.
 */
double Cell::calcPotentialEGaussPoints(const lmx::Vector<data_type>& q)
{
    double potentialEnergy = 0;

    for (auto& point : gPoints)
    {
        potentialEnergy += point->calcPotentialE(q);
    }
    return potentialEnergy;
}


/**
 * @brief Accumulates the kinetic energy across all Gauss points.
 * @param qdot Current global velocity state vector.
 * @return Total kinetic energy contribution from this cell.
 */
double Cell::calcKineticEGaussPoints(const lmx::Vector<data_type>& qdot)
{
    double kineticEnergy = 0;

    for (auto& point : gPoints)
    {
        kineticEnergy += point->calcKineticE(qdot);
    }
    return kineticEnergy;
}


/**
 * @brief Accumulates the elastic strain energy across all Gauss points.
 * @return Total elastic energy contribution from this cell.
 */
double Cell::calcElasticEGaussPoints()
{
    double elasticEnergy = 0;

    for (auto& point : gPoints)
    {
        elasticEnergy += point->calcElasticE();
    }
    return elasticEnergy;
}


/**
 * @brief Writes the node numbers of this cell's body points to an output file.
 * @param outfile Pointer to the open output file stream.
 */
void Cell::outputConnectivityToFile(std::ofstream * outfile)
{
    *outfile << "\t\t\t";
    for (auto& point : bodyPoints)
    {
        *outfile << point->getNumber() << " ";
    }
    *outfile << std::endl;
}


/**
 * @brief Outputs Gauss-point stress values to a file suitable for gnuplot visualization.
 * @param gptension Output file stream for stress/tension data.
 */
void Cell::gnuplotOutStress(std::ofstream& gptension)
{
    int counter = 0;
    for (auto& point : gPoints)
    {
        ++counter;
        point->gnuplotOutStress(gptension);
        if (counter % 4 == 0) gptension << endl;
    }
}

/**
 * @brief Returns the global node numbers of all body points belonging to this cell.
 * @return Vector of node numbers.
 */
std::vector<int> Cell::getNodeNumbers()
{
    std::vector<int> nodeNumbers;
    for (auto& point : bodyPoints)
    {
        nodeNumbers.push_back( point->getNumber() );
    }
    return nodeNumbers;
}


} //Namespace mknix
