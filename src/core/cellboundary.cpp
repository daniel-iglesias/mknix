//-- Licencia --
#include "LMX/lmx.h"
#include "node.h"
#include "cellboundary.h"
#include "gausspointboundary.h"

#include <system/loadthermalboundary1D.h>


namespace mknix
{

/**
 * @brief Default constructor for CellBoundary.
 */
CellBoundary::CellBoundary()
{
}


/**
 * @brief Constructs a CellBoundary with formulation, influence radius factor, and number of Gauss points.
 * @param formulation_in String identifying the formulation (e.g. "EFG", "RPIM").
 * @param alpha_in Influence radius scaling factor.
 * @param nGPoints_in Number of integration Gauss points.
 */
CellBoundary::CellBoundary(std::string formulation_in,
                           double alpha_in,
                           int nGPoints_in)
    : formulation(formulation_in)
    , alpha(alpha_in)
    , nGPoints(nGPoints_in)
    , dc(0)
{
}

/**
 * @brief Destructor. Releases all dynamically allocated boundary Gauss points.
 */
CellBoundary::~CellBoundary()
{
    for (auto& point : gPoints)
    {
        delete point;
    }
}

// Only for Meshfree CellBoundarys, function is specialized for FEM in each cell type
/**
 * @brief Initializes boundary Gauss points by finding support nodes (meshfree formulations only).
 * @param nodes_in Vector of domain nodes used to build the support neighbourhood.
 */
void CellBoundary::initialize(std::vector<Node *>& nodes_in)
{
    if (formulation == "RPIM" || formulation == "EFG")
    {
        // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
        for (auto& point : gPoints)
        {
            point->findSupportNodes(nodes_in);
        }
    }
}


/**
 * @brief Computes shape function values at each boundary Gauss point.
 */
void CellBoundary::computeShapeFunctions()
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
 * @brief Computes the external thermal boundary load contribution at each Gauss point.
 * @param loadThermalBoundary1D_in Pointer to the 1D boundary thermal load object.
 */
void CellBoundary::computeQextGaussPoints(LoadThermalBoundary1D * loadThermalBoundary1D_in)
{
    for (auto& point : gPoints)
    {
        point->computeQext(loadThermalBoundary1D_in);
    }
}

/**
 * @brief Assembles boundary thermal load contributions from Gauss points into the global heat vector.
 * @param globalQext Global external heat load vector to be updated.
 */
void CellBoundary::assembleQextGaussPoints(lmx::Vector<data_type>& globalQext)
{
    for (auto& point : gPoints)
    {
        point->assembleQext(globalQext);
    }
}

// void CellBoundary::outputConnectivityToFile(std::ofstream* outfile)
// {
//   std::vector< Point* >::iterator it_points;
//   *outfile << "\t\t\t";
//   for(it_points=bodyPoints.begin();
//       it_points!=bodyPoints.end();
//       ++it_points){
//     *outfile << (*it_points)->getNumber() << " ";
//   }
//   *outfile << std::endl;
// }
//
//
// void CellBoundary::gnuplotOutStress( std::ofstream & gptension )
// {
//     int counter;
//     for(std::vector<GaussPoint*>::iterator it=gPoints.begin();
//             it!=gPoints.end();
//             ++it)
//     {
//         ++counter;
//         (*it)->gnuplotOutStress( gptension );
//         if (counter%4 == 0) gptension << endl;
//     }
// }

} //Namespace mknix
