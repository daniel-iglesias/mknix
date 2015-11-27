//-- Licencia --
#include "LMX/lmx.h"
#include "cellboundary.h"
#include "gausspointboundary.h"
#include "node.h"
#include "loadthermalboundary1D.h"

namespace mknix {

CellBoundary::CellBoundary()
{
}


CellBoundary::CellBoundary(std::string formulation_in,
                           double alpha_in,
                           int nGPoints_in)
        : formulation(formulation_in)
        , alpha(alpha_in)
        , nGPoints(nGPoints_in)
        , dc(0)
{
}

CellBoundary::~CellBoundary()
{
    for (auto& point : gPoints) {
        delete point;
    }
}

// Only for Meshfree CellBoundarys, function is specialized for FEM in each cell type
void CellBoundary::initialize(std::vector<Node *>& nodes_in)
{
    if (formulation == "RPIM" || formulation == "EFG") {
        // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
        for (auto& point : gPoints) {
            point->findSupportNodes(nodes_in);
        }
    }
}


void CellBoundary::computeShapeFunctions()
{
    for (auto& point : gPoints) {
        if (formulation == "RPIM") {
            point->shapeFunSolve("RBF", 1.03);
        } else if (formulation == "EFG") {
            point->shapeFunSolve("MLS", 1.03);
        }
    }
}


void CellBoundary::computeQextGaussPoints(LoadThermalBoundary1D * loadThermalBoundary1D_in)
{
    for (auto& point : gPoints) {
        point->computeQext(loadThermalBoundary1D_in);
    }
}

void CellBoundary::assembleQextGaussPoints(lmx::Vector<data_type>& globalQext)
{
    for (auto& point : gPoints) {
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
