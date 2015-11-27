//-- Licencia --
#include "gausspointboundary.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctionMLS.h"
#include "shapefunctionlinear-x.h"
#include "simulation.h"
#include "system.h"
#include "loadthermalboundary1D.h"

namespace mknix {

GaussPointBoundary::GaussPointBoundary()
{
}


GaussPointBoundary::GaussPointBoundary(int dim_in,
                                       double alpha_in,
                                       double weight_in,
                                       double jacobian_in,
                                       int num_in,
                                       double coor_x,
                                       double dc_in
)
        : Point(dim_in, num_in, coor_x, 0., 0., jacobian_in, alpha_in, dc_in)
        , weight(weight_in)
{
}

GaussPointBoundary::GaussPointBoundary(int dim_in,
                                       double alpha_in,
                                       double weight_in,
                                       double jacobian_in,
                                       int num_in,
                                       double coor_x,
                                       double coor_y,
                                       double dc_in
)
        : Point(dim_in, num_in, coor_x, coor_y, 0., jacobian_in, alpha_in, dc_in)
        , weight(weight_in)
{
}

GaussPointBoundary::~GaussPointBoundary()
{
}


void GaussPointBoundary::shapeFunSolve(std::string type_in, double q_in)
{
    cout << "INFO AT shapeFunSolve IN GaussPointBoundary: (x, y) = "
    << this->X << ", " << this->Y << endl;
    cout << "\t alphai = " << alphai << ", "
    << "dc = " << dc << ", "
    << "q_in = " << q_in
    << endl;
    cout << "\t Number of Support Nodes = " << supportNodesSize << endl;

    q_in = .5; // TODO: Take into account and validate the variable q_in
    if (!shapeFun) {
        if (type_in == "RBF") {
            this->shapeFun = new ShapeFunctionRBF(supportNodesSize,
                                                  0,
                                                  0, // RBF type
                                                  alphai,
                                                  dc,
                                                  q_in,
                                                  this);
        } else if (type_in == "MLS") {
            this->shapeFun = new ShapeFunctionMLS(supportNodesSize,
                                                  1,
                                                  1, // weight type
                                                  alphai,
                                                  dc,
                                                  this);
        } else if (type_in == "1D-X") {
            shapeFun = new ShapeFunctionLinearX(this);
        }

        shapeFun->calc();
    }
    initializeMatVecs();
}


void GaussPointBoundary::computeQext(LoadThermalBoundary1D * loadThermalBoundary_in)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        // Qi = wg * r * N_I * |Jc|
        Qext.writeElement(weight * shapeFun->getPhi(0, i) * loadThermalBoundary_in->getLoadThermalBoundary1D(this)
                          * std::abs(jacobian), i);
//       M.writeElement( weight * shapeFun->getPhi(0,i) * shapeFun->getPhi(0,j) * jacobian, dim*i+2, dim*j+2 );
    }
//     cout << weight << ", " 
// 	 << shapeFun->getPhi(0,i)  <<  ", " 
// 	 << loadThermalBoundary_in->getLoadThermalBoundary1D( this ) <<  ", " 
// 	 << std::abs(jacobian) << endl; 
}


void GaussPointBoundary::assembleQext(lmx::Vector<data_type>& globalHeat)
{
    for (auto i = 0u; i < supportNodesSize; ++i) {
        globalHeat.addElement(Qext.readElement(i),
                              supportNodes[i]->getThermalNumber()
        );
    }
//     cout << Qext << endl; 
}

/*
void GaussPointBoundary::gnuplotOutStress( std::ofstream & gptension )
{
    gptension << X << " " << Y << " " << tension(0) << endl;
}*/

void GaussPointBoundary::initializeMatVecs()
{
    Qext.resize(supportNodesSize);
}


} //Namespace mknix
