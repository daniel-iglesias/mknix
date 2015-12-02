//-- Licencia --
#include "cellboundarylinear.h"
#include "gausspointboundary.h"

#include <core/node.h>

namespace mknix {

CellBoundaryLinear::CellBoundaryLinear()
{
}


CellBoundaryLinear::CellBoundaryLinear(std::string formulation_in,
                                       double alpha_in,
                                       int nGPoints_in,
                                       Point * n1_in,
                                       Point * n2_in
)
        : CellBoundary(formulation_in, alpha_in, nGPoints_in
)
{
    this->bodyPoints.push_back(n1_in);
    this->bodyPoints.push_back(n2_in);

    points(0, 0) = n1_in->getX();
    points(0, 1) = n1_in->getY();
    points(0, 2) = n1_in->getZ();
    points(1, 0) = n2_in->getX();
    points(1, 1) = n2_in->getY();
    points(1, 2) = n2_in->getZ();

    this->jacobian = points(0, 0) - points(1, 0);
//     this->jacobian = n1_in->distance(*n2_in); // if in X-only: points(0,0) - points(1,0);

//     cout << this->jacobian << ", " << points(0,0) << ", " << points(1,0) << endl;
    dc = this->jacobian;

    this->createGaussPoints();
}

CellBoundaryLinear::CellBoundaryLinear(std::string formulation_in,
                                       double alpha_in,
                                       int nGPoints_in,
                                       Point * n1_in,
                                       Point * n2_in,
                                       double dc_in
)
        : CellBoundary(
        formulation_in, alpha_in, nGPoints_in
)
{
    this->bodyPoints.push_back(n1_in);
    this->bodyPoints.push_back(n2_in);

    points(0, 0) = n1_in->getX();
    points(0, 1) = 0.;
    points(0, 2) = 0.;
    points(1, 0) = n2_in->getX();
    points(1, 1) = 0.;
    points(1, 2) = 0.;

    this->jacobian = points(0, 0) - points(1, 0);

    dc = dc_in;

    this->createGaussPoints();
}

CellBoundaryLinear::~CellBoundaryLinear()
{
}


void CellBoundaryLinear::initialize(std::vector<Node *>& nodes_in)
{
// Meshfree CellBoundarys, 
    CellBoundary::initialize(nodes_in);
// function is specialized for FEM here
    if (formulation != "RPIM" || formulation != "EFG") {
        // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
        int i;
        for (auto& point : gPoints) {
            for (i = 0; i < 2; ++i) {
                point->addSupportNode(dynamic_cast<Node *>(this->bodyPoints[i]));
            }
        }
    }
}


void CellBoundaryLinear::computeShapeFunctions()
{
// Meshfree CellBoundarys, 
    CellBoundary::computeShapeFunctions();
// function is specialized for FEM here
    if (formulation != "RPIM" || formulation != "EFG") {
        for (auto& point : gPoints) {
            point->shapeFunSolve("1D-X", 1.03); // Could be changed to more generic ´´1D´´
        }
    }
}


void CellBoundaryLinear::createGaussPoints()
{
    lmx::DenseMatrix<double> gCoef(size_type(nGPoints), 2);

    // reference: http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
    if (nGPoints == 1) {
        gCoef(0, 0) = 0.;
        gCoef(0, 1) = 2.;
    }
    else if (nGPoints == 2) {
        double wdist1 = 0.5773502691896257645091487805019574556476017512701268760186023264839776723029333456937153955857495252252087138051355676766566483649996508262705518373647912161760310773007685273559916067003615583077550051041144223011076288835574182229739459904090157105534559538626730166621791266197964892168;
        double w1 = 1.;
        gCoef(0, 0) = -wdist1;
        gCoef(0, 1) = w1;
        gCoef(1, 0) = wdist1;
        gCoef(1, 1) = w1;
    }
    else if (nGPoints == 3) {
        double wdist1 = 0.774596669241483377035853079956479922166584341058318165317514753222696618387395806703857475371734703583260441372189929402637908087832729923135978349224240702213750958202698716256783906245777858513169283405612501838634682531972963691092925710263188052523534528101729260090115562126394576188;
        double wdist2 = 0.;
        double w1 = 8. / 9.;
        double w2 = 5. / 9.;
        gCoef(0, 0) = -wdist1;
        gCoef(0, 1) = w1;
        gCoef(1, 0) = wdist2;
        gCoef(1, 1) = w2;
        gCoef(2, 0) = wdist1;
        gCoef(2, 1) = w1;
    }
    else if (nGPoints == 4) { //must be checked
        double wdist1 = 0.8611363115940525752239464888928095050957253796297176376157219209065294714950488657041623398844793052105769209319781763249637438391157919764084938458618855762872931327441369944290122598469710261906458681564745219362114916066097678053187180580268539141223471780870198639372247416951073770551;
        double wdist2 = 0.3399810435848562648026657591032446872005758697709143525929539768210200304632370344778752804355548115489602395207464932135845003241712491992776363684338328221538611182352836311104158340621521124125023821932864240034767086752629560943410821534146791671405442668508151756169732898924953195536;
        double w1 = 0.3399810435848562648026657591032446872005758697709143525929539768210200304632370344778752804355548115489602395207464932135845003241712491992776363684338328221538611182352836311104158340621521124125023821932864240034767086752629560943410821534146791671405442668508151756169732898924953195536;
        double w2 = 0.6521451548625461426269360507780005927646513041661064595074706804812481325340896482780162322677418404902018960952364978455755577496740182191429757016783303751407135229556360801973666260481564013273531860737119707353160256000107787211587578617532049337456560923057986412084590467808124974086;
        gCoef(0, 0) = -wdist1;
        gCoef(0, 1) = w1;
        gCoef(1, 0) = -wdist2;
        gCoef(1, 1) = w2;
        gCoef(2, 0) = wdist2;
        gCoef(2, 1) = w2;
        gCoef(3, 0) = wdist1;
        gCoef(3, 1) = w1;
    }
//   else if(nGPoints == 20){
//   }

    double x0 = (points(0, 0) + points(1, 0)) / 2.;
    double dx0 = (points(0, 0) - points(1, 0)) / 2.;
    double y0 = (points(0, 1) + points(1, 1)) / 2.;
    double dy0 = (points(0, 1) - points(1, 1)) / 2.;

    for (int i = 0; i < nGPoints; ++i) {
        gPoints.push_back(new GaussPointBoundary(1, //dimension
                                                 this->alpha,
                                                 0.5 * gCoef(i, 1),
                                                 jacobian,
                                                 i,
                                                 x0 + dx0 * gCoef(i, 0),
                                                 y0 + dy0 * gCoef(i, 0),
                                                 dc
        )
        );
    }
}

/*
void CellBoundaryLinear::gnuplotOut( std::ofstream& data, std::ofstream& gpdata )
{
    for (int i=0; i<3; ++i) {
        data << points(i,0) << " " << points(i,1) << " 0\n";
    }
    data << points(0,0) << " " << points(0,1) << " 0\n";
    data << std::endl;

    for(std::vector<GaussPoint*>::iterator it=gPoints.begin();
            it!=gPoints.end();
            ++it) {
        (*it)->gnuplotOut( gpdata );
    }
}*/

} //Namespace mknix
