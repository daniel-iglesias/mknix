//-- Licencia --
#include "point.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctionMLS.h"
#include "shapefunctiontriangle.h"
#include "shapefunctiontetrahedron.h"
#include "shapefunctionlinear.h"
#include "shapefunctionlinear-x.h"

#include <simulation/simulation.h>

namespace mknix
{

Point::Point()
{
}

Point::Point(const Point& point_in)
    : dim(point_in.dim)
    , num(point_in.getNumber())
    , X(point_in.X)
    , Y(point_in.Y)
    , Z(point_in.Z)
    , alphai(point_in.alphai)
    , dc(point_in.dc)
    , shapeFunType(point_in.shapeFunType)
    , shapeFun(point_in.shapeFun)
    , supportNodes(point_in.supportNodes)
    , supportNodesSize(point_in.supportNodesSize)
    , jacobian(point_in.jacobian)
{
}

Point::Point(const Point * point_in)
    : dim(point_in->dim)
    , num(point_in->getNumber())
    , X(point_in->X)
    , Y(point_in->Y)
    , Z(point_in->Z)
    , alphai(point_in->alphai)
    , dc(point_in->dc)
    , shapeFunType(point_in->shapeFunType)
    , shapeFun(point_in->shapeFun)
    , supportNodes(point_in->supportNodes)
    , supportNodesSize(point_in->supportNodesSize)
    , jacobian(point_in->jacobian)
{
}

Point::Point(int i,
             double coor_x,
             double coor_y,
             double coor_z
            )
    : dim(Simulation::getDim())
    , num(i)
    , X(coor_x)
    , Y(coor_y)
    , Z(coor_z)
    , alphai(0)
    , dc(0)
    , shapeFun(0)
    , supportNodesSize(0)
    , jacobian(0) // initialized to one for rigid body additional points
{
}

Point::Point(int dim_in,
             int i,
             double coor_x,
             double coor_y,
             double coor_z,
             double alpha_in,
             double dc_in
            )
    : dim(dim_in)
    , num(i)
    , X(coor_x)
    , Y(coor_y)
    , Z(coor_z)
    , alphai(alpha_in)
    , dc(dc_in)
    , shapeFun(0)
    , supportNodesSize(0)
    , jacobian(0)
{
}

Point::Point(int dim_in,
             int i,
             double coor_x,
             double coor_y,
             double coor_z,
             double jacobian_in,
             double alpha_in,
             double dc_in
            )
    : dim(dim_in)
    , num(i)
    , X(coor_x)
    , Y(coor_y)
    , Z(coor_z)
    , alphai(alpha_in)
    , dc(dc_in)
    , shapeFun(0)
    , supportNodesSize(0)
    , jacobian(jacobian_in)
{
}

Point::~Point()
{
    if (shapeFun)
    {
        delete (shapeFun);
    }
}

double Point::getTemp() const
{
    double conf_value(0);
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        conf_value += shapeFun->getPhi(0, i) * (supportNodes[i]->getqt());
    }
    return conf_value;
}

double Point::getConf(int dof) const
{
    double conf_value(0);
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        conf_value += shapeFun->getPhi(0, i) * (supportNodes[i]->getqx(dof));
    }
    return conf_value;
}

double Point::distance(Point& node_in) const
{
    return std::sqrt(std::pow(node_in.getX() - X, 2)
                     + std::pow(node_in.getY() - Y, 2)
                     + std::pow(node_in.getZ() - Z, 2));
}

int Point::getSupportNodeNumber(int deriv, int s_node)
{
    return supportNodes[s_node]->getNumber();
}


void Point::addSupportNode(Node * node_in)
{
    supportNodes.push_back(node_in);
    ++supportNodesSize;
}


void Point::findSupportNodes(std::vector<Node *>& domainNodes)
{
    double di = alphai * dc;
    double smax = 0.;
    double distance;
    double minX = X - di;
    double maxX = X + di;
    double minY = Y - di;
    double maxY = Y + di;

//  cout << alphai << "*"<<dc<<" = " << di << endl;
//  cout << "("<<minX<<"-"<<maxX<<")"<<",("<<minY<<"-"<<maxY<<")"<< endl;

    if (dim == 2)
    {
        for (auto& node : domainNodes)
        {
            // enclose first in a rectangle, seems more efficient:
            if (node->getX() >= minX && node->getX() <= maxX)
            {
                if (node->getY() >= minY && node->getY() <= maxY)
                {
                    // distance criteria:
                    distance = sqrt(pow(node->getX() - X, 2)
                                    + pow(node->getY() - Y, 2));
                    if (distance <= di)
                    {
                        //          cout << "Distance = " << distance << endl;
                        supportNodes.push_back(node);
                        if (distance > smax) smax = distance;
                    }
                }
            }
        }
    }
    else if (dim == 3)
    {
        double maxZ = Z + di;
        double minZ = Z - di;
        for (auto& node : domainNodes)
        {
            // enclose first in a rectangle, seems more efficient:
            if (node->getX() >= minX && node->getX() <= maxX)
            {
                if (node->getY() >= minY && node->getY() <= maxY)
                {
                    if (node->getZ() >= minZ && node->getZ() <= maxZ)
                    {
                        // distance criteria:
                        distance = sqrt(pow(node->getX() - X, 2)
                                        + pow(node->getY() - Y, 2)
                                        + pow(node->getZ() - Z, 2));
                        if (distance <= di)
                        {
                            //          cout << "Distance = " << distance << endl;
                            supportNodes.push_back(node);
                            if (distance > smax) smax = distance;
                        }
                    }
                }
            }
        }
    }
//  cout << "Old alphai = " << alphai;
// TODO: investigate if it's worthy to set automatically alphai
    alphai = smax / dc;
//  cout << ", New alphai = " << alphai << endl;
    supportNodesSize = supportNodes.size();
}


// for rectangular domains
void Point::findSupportNodes(std::vector<Node *>& domainNodes,
                             double domain_minX,
                             double domain_maxX,
                             double domain_minY,
                             double domain_maxY
                            )
{
    double dix = alphai * dc;
    double diy = alphai * dc;

    for (auto& node : domainNodes)
    {
        // First, check if bounding box is outside the cells' patch domain
        // and define the box of searching:
        double minX;
        double maxX;
        if (X - dix < domain_minX)
        {
            minX = domain_minX;
            maxX = minX + 2 * dix;
        }
        else if (X + dix > domain_maxX)
        {
            maxX = domain_maxX;
            minX = maxX - 2 * dix;
        }
        else
        {
            minX = X - dix;
            maxX = X + dix;
        }
        double minY;
        double maxY;
        if (Y - diy < domain_minY)
        {
            minY = domain_minY;
            maxY = minY + 2 * diy;
        }
        else if (Y + diy > domain_maxY)
        {
            maxY = domain_maxY;
            minY = maxY - 2 * diy;
        }
        else
        {
            minY = Y - dix;
            maxY = Y + dix;
        }

        // Search for nodes inside the box:
        if (node->getX() >= minX && node->getX() <= maxX)
        {
            if (node->getY() >= minY && node->getY() <= maxY)
            {
                supportNodes.push_back(node);
            }
        }
    }
    supportNodesSize = supportNodes.size();
}

void Point::shapeFunSolve(std::string type_in, double q_in)
{
    if (type_in == "RBF" || type_in == "MLS")
    {
        cout << "INFO AT shapeFunSolve IN Point: (x, y) = "
             << this->X << ", " << this->Y << endl;
        cout << "\t alphai = " << alphai << ", "
             << "dc = " << dc << ", "
             << "q_in = " << q_in
             << endl;
        cout << "\t Number of Support Nodes = " << supportNodesSize << endl;
        q_in = .5;
    }
    if (!shapeFun)
    {
        if (type_in == "RBF")
        {
            shapeFun = new ShapeFunctionRBF(supportNodesSize,
                                            0,
                                            0, // RBF type
                                            alphai,
                                            dc,
                                            q_in,
                                            this);
        }
        else if (type_in == "MLS")
        {
            shapeFun = new ShapeFunctionMLS(supportNodesSize,
                                            1,
                                            1, // weight type
                                            alphai,
                                            dc,
                                            this);
        }
        else if (type_in == "1D")
        {
            shapeFun = new ShapeFunctionLinear(this);
        }
        else if (type_in == "1D-X")
        {
            shapeFun = new ShapeFunctionLinearX(this);
        }
        else if (type_in == "2D")
        {
            shapeFun = new ShapeFunctionTriangle(this);
        }
        else if (type_in == "3D")
        {
            shapeFun = new ShapeFunctionTetrahedron(this);
        }

        shapeFun->calc();
    }
}

void Point::gnuplotOut(std::ofstream& gpdata)
{
    std::ofstream this_gp("onegp.dat");
    std::ofstream this_gp_sup_nodes("onegpsupportnodes.dat");

    gpdata << X << " " << Y << " " << Z << endl;
    this_gp << X << " " << Y << " " << Z << endl;

    // For output of support nodes.
    for (auto& node : supportNodes)
    {
        this_gp_sup_nodes << node->getX() << " " << node->getY() << " " << node->getZ() << endl;
    }

    // For output of shape function.
    this->shapeFun->gnuplotOut();

    std::ofstream out("gnuplot.plt");

    out << "#!/usr/bin/gnuplot -persist" << endl;
    out << "#     G N U P L O T" << endl;
    out << "set encoding iso_8859_1" << endl;
    out << "set terminal x11" << endl;
    out << "set key nobox right spacing 1.7" << endl;
    out << "set xlabel \"\\textsf{X}\"" << endl;
    out << "set ylabel \"\\textsf{Y}\"" << endl;
    out << "set zlabel \"\\textsf{Z}\"" << endl;
    out << "set title \"\\textsf{title}\"" << endl;
    out << "set xrange [*:*]" << endl;
    out << "set yrange [*:*]" << endl;
    out << "set zrange [*:*]" << endl;
    out << "set samples 25; set isosamples 20" << endl;
//   out << "set pm3d; set hidden" << endl;
    out << "splot 'onegp.dat' with points ps 1 palette" << endl;
    out << "replot 'onegpsupportnodes.dat' with points ps 2 palette" << endl;
    out << "pause -1" << endl;
    out << "set pm3d" << endl;
    out << "splot 'shapefunction.dat' with lines palette linewidth 4" << endl;
    out << "pause -1" << endl;
    out << "splot 'shapefunction_x.dat' with lines palette linewidth 4" << endl;
    out << "pause -1" << endl;
    out << "splot 'shapefunction_y.dat' with lines palette linewidth 4" << endl;
    out << "pause -1" << endl;

//   system("gnuplot -bg white gnuplot.plt");
}


} //Namespace mknix
