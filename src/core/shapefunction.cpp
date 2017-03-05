//-- Licencia --
#include "shapefunction.h"
#include "node.h"

#include <simulation/simulation.h>

namespace mknix {

ShapeFunction::ShapeFunction()
        : dim(Simulation::getDim())
{
}

ShapeFunction::ShapeFunction(const ShapeFunction * sf_in)
        : dim(sf_in->dim)
        , phi(sf_in->phi)
        , gp(sf_in->gp)
{
}


ShapeFunction::ShapeFunction(Point * gp_in)
        : dim(gp_in->getDim())
        , gp(gp_in)
{
}


ShapeFunction::~ShapeFunction()
{
}

void ShapeFunction::outputValues()
{
    // output values:
    double tempx = (*gp->supportNodes.begin())->getX();
    auto counter = 0u;


    cout << endl << gp->X << " " << gp->Y << endl;

    for (std::vector<Node *>::iterator it = gp->supportNodes.begin();
         it != gp->supportNodes.end();
         ++it) {
        if (tempx != (*it)->getX()) {
            tempx = (*it)->getX();
        }
        cout << (*it)->getX() << " "
        << (*it)->getY() << " "
        << phi(0, counter) << " "
        << phi(1, counter) << " "
        << phi(2, counter) << endl;

        counter++;
    }
}


void ShapeFunction::gnuplotOut()
{   // for some reason, can't use a std::vector...
    std::ofstream data("shapefunction.dat");
    std::ofstream data_x("shapefunction_x.dat");
    std::ofstream data_y("shapefunction_y.dat");
    double tempx = (*gp->supportNodes.begin())->getX();
    auto counter = 0u;

    for (std::vector<Node *>::iterator it = gp->supportNodes.begin();
         it != gp->supportNodes.end();
         ++it) {
        if (tempx != (*it)->getX()) {
            data << endl;
            data_x << endl;
            data_y << endl;
            tempx = (*it)->getX();
        }

        data  << (*it)->getX() << " "
              << (*it)->getY() << " "
              << phi(0,counter) << endl;
        data_x  << (*it)->getX() << " "
                << (*it)->getY() << " "
                << phi(1,counter) << endl;
        data_y  << (*it)->getX() << " "
                << (*it)->getY() << " "
                << phi(2,counter) << endl;

        counter++;
    }
}

} //Namespace mknix
