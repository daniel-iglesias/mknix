#include <core/node.h>
#include <simulation/simulation.h>
#include <system/constraintdistance.h>
#include <system/force.h>
#include <system/loadthermal.h>
#include <system/loadthermalcontact.h>

#include <cmath>
#include <stdexcept>
#include <string>

void testLoadsAndConstraints()
{
    mknix::Node nodeA(0, 0.0, 0.0, 0.0);
    mknix::Node nodeB(1, 1.0, 0.0, 0.0);

    lmx::Vector<mknix::data_type> forces(4);
    mknix::Force force(&nodeB, 3.0, -4.0, 0.0);
    force.assembleExternalForces(forces);
    if (std::fabs(forces.readElement(2) - 3.0) > 1e-12 ||
        std::fabs(forces.readElement(3) + 4.0) > 1e-12)
    {
        throw std::runtime_error("Force assembly used incorrect global DOF indices or values");
    }

    nodeA.setqt(100.0);
    nodeB.setqt(20.0);
    lmx::Vector<mknix::data_type> heat(2);
    mknix::LoadThermal thermalLoad(&nodeA, 7.5);
    thermalLoad.assembleExternalHeat(heat);
    if (std::fabs(heat.readElement(0) - 7.5) > 1e-12)
    {
        throw std::runtime_error("LoadThermal did not assemble its heat value");
    }

    mknix::LoadThermalContact contact(&nodeA, &nodeB, 0.5);
    contact.assembleExternalHeat(heat);
    if (std::fabs(heat.readElement(0) + 32.5) > 1e-12 ||
        std::fabs(heat.readElement(1) - 40.0) > 1e-12 ||
        std::fabs(heat.readElement(0) + heat.readElement(1) - 7.5) > 1e-12)
    {
        throw std::runtime_error("thermal contact heat transfer is not equal and opposite");
    }

    double alpha = 10.0;
    std::string method = "PENALTY";
    mknix::ConstraintDistance distance(&nodeA, &nodeB, alpha, method);
    distance.calcInternalForces();
    if (std::fabs(distance.getInternalForces().readElement(0)) > 1e-12)
    {
        throw std::runtime_error("distance constraint force should be zero at its reference length");
    }
    nodeB.setUx(1.0);
    distance.calcInternalForces();
    if (std::fabs(distance.getInternalForces().readElement(0) + 120.0) > 1e-12 ||
        std::fabs(distance.getInternalForces().readElement(2) - 120.0) > 1e-12)
    {
        throw std::runtime_error("distance constraint penalty force does not match the stretched residual");
    }
}