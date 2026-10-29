#include <simulation/simulation.h>

#include <cmath>
#include <stdexcept>

void testSimulationSettings()
{
    if (mknix::Simulation::getDim() != 2)
    {
        throw std::runtime_error("default simulation dimension should be 2");
    }
    if (mknix::Simulation::getTime() != 0.0)
    {
        throw std::runtime_error("default simulation time should be zero");
    }
    if (mknix::Simulation::getConstraintMethod() != "PENALTY")
    {
        throw std::runtime_error("default constraint method should be PENALTY");
    }
    if (mknix::Simulation::getSmoothingType() != "GLOBAL")
    {
        throw std::runtime_error("default smoothing type should be GLOBAL");
    }
    if (mknix::Simulation::getAlpha() <= 0.0)
    {
        throw std::runtime_error("default constraint alpha should be positive");
    }

    mknix::Simulation::setAugmentedTolerance(3.25);
    if (std::fabs(mknix::Simulation::getAugmentedTolerance() - 3.25) > 1e-12)
    {
        throw std::runtime_error("augmented tolerance setter/getter mismatch");
    }
    mknix::Simulation::setAugmentedTolerance(5.0);
}