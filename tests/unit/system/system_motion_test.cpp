#include <core/node.h>
#include <system/motion.h>
#include <system/system.h>

#include <cmath>
#include <stdexcept>

void testSystemAndMotion()
{
    mknix::System system("system_unit");
    if (system.getTitle() != "system_unit" || system.getNumberOfNodes() != 0)
    {
        throw std::runtime_error("System title or empty interface-node count is incorrect");
    }
    if (system.getConstraint("missing") != nullptr || !system.flexBodyNames().empty())
    {
        throw std::runtime_error("empty System query returned an unexpected result");
    }

    mknix::Node node(0, 1.0, 2.0, 3.0);
    mknix::Motion motion(&node);
    std::map<double, double> x{{0.0, 0.0}, {10.0, 2.0}};
    std::map<double, double> y{{0.0, 1.0}, {10.0, 3.0}};
    std::map<double, double> z{{0.0, -1.0}, {10.0, 1.0}};
    motion.setTimeUx(x);
    motion.setTimeUy(y);
    motion.setTimeUz(z);
    motion.update(5.0);

    if (std::fabs(node.getUx() - 1.0) > 1e-12 ||
        std::fabs(node.getUy() - 2.0) > 1e-12 ||
        std::fabs(node.getUz()) > 1e-12)
    {
        throw std::runtime_error("Motion did not interpolate displacement at the requested time");
    }
}