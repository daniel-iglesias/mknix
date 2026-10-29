#ifndef MKNIX_CORE_TEST_SUPPORT_H
#define MKNIX_CORE_TEST_SUPPORT_H

#include <cmath>
#include <stdexcept>
#include <string>

namespace core_test
{

inline void expect(bool condition, const std::string& description)
{
    if (!condition)
    {
        throw std::runtime_error(description);
    }
}

inline void expectNear(double actual, double expected, const std::string& description)
{
    const double tolerance = 1e-9;
    if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance)
    {
        throw std::runtime_error(description + ": expected " +
                                 std::to_string(expected) + ", got " +
                                 std::to_string(actual));
    }
}

}

void testMaterial();
void testPointAndNode();
void testShapeFunctions();
void testCellsAndGaussPoints();

#endif