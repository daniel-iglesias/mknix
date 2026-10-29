#include "core_test_support.h"

#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"Material properties and mechanics", testMaterial},
        {"Point and Node", testPointAndNode},
        {"Shape functions", testShapeFunctions},
        {"Cells and Gauss points", testCellsAndGaussPoints}
    };

    int failures = 0;
    for (std::vector<std::pair<std::string, std::function<void()>>>::const_iterator test = tests.begin();
         test != tests.end(); ++test)
    {
        try
        {
            test->second();
            std::cout << "[PASS] " << test->first << std::endl;
        }
        catch (const std::exception& error)
        {
            ++failures;
            std::cerr << "[FAIL] " << test->first << ": " << error.what() << std::endl;
        }
    }

    return failures == 0 ? 0 : 1;
}