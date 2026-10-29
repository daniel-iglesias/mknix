#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void testSystemAndMotion();
void testLoadsAndConstraints();

int main()
{
    const std::vector<std::pair<std::string, void (*)()>> tests = {
        {"System queries and Motion interpolation", testSystemAndMotion},
        {"Force, thermal load, contact, and distance constraint", testLoadsAndConstraints}
    };

    int failures = 0;
    for (std::vector<std::pair<std::string, void (*)()>>::const_iterator test = tests.begin();
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