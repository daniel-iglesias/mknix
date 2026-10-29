#include "reader_test_support.h"

#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

int main()
{
    const std::vector<std::pair<std::string, void (*)()>> tests = {
        {"Reader dispatch and simulation settings", testReader},
        {"ReaderConstraints settings", testReaderConstraints},
        {"ReaderFlex body configuration", testReaderFlex},
        {"ReaderRigid settings", testReaderRigid},
        {"SectionReader fields and nested sections", testSectionReader}
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