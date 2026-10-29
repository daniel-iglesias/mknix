#include "reader_test_support.h"

#include <simulation/simulation.h>

#include <cstdio>
#include <fstream>
#include <sstream>

void testReader()
{
    reader_test::TemporaryInput input(
        "TITLE reader_unit_case\n"
        "DIMENSION 3\n"
        "GRAVITY 1.5 -2.0 0.25\n"
        "SYSTEM child\n"
        "ENDSYSTEM\n");

    mknix::Simulation simulation;
    simulation.inputFromFile(input.path());

    reader_test::expect(mknix::Simulation::getDim() == 3, "dimension keyword updates simulation");
    reader_test::expectNear(mknix::Simulation::getGravity(0), 1.5, "gravity x");
    reader_test::expectNear(mknix::Simulation::getGravity(1), -2.0, "gravity y");
    reader_test::expectNear(mknix::Simulation::getGravity(2), 0.25, "gravity z");
    reader_test::expect(simulation.getInterfaceNumberOfNodes("child") == 0,
                        "SYSTEM keyword creates an empty child system");

    std::ifstream log("output.reader");
    std::ostringstream contents;
    contents << log.rdbuf();
    log.close();
    std::remove("output.reader");
    reader_test::expectContains(contents.str(), "TITLE: reader_unit_case", "title directive log");
    reader_test::expectContains(contents.str(), "SYSTEM: baseSystem.child", "system directive log");
}