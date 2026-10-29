#include "reader_test_support.h"

#include <reader/readerconstraints.h>
#include <simulation/simulation.h>
#include <system/system.h>

void testReaderConstraints()
{
    reader_test::TemporaryInput input("PENALTY ALPHA 12.5 AUGMENTED 0.025 ENDJOINTS");
    reader_test::TemporaryOutput log;
    mknix::Simulation simulation;
    mknix::System system("constraints_test");

    mknix::ReaderConstraints reader(&simulation, log.stream(), input.stream());
    reader.readConstraints(&system);

    reader_test::expectEqual(mknix::Simulation::getConstraintMethod(), "AUGMENTED",
                             "constraint method");
    reader_test::expectNear(mknix::Simulation::getAlpha(), 12.5, "constraint alpha");
    reader_test::expectNear(mknix::Simulation::getAugmentedTolerance(), 0.025,
                            "augmented constraint tolerance");
    reader_test::expectContains(log.contents(), "AUGMENTED TOLERANCE: 0.025",
                                "augmented tolerance log");
}