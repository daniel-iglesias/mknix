#include "reader_test_support.h"

#include <reader/readerrigid.h>
#include <simulation/simulation.h>
#include <system/system.h>

void testReaderRigid()
{
    reader_test::TemporaryInput input("PENALTY ALPHA 8.75 AUGMENTED ENDRIGIDBODIES");
    reader_test::TemporaryOutput log;
    mknix::Simulation simulation;
    mknix::System system("rigid_test");

    mknix::ReaderRigid reader(&simulation, log.stream(), input.stream());
    reader.readRigidBodies(&system);

    reader_test::expectEqual(mknix::Simulation::getConstraintMethod(), "AUGMENTED",
                             "rigid-body constraint method");
    reader_test::expectNear(mknix::Simulation::getAlpha(), 8.75, "rigid-body alpha");
    const std::string contents = log.contents();
    reader_test::expectContains(contents, "PENALTY set", "penalty directive log");
    reader_test::expectContains(contents, "AUGMENTED set", "augmented directive log");
}