#include "reader_test_support.h"

#include <reader/readerflex.h>
#include <system/bodyflexglobalgalerkin.h>
#include <system/system.h>

namespace
{

class InspectableSystem : public mknix::System
{
public:
    using System::System;

    mknix::FlexBody* flexBody(const std::string& name)
    {
        return flexBodies.at(name);
    }

    mknix::Body* thermalBody(const std::string& name)
    {
        return thermalBodies.at(name);
    }
};

}

void testReaderFlex()
{
    reader_test::TemporaryInput input(
        "MESHFREE specimen\n"
        "FORMULATION EFG\n"
        "METHOD EFG\n"
        "INITIALTEMPERATURE 27\n"
        "OUTPUT ENERGY\n"
        "ENDMESHFREE\n"
        "SHARENODES missing peer\n");
    reader_test::TemporaryOutput log;
    InspectableSystem system("flex_test");

    mknix::ReaderFlex reader(nullptr, log.stream(), input.stream());
    reader.readFlexBodies(&system);

    reader_test::expect(system.flexBodyNames().size() == 1,
                        "MESHFREE creates the named flex body");
    reader_test::expectEqual(system.flexBody("specimen")->getType(), "MESHFREE",
                             "flex-body type");
    reader_test::expect(system.thermalBody("specimen") == system.flexBody("specimen"),
                        "flex body is registered as a thermal body");

    const std::string contents = log.contents();
    reader_test::expectContains(contents, "FORMULATION: EFG", "formulation log");
    reader_test::expectContains(contents, "INITIAL TEMPERATURE SET TO: 27",
                                "initial temperature log");
    reader_test::expectContains(contents, "OUTPUT: ENERGY", "output selection log");
    reader_test::expectContains(contents, "SHARENODES: missing peer",
                                "shared-node directive log");
}