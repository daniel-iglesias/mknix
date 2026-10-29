#include "reader_test_support.h"

#include <reader/sectionreader.h>

#include <sstream>

void testSectionReader()
{
    reader_test::TemporaryInput input(
        "TITLE example system\n"
        "BODY\n"
        "MASS 12 kg\n"
        "ENDBODY\n"
        "ENDSYSTEM\n"
        "AFTER\n");
    std::ofstream log;

    mknix::SectionReader root("SYSTEM");
    root.addField("TITLE");
    root.addSubSection(mknix::SectionReader("BODY").addField("MASS"));

    size_t lineNumber = 0;
    root.read(input.stream(), log, lineNumber);

    reader_test::expect(lineNumber == 5, "line count includes section terminators");
    reader_test::expect(root.getFields().size() == 1, "top-level field count");
    reader_test::expectEqual(root.getFields().at(0).first, "TITLE", "top-level field name");
    reader_test::expectEqual(root.getFields().at(0).second, "example system", "top-level field value");
    reader_test::expect(root.getSubSections().size() == 1, "sub-section count");
    mknix::SectionReader nested = root.getSubSections().at(0);
    reader_test::expectEqual(nested.getFields().at(0).second, "12 kg", "nested field value");

    std::string remaining;
    input.stream() >> remaining;
    reader_test::expectEqual(remaining, "AFTER", "parser stops at matching section terminator");
}
