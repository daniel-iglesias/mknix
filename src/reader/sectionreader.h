//
// Created by jholloc on 25/02/16.
//

#ifndef MKNIX_SECTIONREADER_H
#define MKNIX_SECTIONREADER_H

#include <string>
#include <vector>
#include <fstream>

namespace mknix
{

class SectionReader
{
public:
    SectionReader(const std::string& sectionName);
    void read(std::ifstream& input, std::ofstream& log, size_t& line_no);
    const std::vector<std::pair<std::string, std::string>>& getFields()
    {
        return fields;
    };
    const std::vector<SectionReader>& getSubSections()
    {
        return subSections;
    };


    SectionReader& addField(const std::string& fieldName);
    SectionReader& addSubSection(SectionReader subSection);

private:
    const std::string sectionName;
    std::vector<std::string> fieldNames;
    std::vector<SectionReader> subSections;
    std::vector<std::pair<std::string, std::string>> fields;
};

}

#endif //MKNIX_SECTIONREADER_H
