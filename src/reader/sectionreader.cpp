//
// Created by jholloc on 25/02/16.
//

#include "sectionreader.h"

#include <sstream>
#include <algorithm>
#include <iostream>

mknix::SectionReader::SectionReader(const std::string& name)
        : sectionName(name)
{}

void mknix::SectionReader::read(std::ifstream& input, std::ofstream& log, size_t& line_no)
{
    std::string line;
    while (std::getline(input, line)) {
        ++line_no;

        std::cerr << line_no << " [" << sectionName << "]: " << line << std::endl;

        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string name;
        ss >> name;

        if (name == ("END" + sectionName)) {
            break;
        }

        std::ostringstream rem;
        for (std::string word, delim; ss >> word;) {
            rem << delim << word;
            delim = " ";
        }

        auto fieldName = std::find(fieldNames.begin(), fieldNames.end(), name);
        auto subSection = std::find_if(subSections.begin(), subSections.end(), [name](SectionReader& subSec) {
            return subSec.sectionName == name;
        });

        if (fieldName != fieldNames.end()) {
            fields.emplace_back(name, rem.str());
        } else if (subSection != subSections.end()) {
            subSection->read(input, log, line_no);
        } else {
            throw std::logic_error("unexpected field name " + name + " found at line " + std::to_string(line_no));
        }
    }
}

mknix::SectionReader& mknix::SectionReader::addField(const std::string& fieldName)
{
    fieldNames.push_back(fieldName);
    return *this;
}

mknix::SectionReader& mknix::SectionReader::addSubSection(mknix::SectionReader subSection)
{
    subSections.push_back(subSection);
    return *this;
}