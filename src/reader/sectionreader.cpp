//
// Created by jholloc on 25/02/16.
//

#include "sectionreader.h"

#include <sstream>
#include <algorithm>
#include <iostream>

/**
 * @brief Constructs a SectionReader with the given section name.
 * @param name The keyword name of the section (e.g., "SYSTEM").
 */
mknix::SectionReader::SectionReader(const std::string& name)
    : sectionName(name)
{}

/**
 * @brief Reads and parses the section from the input stream, dispatching tokens to fields or sub-sections.
 * @param input Input file stream to read tokens from.
 * @param log Output log stream for diagnostic messages.
 * @param line_no Current line number, incremented as lines are consumed.
 */
void mknix::SectionReader::read(std::ifstream& input, std::ofstream& log, size_t& line_no)
{
    std::string line;
    while (std::getline(input, line))
    {
        ++line_no;

        std::cerr << line_no << " [" << sectionName << "]: " << line << std::endl;

        if (line.empty())
        {
            continue;
        }

        std::stringstream ss(line);
        std::string name;
        ss >> name;

        if (name == ("END" + sectionName))
        {
            break;
        }

        std::ostringstream rem;
        for (std::string word, delim; ss >> word;)
        {
            rem << delim << word;
            delim = " ";
        }

        auto fieldName = std::find(fieldNames.begin(), fieldNames.end(), name);
        auto subSection = std::find_if(subSections.begin(), subSections.end(), [name](SectionReader& subSec)
        {
            return subSec.sectionName == name;
        });

        if (fieldName != fieldNames.end())
        {
            fields.emplace_back(name, rem.str());
        }
        else if (subSection != subSections.end())
        {
            subSection->read(input, log, line_no);
        }
        else
        {
            throw std::logic_error("unexpected field name " + name + " found at line " + std::to_string(line_no));
        }
    }
}

/**
 * @brief Registers a valid field name for this section.
 * @param fieldName The name of the field to accept.
 * @return Reference to this SectionReader for chaining.
 */
mknix::SectionReader& mknix::SectionReader::addField(const std::string& fieldName)
{
    fieldNames.push_back(fieldName);
    return *this;
}

/**
 * @brief Registers a nested sub-section reader.
 * @param subSection The SectionReader instance for the sub-section.
 * @return Reference to this SectionReader for chaining.
 */
mknix::SectionReader& mknix::SectionReader::addSubSection(mknix::SectionReader subSection)
{
    subSections.push_back(subSection);
    return *this;
}