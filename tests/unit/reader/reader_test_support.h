#ifndef MKNIX_READER_TEST_SUPPORT_H
#define MKNIX_READER_TEST_SUPPORT_H

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace reader_test
{

class TemporaryInput
{
public:
    explicit TemporaryInput(const std::string& contents)
    {
        static unsigned long sequence = 0;
        path_ = "mknix_reader_unit_" + std::to_string(++sequence) + ".tmp";
        std::ofstream output(path_.c_str());
        if (!output)
        {
            throw std::runtime_error("unable to create temporary reader input");
        }
        output << contents;
        output.close();
        input_.open(path_.c_str());
        if (!input_)
        {
            throw std::runtime_error("unable to open temporary reader input");
        }
    }

    ~TemporaryInput()
    {
        input_.close();
        std::remove(path_.c_str());
    }

    const std::string& path() const
    {
        return path_;
    }

    std::ifstream& stream()
    {
        return input_;
    }

private:
    std::string path_;
    std::ifstream input_;
};

class TemporaryOutput
{
public:
    TemporaryOutput()
    {
        static unsigned long sequence = 0;
        path_ = "mknix_reader_unit_log_" + std::to_string(++sequence) + ".tmp";
        output_.open(path_.c_str());
        if (!output_)
        {
            throw std::runtime_error("unable to create temporary reader log");
        }
    }

    ~TemporaryOutput()
    {
        output_.close();
        std::remove(path_.c_str());
    }

    std::ofstream& stream()
    {
        return output_;
    }

    std::string contents()
    {
        output_.flush();
        output_.close();
        std::ifstream input(path_.c_str());
        std::ostringstream result;
        result << input.rdbuf();
        return result.str();
    }

private:
    std::string path_;
    std::ofstream output_;
};

inline void expect(bool condition, const std::string& description)
{
    if (!condition)
    {
        throw std::runtime_error(description);
    }
}

inline void expectEqual(const std::string& actual,
                        const std::string& expected,
                        const std::string& description)
{
    if (actual != expected)
    {
        throw std::runtime_error(description + ": expected '" + expected +
                                 "', got '" + actual + "'");
    }
}

inline void expectNear(double actual, double expected, const std::string& description)
{
    const double tolerance = 1e-10;
    if (actual < expected - tolerance || actual > expected + tolerance)
    {
        throw std::runtime_error(description + ": expected " +
                                 std::to_string(expected) + ", got " +
                                 std::to_string(actual));
    }
}

inline void expectContains(const std::string& actual,
                           const std::string& expected,
                           const std::string& description)
{
    if (actual.find(expected) == std::string::npos)
    {
        throw std::runtime_error(description + ": missing '" + expected + "'");
    }
}

}

void testReader();
void testReaderConstraints();
void testReaderFlex();
void testReaderRigid();
void testSectionReader();

#endif