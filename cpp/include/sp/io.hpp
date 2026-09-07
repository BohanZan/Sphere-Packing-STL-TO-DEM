#pragma once
#include "sp/packing.hpp"
#include <filesystem>
#include <iosfwd>
namespace sp {
Mesh readStl(const std::filesystem::path &);
std::vector<double> readNumbers(const std::filesystem::path &);
void writeOutputs(const std::filesystem::path &directory,
                  const std::string &prefix, const Context &, Result &);
void printSummary(std::ostream &, const Context &, const Result &);
void writeReportJson(const std::filesystem::path &, const Result &,
                     double importSeconds, double preprocessSeconds,
                     double packingSeconds, double outputSeconds,
                     const Random &);
} // namespace sp
