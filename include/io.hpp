#pragma once
#include "packing.hpp"
#include "radii.hpp"
#include <filesystem>
#include <iosfwd>
namespace sp {
// READ STL Accept binary files and strict ASCII solid/facet syntax while
// retaining finite triangle coordinates and their input order.
Mesh readStl(const std::filesystem::path &);
// Read finite scalar values separated by commas or whitespace. This helper is
// intentionally single-column/flat; structured particle records need a parser
// with explicit columns rather than silently flattening them.
std::vector<double> readNumbers(const std::filesystem::path &);
// Persist the established sphere, summary, grid-point and grid-cell CSV files.
void writeOutputs(const std::filesystem::path &directory,
                  const std::string &prefix, const Context &, Result &);
// Print the academic end-of-run summary while preserving the reported frame.
void printSummary(std::ostream &, const Context &, const Result &);
// Write machine-readable report fields and stage timings with JSON escaping.
void writeReportJson(const std::filesystem::path &, const Result &,
                     double importSeconds, double preprocessSeconds,
                     double packingSeconds, double outputSeconds,
                     const Random &, const RadiiGenerationReport * = nullptr);
} // namespace sp
