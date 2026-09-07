#pragma once
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace sp {
// INDEX-BASED RADIUS INPUT
// Evaluate r(i) = C*i^(-p) for ALL particle indices i=1,...,count.
// Median and dispersion describe the sequence BEFORE clamping. Specify at
// most one of exponent, variance or standardDeviation; the default is p=1.
struct RadiiOptions {
  std::size_t count = 100;
  double median = 0.3, minimum = 0.1, maximum = 0.6;
  std::optional<double> exponent, variance, standardDeviation;
};

struct RadiiStatistics {
  std::size_t count = 0;
  double minimum = 0, maximum = 0, mean = 0, median = 0;
  double variance = 0, standardDeviation = 0;
};

struct RadiiGenerationReport {
  RadiiOptions options;
  double exponent = 1;
  double generationSeconds = 0; // Filled by the command-line driver.
  std::size_t clampedBelow = 0, clampedAbove = 0;
  // Extremely steep curves can exceed double precision before clamping.
  // The bounded output is still valid; unavailable raw statistics are omitted.
  std::optional<RadiiStatistics> beforeClamping;
  RadiiStatistics afterClamping;
};

// GENERATERADII Return radii[j] for particle index i=j+1, without reordering.
// Scale to the requested median, optionally solve p from population variance
// (denominator N), then clamp each radius to [minimum, maximum]. No RNG is used.
std::vector<double> generateRadii(const RadiiOptions &,
                                 RadiiGenerationReport *report = nullptr);

// RADIUSSTATISTICS Summarise a nonempty positive sequence without changing it.
// Even-count medians average the two central values. Variance uses N, not N-1;
// an unrepresentable squared standard deviation is returned as +infinity.
RadiiStatistics radiusStatistics(std::span<const double>);
} // namespace sp
