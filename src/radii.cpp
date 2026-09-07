#include "radii.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace sp {
namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
bool positive(double x) { return std::isfinite(x) && x > 0; }

// LOGRATIO Evaluate log(a/b) without subtracting two nearly equal logarithms.
// Integer subtraction retains the separation of nearby, large input indices.
double logRatio(std::size_t a, std::size_t b) {
  return a >= b ? std::log1p(double(a-b) / double(b))
                : -std::log1p(double(b-a) / double(a));
}

struct Curve {
  std::size_t count, middle;
  double logMedian, middleGap;
  explicit Curve(const RadiiOptions &o)
      : count(o.count), middle(1+(count-1)/2), logMedian(std::log(o.median)),
        middleGap(count % 2 ? 0 : std::log1p(1.0/double(middle))) {}

  // For even N, divide the shape by the average of its two central values.
  // This correction is bounded by log(2), even for very steep power curves.
  double correction(double p) const {
    return count % 2 ? 0 : std::log(2.0)-std::log1p(std::exp(-p*middleGap));
  }
};

// FITEXPONENT Keep the finite-sequence median fixed and solve its dispersion.
// Radius ratios are scaled to the first (largest) value. expm1 preserves very
// small spreads, and the logarithmic standard deviation avoids squaring large
// radii. No random sampling or array reordering enters the fit.
double fitExponent(const RadiiOptions &o, const Curve &curve, double target) {
  if (target == 0) return 0;
  if (curve.count == 1)
    throw std::invalid_argument("A single radius can only have zero variance");
  if (curve.count == 2 && target >= o.median)
    throw std::invalid_argument(
        "Two radii require standard deviation < median (variance < median^2)");

  std::vector<double> logRatios(curve.count);
  for (std::size_t j = 0; j < curve.count; ++j)
    logRatios[j] = -std::log1p(double(j));
  const double peak = logRatio(curve.middle, 1);
  auto logStd = [&](double p) {
    if (p == 0) return -infinity;
    const double spread = -std::expm1(p*logRatios.back());
    if (spread == 0) return -infinity;
    double mean = 0, m2 = 0;
    std::size_t n = 0;
    for (double ratio : logRatios) {
      const double value = std::expm1(p*ratio)/spread;
      const double delta = value-mean;
      mean += delta/double(++n);
      m2 += delta*(value-mean);
    }
    return curve.logMedian+p*peak+curve.correction(p)+std::log(spread)
           +0.5*std::log(m2/double(curve.count));
  };
  const double wanted = std::log(target);
  double lo = 0, hi = 1;
  while (logStd(hi) < wanted) {
    lo = hi;
    if (hi > std::numeric_limits<double>::max()/2)
      throw std::invalid_argument("Requested dispersion needs an unrepresentable exponent");
    hi *= 2;
  }
  // Bisection preserves the inverse-power family; adding an offset to force
  // the statistics would produce a different law. The residual is checked.
  for (int iteration = 0; iteration < 1100; ++iteration) {
    double p = std::midpoint(lo, hi);
    const double actual = logStd(p);
    if (std::abs(actual-wanted) <= 2e-13) return p;
    if (p == lo || p == hi) break;
    if (actual < wanted) lo = p;
    else hi = p;
  }
  throw std::invalid_argument("Requested dispersion cannot be resolved in double precision");
}
} // namespace

// RADIUSSTATISTICS Use a shifted, scaled Welford accumulation. This avoids
// overflow of the sum and cancellation when radii are large but nearly equal.
// Only a temporary copy is partitioned to obtain the median; input IDs stay put.
RadiiStatistics radiusStatistics(std::span<const double> radii) {
  if (radii.empty() || !std::all_of(radii.begin(), radii.end(), positive))
    throw std::invalid_argument("Radius statistics require a nonempty positive finite sequence");
  RadiiStatistics s;
  s.count = radii.size();
  const auto [lo, hi] = std::minmax_element(radii.begin(), radii.end());
  s.minimum = *lo; s.maximum = *hi;
  const double spread = s.maximum-s.minimum;
  double mean = 0, m2 = 0;
  std::size_t n = 0;
  if (spread != 0) {
    for (double r : radii) {
      const double value = (r-s.minimum)/spread, delta = value-mean;
      mean += delta/double(++n);
      m2 += delta*(value-mean);
    }
  }
  s.mean = std::clamp(s.minimum+spread*mean, s.minimum, s.maximum);
  s.standardDeviation = spread*std::sqrt(std::max(0.0, m2/double(s.count)));
  s.variance = s.standardDeviation*s.standardDeviation;
  std::vector<double> ordered(radii.begin(), radii.end());
  auto mid = ordered.begin()+ordered.size()/2;
  std::nth_element(ordered.begin(), mid, ordered.end());
  s.median = s.count % 2 ? *mid
      : std::midpoint(*std::max_element(ordered.begin(), mid), *mid);
  return s;
}

// GENERATERADII Scale first, then clamp each individual entry. In particular,
// bounds never feed back into the fit and never cause rejection/resampling.
std::vector<double> generateRadii(const RadiiOptions &o,
                                 RadiiGenerationReport *report) {
  // Above 2^53-1 not every integer has a distinct double representation.
  constexpr std::size_t largestExactIndex = 9007199254740991ULL;
  if (!o.count || o.count > largestExactIndex)
    throw std::invalid_argument("Radius count must satisfy 1 <= count <= 2^53-1");
  if (!positive(o.median) || !positive(o.minimum) || !positive(o.maximum)
      || o.minimum > o.maximum)
    throw std::invalid_argument("Radius median and bounds must be positive finite; minimum <= maximum");
  if (int(o.exponent.has_value())+int(o.variance.has_value())
      +int(o.standardDeviation.has_value()) > 1)
    throw std::invalid_argument("Specify only one of radius exponent, variance or standard deviation");
  for (auto value : {o.exponent, o.variance, o.standardDeviation})
    if (value && (!std::isfinite(*value) || *value < 0))
      throw std::invalid_argument("Radius exponent and dispersion must be finite and nonnegative");

  const Curve curve(o);
  const double p = o.variance ? fitExponent(o, curve, std::sqrt(*o.variance))
      : o.standardDeviation ? fitExponent(o, curve, *o.standardDeviation)
                            : o.exponent.value_or(1);
  const double correction = curve.correction(p),
               logMin = std::log(o.minimum), logMax = std::log(o.maximum);
  std::vector<double> radii(curve.count), raw;
  const bool fitted = o.variance.has_value() || o.standardDeviation.has_value();
  const bool needRaw = report || fitted;
  if (needRaw) raw.resize(curve.count);
  std::size_t below = 0, above = 0;
  bool rawAvailable = true;
  for (std::size_t j = 0; j < curve.count; ++j) {
    const double logR = curve.logMedian
        +p*logRatio(curve.middle, j+1)+correction;
    const double r = (p == 0 || (curve.count % 2 && j+1 == curve.middle))
        ? o.median : std::exp(logR);
    if (needRaw) {
      raw[j] = r;
      rawAvailable = rawAvailable && positive(r);
    }
    if (logR < logMin) { radii[j] = o.minimum; ++below; }
    else if (logR > logMax) { radii[j] = o.maximum; ++above; }
    else radii[j] = std::clamp(r, o.minimum, o.maximum);
  }
  std::optional<RadiiStatistics> before;
  if (needRaw && rawAvailable) before = radiusStatistics(raw);
  if (fitted && before) {
    // Solving the mathematical curve is not enough when the spread is smaller
    // than one ULP of a stored radius. Check the actual finite double sequence.
    const double target = o.variance ? std::sqrt(*o.variance) : *o.standardDeviation;
    if (target > 0 && (before->standardDeviation == 0 ||
        std::abs(std::log(before->standardDeviation)-std::log(target)) > 1e-8))
      throw std::invalid_argument("Requested dispersion cannot be represented by the stored double radii");
  }
  if (report) {
    RadiiGenerationReport result;
    result.options = o; result.exponent = p;
    result.clampedBelow = below; result.clampedAbove = above;
    result.beforeClamping = before;
    result.afterClamping = radiusStatistics(radii);
    *report = std::move(result);
  }
  return radii;
}
} // namespace sp
