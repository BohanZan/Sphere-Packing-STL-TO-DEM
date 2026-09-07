#include "exception.hpp"
#include "io.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>

// COMMAND-LINE DRIVER
// Parse, validate and time the same high-level stages as the MATLAB entry
// point: import the model, build immutable context, pack, then persist output.
namespace {
double real(std::string_view s) {
  double v;
  auto p = std::from_chars(s.data(), s.data() + s.size(), v);
  if (p.ec != std::errc{} || p.ptr != s.data() + s.size() || !std::isfinite(v))
    throw std::invalid_argument("Invalid number: " + std::string(s));
  return v;
}
std::size_t integer(std::string_view s) {
  std::size_t n;
  auto p = std::from_chars(s.data(), s.data() + s.size(), n);
  if (p.ec != std::errc{} || p.ptr != s.data() + s.size())
    throw std::invalid_argument("Invalid nonnegative integer: " +
                                std::string(s));
  return n;
}
// Keep the command-line contract close to the MATLAB option names while
// making the executable self-describing for ordinary PowerShell use.
constexpr const char *help =
    R"(SpherePacking: standalone C++20 CPU packing, Intel oneAPI build
Usage: sphere_packing --model mesh.stl --radius 0.3 --count 1000 [options]
       sphere_packing --model mesh.stl --radii radii.csv [options]
       sphere_packing --model mesh.stl --generate-radii inverse-power [options]
       sphere_packing --preset buddha [overrides]

  --model PATH             ASCII or binary closed STL surface
  --radius X               Use the same positive radius for every particle (0.3)
  --count N                Number of particles (100; buddha preset: 44000)
  --radii PATH              Ordered positive values, comma/whitespace separated
  --generate-radii inverse-power  r(i)=C*i^(-p) for ALL i=1..count; excludes --radius
  --radius-median X          Median BEFORE clamping (0.3)
  --radius-exponent P        Nonnegative inverse-power exponent (default 1)
  --radius-variance V        Fit exponent to population variance BEFORE clamping
  --radius-stddev S          Alternative to variance; specify only one dispersion
  --radius-min X --radius-max X  Positive clamp bounds (0.1, 0.6)
  --attempts N             Attempts per prescribed radius (default 1000)
  --buffer X               Positive per-move limit and grid padding (auto: max radius)
  --gravity GX GY GZ        Any finite nonzero gravity vector (0 0 -1)
  --sweeps N               Optional sweep guard; failure if not converged (0 disables motion)
  --shake-sweeps N          Lateral shake sweeps (2)
  --initial-relaxation batch|layer|all  Initial moving range (layer, candidate A)
  --refill-passes N         Cumulative empty refill sweep limit (3; 0 disables)
  --tolerance X            Relative geometric tolerance (1e-9)
  --compression-tolerance X  Relative gravitational-energy convergence (1e-5)
  --density X              Positive density (1)
  --frame world|center_of_mass (default center_of_mass)
  --seed N                 Native C++ seed (42); normal stream differs from MATLAB
  --occupancy on|off        Conservative parity acceleration (on)
  --occupancy-cell-size X   Optional width; default radius/2
  --occupancy-max-cells N   Maximum occupancy allocation (2000000)
  --output DIR             CSV directory (default STL directory)
  --prefix NAME            CSV filename stem (default STL stem)
  --no-output              Omit CSV writing for compute-only timing
  --report PATH            Full report + stage times as JSON
  --uniform-tape PATH --normal-tape PATH  Strict MATLAB random-draw replay
  --quiet                  Suppress progress and summary (timing line remains)
  --preset buddha           44000 x 0.3, gravity -Z, attempts60, shake4,
                           refill5, tolerance1e-8, world frame, layer scope; no sweep cap
  --help                   Show this help

Partial packing is a successful result with capacity_reached in the report.
Malformed inputs or failed writes return exit code 1. All numbers use dot decimals.
)";
} // namespace
int main(int argc, char **argv) {
  try {
    // Parse options before importing the STL or allocating packing grids. This
    // keeps malformed requests side-effect free and makes failures diagnosable.
    sp::Options o;
    sp::RadiiOptions radiusOptions;
    std::optional<sp::RadiiGenerationReport> generation;
    std::filesystem::path model, radiiPath, directory, reportPath, uniformTape,
        normalTape;
    std::string prefix;
    double radius = .3;
    std::size_t count = 100;
    bool noOutput = false;
    bool generate = false, radiusParameters = false, explicitRadius = false;
    auto argument = [&](int &i) -> std::string_view {
      if (++i >= argc)
        throw std::invalid_argument("Missing argument value");
      return argv[i];
    };
    for (int i = 1; i < argc; ++i) {
      std::string_view key = argv[i];
      if (key == "--help" || key == "-h") {
        std::cout << help;
        return 0;
      } else if (key == "--preset") {
        if (argument(i) != "buddha")
          throw std::invalid_argument("Unknown preset");
        model = "inputs/greatBudda/greatBudda.stl";
        count = 44000;
        radius = .3;
        // Defaults do not erase the fact that --radius was explicitly supplied.
        // Mutual exclusion with the generator must hold in every argument order.
        o.gravity = {0, 0, -1};
        o.maxAttempts = 60;
        o.compressionTolerance = 1e-8;
        o.shakeSweeps = 4;
        o.maxRefillPasses = 5;
        o.coordinateFrame = "world";
        directory = "results/Great_Budda_cpp";
        prefix = "Great_Budda_packing";
      } else if (key == "--model")
        model = argument(i);
      else if (key == "--radii")
        radiiPath = argument(i);
      else if (key == "--generate-radii") {
        if (argument(i) != "inverse-power")
          throw std::invalid_argument("generate-radii expects inverse-power");
        generate = true;
      } else if (key == "--radius-median") {
        radiusOptions.median = real(argument(i)); radiusParameters = true;
      } else if (key == "--radius-exponent") {
        radiusOptions.exponent = real(argument(i)); radiusParameters = true;
      } else if (key == "--radius-variance") {
        radiusOptions.variance = real(argument(i)); radiusParameters = true;
      } else if (key == "--radius-stddev") {
        radiusOptions.standardDeviation = real(argument(i)); radiusParameters = true;
      } else if (key == "--radius-min") {
        radiusOptions.minimum = real(argument(i)); radiusParameters = true;
      } else if (key == "--radius-max") {
        radiusOptions.maximum = real(argument(i)); radiusParameters = true;
      } else if (key == "--radius") {
        radius = real(argument(i));
        explicitRadius = true;
      } else if (key == "--count") {
        count = integer(argument(i));
      } else if (key == "--attempts")
        o.maxAttempts = integer(argument(i));
      else if (key == "--buffer")
        o.buffer = real(argument(i));
      else if (key == "--gravity") {
        o.gravity.x = real(argument(i));
        o.gravity.y = real(argument(i));
        o.gravity.z = real(argument(i));
      } else if (key == "--sweeps")
        o.maxCompressionSweeps = integer(argument(i));
      else if (key == "--shake-sweeps")
        o.shakeSweeps = integer(argument(i));
      else if (key == "--initial-relaxation") {
        auto scope = argument(i);
        if (scope == "batch") o.initialRelaxation = sp::InitialRelaxation::Batch;
        else if (scope == "layer") o.initialRelaxation = sp::InitialRelaxation::Layer;
        else if (scope == "all") o.initialRelaxation = sp::InitialRelaxation::All;
        else throw std::invalid_argument("initial-relaxation expects batch, layer or all");
      }
      else if (key == "--refill-passes")
        o.maxRefillPasses = integer(argument(i));
      else if (key == "--tolerance")
        o.tolerance = real(argument(i));
      else if (key == "--compression-tolerance")
        o.compressionTolerance = real(argument(i));
      else if (key == "--density")
        o.density = real(argument(i));
      else if (key == "--frame")
        o.coordinateFrame = argument(i);
      else if (key == "--seed") {
        auto n = integer(argument(i));
        if (n > UINT32_MAX)
          throw std::invalid_argument("Seed exceeds uint32");
        o.randomSeed = static_cast<std::uint32_t>(n);
      } else if (key == "--occupancy") {
        auto value = argument(i);
        if (value != "on" && value != "off")
          throw std::invalid_argument("occupancy expects on or off");
        o.geometry.occupancyAcceleration = value == "on";
      } else if (key == "--occupancy-cell-size")
        o.geometry.occupancyCellSize = real(argument(i));
      else if (key == "--occupancy-max-cells")
        o.geometry.occupancyMaxCells = integer(argument(i));
      else if (key == "--output")
        directory = argument(i);
      else if (key == "--prefix")
        prefix = argument(i);
      else if (key == "--report")
        reportPath = argument(i);
      else if (key == "--uniform-tape")
        uniformTape = argument(i);
      else if (key == "--normal-tape")
        normalTape = argument(i);
      else if (key == "--no-output")
        noOutput = true;
      else if (key == "--quiet")
        o.quiet = true;
      else
        throw std::invalid_argument("Unknown option: " + std::string(key));
    }
    // A model is mandatory; the preset only fills the same ordinary options
    // that explicit flags can override.
    if (model.empty())
      throw std::invalid_argument(
          "--model or --preset is required; use --help");
    if (count == 0 || !std::isfinite(radius) || radius <= 0)
      throw std::invalid_argument("count and radius must be positive");
    if (radiusParameters && !generate)
      throw std::invalid_argument("Radius curve parameters require --generate-radii inverse-power");
    if (generate && (!radiiPath.empty() || explicitRadius))
      throw std::invalid_argument("Cannot combine --generate-radii with --radii or --radius");
    // The packing core consumes an ordered radius sequence. Equal-radius runs
    // are expanded here, while --radii preserves the caller's exact order.
    std::vector<double> radii;
    if (generate) {
      // The curve covers the complete input array. Bounds apply to radii only.
      radiusOptions.count = count;
      generation.emplace();
      auto began = std::chrono::steady_clock::now();
      radii = sp::generateRadii(radiusOptions, &*generation);
      generation->generationSeconds = std::chrono::duration<double>(
          std::chrono::steady_clock::now()-began).count();
      if (!o.quiet) {
        const auto &s = generation->afterClamping;
        const auto precision = std::cout.precision();
        std::cout << std::setprecision(12) << "Radius input: count=" << s.count
                  << ", exponent=" << generation->exponent
                  << "\nAfter clamping: min=" << s.minimum << ", max=" << s.maximum
                  << ", median=" << s.median << ", mean=" << s.mean
                  << ", variance=" << s.variance << ", stddev=" << s.standardDeviation
                  << ", clamped below=" << generation->clampedBelow
                  << ", above=" << generation->clampedAbove << '\n';
        std::cout.precision(precision);
      }
    } else radii = radiiPath.empty() ? std::vector<double>(count, radius)
                                    : sp::readNumbers(radiiPath);
    if (radii.empty() || std::any_of(radii.begin(), radii.end(), [](double r) {
          return !std::isfinite(r) || r <= 0;
        }))
      throw std::invalid_argument(
          "radii must be a nonempty ordered positive sequence");
    // Zero means "derive the safe default" rather than "disable movement".
    // Resolve it before context construction because the value sets the grid
    // spacing as well as the per-move displacement limit.
    if (o.buffer == 0)
      o.buffer = *std::max_element(radii.begin(), radii.end());
    if (uniformTape.empty() != normalTape.empty())
      throw std::invalid_argument("Supply both random replay tapes");
    if (prefix.empty())
      prefix = model.stem().string();
    if (directory.empty())
      directory = model.has_parent_path() ? model.parent_path()
                                          : std::filesystem::path(".");
    sp::Random random(o.randomSeed);
    if (!uniformTape.empty())
      random.loadReplay(sp::readNumbers(uniformTape),
                        sp::readNumbers(normalTape));
    // Keep import, preprocessing, packing and output timings separate so a
    // benchmark can compare the algorithm without conflating file I/O.
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto mesh = sp::readStl(model);
    auto imported = Clock::now();
    auto context = sp::buildContext(
        std::move(mesh), *std::max_element(radii.begin(), radii.end()),
        o.buffer, o.tolerance, o.geometry);
    auto preprocessed = Clock::now();
    auto result = sp::pack(context, radii, o, random);
    auto packed = Clock::now();
    if (!noOutput)
      sp::writeOutputs(directory, prefix, context, result);
    auto written = Clock::now();
    auto seconds = [](auto a, auto b) {
      return std::chrono::duration<double>(b - a).count();
    };
    const double ti = seconds(start, imported),
                 tp = seconds(imported, preprocessed),
                 tc = seconds(preprocessed, packed),
                 to = seconds(packed, written);
    if (!reportPath.empty())
      sp::writeReportJson(reportPath, result, ti, tp, tc, to, random,
                          generation ? &*generation : nullptr);
    if (!o.quiet)
      sp::printSummary(std::cout, context, result);
    std::cout << "TIMING import=" << ti << " preprocess=" << tp
              << " packing=" << tc << " output=" << to << " compute=" << tp + tc
              << " total=" << ti + tp + tc + to
              << " accepted=" << result.report.acceptedCount;
    if (generation) std::cout << " radii=" << generation->generationSeconds;
    std::cout << '\n';
    return 0;
  } catch (...) {
    // Convert every expected and unexpected failure into one stable CLI error
    // line; no partial output is claimed as a successful run.
    std::cerr << "SpherePacking error: " << sp::currentExceptionMessage() << '\n';
    return 1;
  }
}
