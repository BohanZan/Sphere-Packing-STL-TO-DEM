#pragma once
#include "geometry.hpp"
#include <limits>
#include <random>
#include <string>
namespace sp {
// Select which already generated spheres are movable during initial packing.
// Layer is the default candidate A: all successful insertions in the current
// generation layer are settled together before the layer is closed.
enum class InitialRelaxation { Batch, Layer, All };
struct Options {
  std::size_t maxAttempts = 1000,
              maxCompressionSweeps = std::numeric_limits<std::size_t>::max(),
              shakeSweeps = 2,
              maxRefillPasses = 3;
  double buffer = 0, tolerance = 1e-9, compressionTolerance = 1e-5, density = 1;
  Vec3 gravity{0, 0, -1};
  std::string coordinateFrame = "center_of_mass";
  std::uint32_t randomSeed = 42;
  GeometryOptions geometry;
  InitialRelaxation initialRelaxation = InitialRelaxation::Layer;
  bool quiet = false;
};
class Random {
  // The replay path is deliberately strict: an exhausted tape is an error,
  // which prevents a partial random trace from being mistaken for a valid run.
  std::mt19937 engine_;
  std::normal_distribution<double> normal_{0, 1};
  bool replay_ = false;

public:
  std::vector<double> uniformReplay, normalReplay;
  std::size_t uniformCount = 0, normalCount = 0;
  explicit Random(std::uint32_t seed = 42) : engine_(seed) {}
  void loadReplay(std::vector<double> uniforms, std::vector<double> normals);
  double uniform();
  double normal();
};
struct State {
  // All vectors use insertion order. A sphere ID is stable while its sparse
  // cell membership changes after movement.
  std::vector<Vec3> centres;
  std::vector<double> radii;
  Grid sphereCells;
  std::vector<Cell> cellIndices;
};
struct GravityFrame {
  // direction points with gravity; up is the opposite direction used for
  // candidate layers and relative-height potential energy.
  Vec3 direction, up, origin;
  double height = 0;
};
struct Report {
  // Report fields describe both complete and capacity-limited packings. A
  // capacity warning retains valid spheres instead of discarding the assembly.
  std::size_t requestedCount = 0, acceptedCount = 0, unplacedCount = 0,
              nextUnplacedRadiusIndex = 1, initialFailures = 0,
              refillPasses = 0;
  std::string stopReason = "completed";
  std::string initialRelaxationScope = "batch";
  bool capacityWarning = false;
  Vec3 boundingBoxDimensions, centreOfMass, coordinateShift,
      centreOfMassAfterShift;
  double stlVolume = 0, sphereAssemblyVolume = 0, totalMass = 0;
  std::string coordinateFrame = "center_of_mass";
  std::vector<std::string> outputFiles;
};
struct Result {
  State state;
  std::vector<Vec3> outputCentres;
  std::vector<double> masses;
  double totalVolume = 0;
  std::array<double, 9> inertia{};
  Report report;
};
// NORMALISE GRAVITY and derive a model-relative height frame without rotating
// the mesh or changing the world coordinates.
GravityFrame gravityFrame(const Context &, Vec3);
// Reject a candidate against the AABB, accepted spheres, finite STL facets and
// the exact inside test. This is the one placement gate used by all phases.
bool canPlace(const Context &, const State &, Vec3, double,
              std::size_t ignoreId = std::numeric_limits<std::size_t>::max());
// Append an accepted sphere and register it in its current sparse cell.
void addSphere(const Context &, State &, Vec3, double);
// Move one sphere ID between sparse cells after its coordinates change.
void reindex(const Context &, State &, std::size_t);
// Query the first contact along a unit translation direction.
double firstContactDistance(const Context &, const State &, std::size_t, Vec3);
struct Relaxation {
  std::size_t sweeps = 0;
  bool converged = true;
};
Relaxation relax(const Context &, State &, const Options &, bool globalPass,
                 std::size_t batchStart, Random &);
// Algorithm 1: sample the ordered radius sequence in a bL-high generation
// layer, accepting each candidate only after all geometric tests pass.
void initialPlacement(const Context &, State &, const std::vector<double> &,
                      std::size_t &nextRadius, const Options &, Random &);
// Algorithm 4: retry unresolved radii on inward-facing STL facets. Only a
// completely empty traversal consumes the cumulative refill-pass allowance.
void refill(const Context &, State &, const std::vector<double> &,
            std::size_t &nextRadius, const Options &, Random &);
// Run the complete packing pipeline and assemble mass properties in the
// original world frame before applying the requested output-coordinate shift.
Result pack(const Context &, const std::vector<double> &, const Options &,
            Random &);
} // namespace sp
