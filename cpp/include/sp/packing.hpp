#pragma once
#include "sp/geometry.hpp"
#include <limits>
#include <random>
#include <string>
namespace sp {
struct Options {
  std::size_t maxAttempts = 1000, maxCompressionSweeps = 100, shakeSweeps = 2,
              maxRefillPasses = 3;
  double buffer = 0, tolerance = 1e-9, compressionTolerance = 1e-5, density = 1;
  Vec3 gravity{0, 0, -1};
  std::string coordinateFrame = "center_of_mass";
  std::uint32_t randomSeed = 42;
  GeometryOptions geometry;
  bool quiet = false;
};
class Random {
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
  std::vector<Vec3> centres;
  std::vector<double> radii;
  Grid sphereCells;
  std::vector<Cell> cellIndices;
};
struct GravityFrame {
  Vec3 direction, up, origin;
  double height = 0;
};
struct Report {
  std::size_t requestedCount = 0, acceptedCount = 0, unplacedCount = 0,
              nextUnplacedRadiusIndex = 1, initialFailures = 0,
              refillPasses = 0;
  std::string stopReason = "completed";
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
GravityFrame gravityFrame(const Context &, Vec3);
bool canPlace(const Context &, const State &, Vec3, double,
              std::size_t ignoreId = std::numeric_limits<std::size_t>::max());
void addSphere(const Context &, State &, Vec3, double);
void reindex(const Context &, State &, std::size_t);
double firstContactDistance(const Context &, const State &, std::size_t, Vec3);
void relax(const Context &, State &, const Options &, bool globalPass,
           std::size_t batchStart, Random &);
void initialPlacement(const Context &, State &, const std::vector<double> &,
                      std::size_t &nextRadius, const Options &, Random &);
void refill(const Context &, State &, const std::vector<double> &,
            std::size_t &nextRadius, const Options &, Random &);
Result pack(const Context &, const std::vector<double> &, const Options &,
            Random &);
} // namespace sp
