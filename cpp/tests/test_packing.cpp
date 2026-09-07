#include "sp/packing.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace sp;
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
void near(double a, double b, double e, const char *m) {
  check(std::abs(a - b) <= e, m);
}
Mesh cube() {
  Mesh m;
  Vec3 p[] = {{0, 0, 0},  {10, 0, 0},  {10, 10, 0},  {0, 10, 0},
              {0, 0, 10}, {10, 0, 10}, {10, 10, 10}, {0, 10, 10}};
  int f[][3] = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
                {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
  for (auto &x : f)
    m.triangles.push_back({p[x[0]], p[x[1]], p[x[2]]});
  return m;
}
int main() {
  try {
    GeometryOptions go;
    go.occupancyAcceleration = false;
    auto c = buildContext(cube(), 1, 0, 1e-9, go);
    State s;
    addSphere(c, s, {5, 5, 8}, 1);
    addSphere(c, s, {5, 5, 3}, 1);
    near(firstContactDistance(c, s, 0, {0, 0, -1}), 3, 1e-8,
         "swept sphere DDA contact");
    s.centres[0] = {5, 5, 5};
    reindex(c, s, 0);
    near(firstContactDistance(c, s, 0, {0, 0, -1}), 0, 1e-8,
         "contact blocks approach");
    check(firstContactDistance(c, s, 0, {0, 0, 1}) > 3.9,
          "contact allows separation");
    check(!canPlace(c, s, {5, 5, 3}, 0.5), "reject overlap");
    check(!canPlace(c, s, {0.5, 5, 5}, 1), "reject bounds");
    Random rng;
    Options o;
    o.quiet = true;
    o.shakeSweeps = 0;
    relax(c, s, o, true, 1, rng);
    near(s.centres[0].z, 5, 1e-8, "old batch fixed");
    near(s.centres[1].z, 1, 1e-8, "new batch settled");
    auto frame = gravityFrame(c, {1e300, 0, 0});
    near(frame.direction.x, 1, 1e-15, "scaled gravity");
    near(frame.height, 10, 1e-12, "gravity height");
    Random replay;
    replay.loadReplay({0.5, 0.5, 0.5}, {});
    o.shakeSweeps = 0;
    auto r = pack(c, {1}, o, replay);
    check(r.report.acceptedCount == 1, "accepted count");
    near(r.state.centres[0].z, 1, 1e-8, "settled centre");
    near(r.totalVolume, 4 * std::acos(-1.0) / 3, 1e-12, "volume");
    near(r.inertia[0], 0.4 * r.masses[0], 1e-12, "sphere inertia");
    near(norm(r.outputCentres[0]), 0, 1e-12, "CoM frame");
    Random impossible;
    impossible.loadReplay({}, {});
    o.maxRefillPasses = 0;
    auto big = buildContext(cube(), 6, 0, 1e-9, go);
    auto cap = pack(big, {6}, o, impossible);
    check(cap.report.initialFailures == 3 && cap.report.capacityWarning,
          "three empty initial attempts");
    check(cap.report.nextUnplacedRadiusIndex == 1, "one based next index");
    bool exhausted = false;
    try {
      replay.uniform();
    } catch (const std::runtime_error &) {
      exhausted = true;
    }
    check(exhausted, "strict replay exhaustion");
    Options batchOptions;
    batchOptions.quiet = true;
    batchOptions.maxAttempts = 1;
    batchOptions.maxCompressionSweeps = 0;
    batchOptions.maxRefillPasses = 0;
    Random batches;
    std::vector<double> tape;
    for (double v : {0.125, 0.125, 0.125, 0.875, 0.875, 0.875, 0.875})
      for (int a = 0; a < 3; ++a)
        tape.push_back(v);
    batches.loadReplay(tape, {});
    auto batchResult = pack(c, {1, 1, 1}, batchOptions, batches);
    check(batchResult.report.acceptedCount == 2 &&
              batchResult.report.initialFailures == 3,
          "initial empty attempts remain cumulative across success");
    check(batches.uniformCount == 21, "batch has exact recorded trial count");
    Options refillOptions = batchOptions;
    refillOptions.maxRefillPasses = 3;
    Random emptyRefill;
    emptyRefill.loadReplay(std::vector<double>(12, 0.3), {});
    auto refilled = pack(big, {6}, refillOptions, emptyRefill);
    check(emptyRefill.uniformCount == 12 &&
              refilled.report.initialFailures == 3,
          "refill traverses all active faces for each empty sweep");
    State lateral;
    addSphere(c, lateral, {8, 4, 4}, 1);
    Options lateralOptions = o;
    lateralOptions.gravity = {-100, 0, 0};
    relax(c, lateral, lateralOptions, true, 0, rng);
    near(lateral.centres[0].x, 1, 1e-8, "arbitrary gravity settles to support");
    Options manyOptions;
    manyOptions.quiet = true;
    manyOptions.maxAttempts = 30;
    manyOptions.maxCompressionSweeps = 8;
    manyOptions.maxRefillPasses = 1;
    Random manyRandom(122);
    std::vector<double> manyRadii(35, 0.6);
    auto many = pack(c, manyRadii, manyOptions, manyRandom);
    check(many.report.acceptedCount == manyRadii.size(), "small ensemble fits");
    for (std::size_t i = 0; i < many.state.radii.size(); ++i) {
      check(pointInside(c, many.state.centres[i]), "packed centres inside");
      for (auto &t : c.mesh.triangles)
        check(pointTriangleDistance(many.state.centres[i], t) >=
                  0.6 - 2 * c.tolerance,
              "packed sphere wall gap");
      for (std::size_t j = 0; j < i; ++j)
        check(norm(many.state.centres[i] - many.state.centres[j]) >=
                  1.2 - 2 * c.tolerance,
              "packed sphere pair gap");
    }
    bool invalidGravity = false;
    try {
      gravityFrame(c, {0, 0, 0});
    } catch (const std::invalid_argument &) {
      invalidGravity = true;
    }
    check(invalidGravity, "zero gravity rejected");
    Options extreme = batchOptions;
    extreme.density = std::numeric_limits<double>::max();
    Random extremeRandom;
    extremeRandom.loadReplay({0.5, 0.5, 0.5}, {});
    bool overflow = false;
    try {
      pack(c, {1}, extreme, extremeRandom);
    } catch (const std::overflow_error &) {
      overflow = true;
    }
    check(overflow, "mass overflow is diagnosed");
    auto largeMesh = cube();
    for (auto &triangle : largeMesh.triangles) {
      triangle.a = triangle.a * 1e8;
      triangle.b = triangle.b * 1e8;
      triangle.c = triangle.c * 1e8;
    }
    auto largeContext = buildContext(largeMesh, 1, 1e9, 1e-9, go);
    Random separatedReplay;
    separatedReplay.loadReplay({0.1, 0.5, 0.5, 0.9, 0.5, 0.5}, {});
    auto separated = pack(largeContext, {1, 1}, batchOptions, separatedReplay);
    near(separated.inertia[0], 2 * 0.4 * (4.0 / 3) * std::acos(-1.0), 1e-12,
         "axial inertia retains intrinsic sphere inertia at large separation");
    std::cout << "packing tests passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
