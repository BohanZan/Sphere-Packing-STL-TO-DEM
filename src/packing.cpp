#include "packing.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// PACKING IMPLEMENTATION
// Keep insertion order as the sphere ID order, use the sparse grid only as a
// candidate accelerator, and leave exact geometry checks in the common gate.
// This mirrors the MATLAB pipeline while replacing temporary array operations
// with RAII containers and continuous C++ sweeps.
namespace sp {
namespace {
constexpr double eps = std::numeric_limits<double>::epsilon();
bool finite(Vec3 p) {
  return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
// SPHERE/TRIANGLE NEIGHBOURS Enumerate a clipped block around INDEX.
// The halo is derived from the permitted move and sphere diameter, so this is
// a complete candidate query even when a path crosses several grid cells.
template <class F> void neighbours(const Context &c, Cell index, F f,
                                   std::int64_t halo = 1) {
  for (auto x = std::max<std::int64_t>(0, index[0] - halo);
       x <= std::min(c.cellCount[0] - 1, index[0] + halo); ++x)
    for (auto y = std::max<std::int64_t>(0, index[1] - halo);
         y <= std::min(c.cellCount[1] - 1, index[1] + halo); ++y)
      for (auto z = std::max<std::int64_t>(0, index[2] - halo);
           z <= std::min(c.cellCount[2] - 1, index[2] + halo); ++z)
        f(Cell{x, y, z});
}
// FACE STAMPS Deduplicate triangle IDs during one continuous path query
// without allocating a temporary set for every moved sphere.
struct Stamps {
  std::vector<std::uint32_t> marks;
  std::uint32_t generation = 0;
  explicit Stamps(std::size_t n) : marks(n, 0) {}
  void next() {
    if (++generation == 0) {
      std::fill(marks.begin(), marks.end(), 0);
      generation = 1;
    }
  }
  bool fresh(std::size_t id) {
    if (marks[id] == generation)
      return false;
    marks[id] = generation;
    return true;
  }
};
// FIRSTCONTACTDISTANCE Internal swept-path query.
// Traverse every cell crossed by the proposed translation (DDA), solve sphere
// and finite-triangle contacts, and return the earliest valid distance.
double contact(const Context &c, const State &s, std::size_t id, Vec3 u,
               Stamps &seen,
               double limit = std::numeric_limits<double>::infinity(),
               std::int64_t halo = 1) {
  const Vec3 p = s.centres.at(id);
  const double r = s.radii.at(id);
  double distance = limit;
  for (std::size_t a = 0; a < 3; ++a) {
    if (u[a] > 0)
      distance = std::min(distance, (c.upper[a] - r - p[a]) / u[a]);
    else if (u[a] < 0)
      distance = std::min(distance, (c.lower[a] + r - p[a]) / u[a]);
  }
  distance = std::max(0.0, distance);
  if (distance == 0 || !std::isfinite(distance))
    return distance;
  Cell index = cellIndex(c, p), step{};
  Vec3 crossing, increments;
  for (std::size_t a = 0; a < 3; ++a) {
    crossing[a] = increments[a] = std::numeric_limits<double>::infinity();
    if (u[a] == 0)
      continue;
    step[a] = u[a] > 0 ? 1 : -1;
    double boundary =
        c.lower[a] +
        static_cast<double>(index[a] + (u[a] > 0 ? 1 : 0)) * c.cellSize;
    crossing[a] = std::max(0.0, (boundary - p[a]) / u[a]);
    increments[a] = c.cellSize / std::abs(u[a]);
  }
  seen.next();
  while (true) {
    neighbours(c, index, [&](Cell key) {
      for (auto j : lookup(s.sphereCells, key)) {
        if (j == id || distance == 0)
          continue;
        Vec3 delta = s.centres[j] - p;
        double along = dot(delta, u), radii = r + s.radii[j];
        if (dot(delta, delta) <=
            (radii + c.tolerance) * (radii + c.tolerance)) {
          double roundoff = 64 * eps *
                            (std::abs(delta.x * u.x) + std::abs(delta.y * u.y) +
                             std::abs(delta.z * u.z));
          if (along > roundoff)
            distance = 0;
          continue;
        }
        if (along <= 0)
          continue;
        Vec3 projection = u * along, separation = delta - projection;
        double discriminant = radii * radii - dot(separation, separation),
               error = 32 * eps * radii * radii;
        for (std::size_t a = 0; a < 3; ++a) {
          double component =
              16 * eps * (std::abs(delta[a]) + std::abs(projection[a]));
          error +=
              2 * std::abs(separation[a]) * component + component * component;
        }
        if (discriminant >= -error)
          distance = std::min(
              distance,
              std::max(0.0, along - std::sqrt(std::max(0.0, discriminant))));
      }
    }, halo);
    if (distance == 0)
      return 0;
    for (auto j : lookup(c.triangleCells, index)) {
      if (distance == 0 || !seen.fresh(j))
        continue;
      const auto &t = c.mesh.triangles[j];
      Vec3 finish = p + u * distance, lo = min(p, finish), hi = max(p, finish),
           tl = min(min(t.a, t.b), t.c), th = max(max(t.a, t.b), t.c);
      bool overlap = true;
      for (std::size_t a = 0; a < 3; ++a)
        if (tl[a] > hi[a] + r + c.tolerance || th[a] < lo[a] - r - c.tolerance)
          overlap = false;
      if (overlap)
        distance = std::min(distance, triangleContactDistance(
                                          p, r, u, t, distance, c.tolerance));
    }
    double next = std::min({crossing.x, crossing.y, crossing.z});
    if (distance <= next)
      return distance;
    for (std::size_t a = 0; a < 3; ++a)
      if (crossing[a] == next) {
        index[a] += step[a];
        if (index[a] < 0 || index[a] >= c.cellCount[a])
          return distance;
        crossing[a] += increments[a];
      }
  }
}
// VALIDATEOPTIONS Check public packing controls before any mutable state is
// changed. Explicitly bounded runs fail if they do not converge.
void validateOptions(const Options &o) {
  if (o.initialRelaxation != InitialRelaxation::Batch &&
      o.initialRelaxation != InitialRelaxation::Layer &&
      o.initialRelaxation != InitialRelaxation::All)
    throw std::invalid_argument("Invalid initial relaxation scope");
  if (!o.maxAttempts || !std::isfinite(o.density) || o.density <= 0 ||
      !std::isfinite(o.compressionTolerance) || o.compressionTolerance < 0 ||
      !std::isfinite(o.buffer) || o.buffer < 0 || !std::isfinite(o.tolerance) ||
      o.tolerance < 0)
    throw std::invalid_argument("Invalid packing options");
  if (o.coordinateFrame != "world" && o.coordinateFrame != "center_of_mass")
    throw std::invalid_argument(
        "coordinateFrame must be world or center_of_mass");
}
// SETTLEBATCH Apply one gravity phase followed by the requested independent
// shaking/compression phases to the same active range.
void settleBatch(const Context &c, State &s, const Options &o,
                 std::size_t start, Random &random) {
  for (bool gravity : {true, false}) {
    auto status = relax(c, s, o, gravity, start, random);
    if (!status.converged)
      throw std::runtime_error(
          std::string(gravity ? "Gravity compression" : "Shaking/compression") +
          " did not converge before --sweeps limit; no packing output written. "
          "Omit --sweeps to run until energy convergence.");
    if (!o.quiet && start < s.radii.size())
      std::cout << '[' << (gravity ? "Gravity compression" : "Shaking/compression")
                << "] batch=" << s.radii.size() - start
                << "; sweeps=" << status.sweeps
                << (o.maxCompressionSweeps ? "; converged\n" : "; disabled\n");
  }
}
} // namespace
// LOADREPLAY Install strict uniform and normal traces for MATLAB/C++ parity
// diagnostics. The counters are reset so a tape always starts at draw zero.
void Random::loadReplay(std::vector<double> u, std::vector<double> n) {
  for (double x : u)
    if (!std::isfinite(x) || x < 0 || x > 1)
      throw std::invalid_argument("Invalid uniform replay value");
  for (double x : n)
    if (!std::isfinite(x))
      throw std::invalid_argument("Invalid normal replay value");
  uniformReplay = std::move(u);
  normalReplay = std::move(n);
  uniformCount = normalCount = 0;
  replay_ = true;
}
// UNIFORM Draw from the replay tape when present, otherwise from the native
// engine. Count every draw so reports can explain trajectory divergence.
double Random::uniform() {
  if (replay_ || !uniformReplay.empty()) {
    if (uniformCount >= uniformReplay.size())
      throw std::runtime_error("Uniform replay exhausted");
    return uniformReplay[uniformCount++];
  }
  ++uniformCount;
  return std::generate_canonical<double, 53>(engine_);
}
// NORMAL Draw the independent shake variate or its recorded replay value.
double Random::normal() {
  if (replay_ || !normalReplay.empty()) {
    if (normalCount >= normalReplay.size())
      throw std::runtime_error("Normal replay exhausted");
    return normalReplay[normalCount++];
  }
  ++normalCount;
  return normal_(engine_);
}
// GRAVITYFRAME Normalize gravity without rotating the mesh and derive relative
// height along -gravity. The frame origin is the lowest mesh vertex in that
// direction, so energy does not depend on the STL's world-coordinate offset.
GravityFrame gravityFrame(const Context &c, Vec3 g) {
  if (!finite(g))
    throw std::invalid_argument("Gravity must be finite and nonzero");
  double scale = std::max({std::abs(g.x), std::abs(g.y), std::abs(g.z)});
  if (scale == 0)
    throw std::invalid_argument("Gravity must be nonzero");
  g = g / scale;
  g = g / norm(g);
  GravityFrame f;
  f.direction = g;
  f.up = -g;
  if (c.mesh.triangles.empty())
    throw std::invalid_argument("Mesh must be nonempty");
  Vec3 anchor = c.mesh.triangles.front().a;
  double bottom = 0, top = 0;
  f.origin = anchor;
  for (auto &t : c.mesh.triangles)
    for (auto v : {t.a, t.b, t.c}) {
      double h = dot(v - anchor, f.up);
      if (h < bottom) {
        bottom = h;
        f.origin = v;
      }
      top = std::max(top, h);
    }
  f.height = top - bottom;
  return f;
}
// CANPLACE Test one candidate against the AABB, accepted spheres, finite STL
// facets and exact/occupancy inside classification. Every insertion path uses
// this gate, so initial placement and refilling share the same safety rules.
bool canPlace(const Context &c, const State &s, Vec3 p, double r,
              std::size_t ignoreId) {
  if (!finite(p) || !std::isfinite(r) || r <= 0)
    return false;
  for (std::size_t a = 0; a < 3; ++a)
    if (p[a] - r < c.lower[a] || p[a] + r > c.upper[a])
      return false;
  bool valid = true;
  Cell index = cellIndex(c, p);
  neighbours(c, index, [&](Cell key) {
    for (auto j : lookup(s.sphereCells, key)) {
      if (j == ignoreId)
        continue;
      double limit = r + s.radii[j] - c.tolerance;
      Vec3 d = s.centres[j] - p;
      if (limit > 0 && dot(d, d) < limit * limit)
        valid = false;
    }
  });
  if (!valid)
    return false;
  // Geometry stores the complete finite contact band for this centre cell.
  if (r > c.tolerance)
    for (auto j : lookup(c.triangleCells, index))
      if (pointTriangleDistance(p, c.mesh.triangles[j]) < r - c.tolerance) {
        valid = false;
        break;
      }
  return valid && pointInside(c, p);
}
// ADDSPHERE Append one accepted sphere in stable insertion order and register
// its ID in the sparse occupancy map at the current centre cell.
void addSphere(const Context &c, State &s, Vec3 p, double r) {
  if (!finite(p) || !std::isfinite(r) || r <= 0)
    throw std::invalid_argument("Invalid sphere");
  auto key = cellIndex(c, p);
  auto id = s.radii.size();
  s.centres.push_back(p);
  s.radii.push_back(r);
  s.cellIndices.push_back(key);
  s.sphereCells[key].push_back(id);
}
// REINDEX Move one sphere ID between sparse cells after its centre changes.
// Empty buckets are removed so long runs do not accumulate dead map entries.
void reindex(const Context &c, State &s, std::size_t id) {
  auto old = s.cellIndices.at(id), key = cellIndex(c, s.centres.at(id));
  if (old == key)
    return;
  auto it = s.sphereCells.find(old);
  if (it == s.sphereCells.end())
    throw std::logic_error("Missing sphere cell");
  auto &members = it->second;
  auto p = std::find(members.begin(), members.end(), id);
  if (p == members.end())
    throw std::logic_error("Missing sphere index");
  members.erase(p);
  if (members.empty())
    s.sphereCells.erase(it);
  s.sphereCells[key].push_back(id);
  s.cellIndices[id] = key;
}
// FIRSTCONTACTDISTANCE Public checked wrapper for the continuous swept query.
double firstContactDistance(const Context &c, const State &s, std::size_t id,
                            Vec3 direction) {
  if (!finite(direction) || std::abs(norm(direction) - 1) > 1e-10)
    throw std::invalid_argument("Movement direction must have unit length");
  Stamps seen(c.mesh.triangles.size());
  return contact(c, s, id, direction, seen);
}
// RELAX Implement Algorithms 2 and 3: continuous gravity motion and lateral
// shakes. Each successful move is limited by min(contact distance, b_u), while
// the selected active range supplies the same particles to motion, energy and
// reindexing. Older particles remain collision obstacles.
Relaxation relax(const Context &c, State &s, const Options &o, bool globalPass,
                 std::size_t start, Random &random) {
  if (start >= s.radii.size())
    return {};
  validateOptions(o);
  const double maxRadius = *std::max_element(s.radii.begin(), s.radii.end());
  const double buffer = o.buffer > 0 ? o.buffer : maxRadius;
  // Lists remain at sweep-start cells while coordinates change. A previously
  // moved sphere can be buffer farther from its recorded centre. Eq. (8)
  // gives halo=1; support directly constructed finer contexts conservatively.
  const double haloWidth = std::ceil((2 * maxRadius + buffer) / c.cellSize);
  const auto halo = static_cast<std::int64_t>(std::min(
      haloWidth, static_cast<double>(*std::max_element(c.cellCount.begin(),
                                                      c.cellCount.end()))));
  auto frame = gravityFrame(c, o.gravity);
  auto energy = [&]() {
    double e = 0;
    for (auto i = start; i < s.radii.size(); ++i) {
      double relativeRadius = s.radii[i] / maxRadius;
      // The common height*radius^3 factor cancels in Eq. (30). Normalize
      // first to avoid false convergence/never-ending loops from under/overflow.
      e += (dot(s.centres[i] - frame.origin, frame.up) / frame.height) *
           relativeRadius * relativeRadius * relativeRadius;
    }
    if (!std::isfinite(e))
      throw std::overflow_error("Nonfinite relaxation energy");
    return e;
  };
  double previous = energy();
  Stamps seen(c.mesh.triangles.size());
  auto sweep = [&](bool shake) {
    for (auto i = start; i < s.radii.size(); ++i) {
      Vec3 u = frame.direction;
      if (shake) {
        Vec3 v{random.normal(), random.normal(), random.normal()};
        v = v - frame.direction * dot(v, frame.direction);
        if (norm(v) < eps) {
          std::size_t a = 0;
          for (std::size_t k = 1; k < 3; ++k)
            if (std::abs(frame.direction[k]) < std::abs(frame.direction[a]))
              a = k;
          v = {};
          v[a] = 1;
          v = v - frame.direction * dot(v, frame.direction);
        }
        u = v / norm(v);
      }
      double d = contact(c, s, i, u, seen, buffer, halo);
      if (d > 0) {
        s.centres[i] = s.centres[i] + u * d;
      }
    }
    // Algorithms 2/3 update the cell lists once the complete sweep is done.
    for (auto i = start; i < s.radii.size(); ++i)
      reindex(c, s, i);
  };
  for (std::size_t k = 0; k < o.maxCompressionSweeps; ++k) {
    if (!globalPass)
      for (std::size_t j = 0; j < o.shakeSweeps; ++j)
        sweep(true);
    sweep(false);
    double e = energy();
    if (std::abs(e - previous) <=
        o.compressionTolerance *
            std::max(std::abs(previous), std::numeric_limits<double>::min()))
      return {k + 1, true};
    previous = e;
  }
  return {o.maxCompressionSweeps, o.maxCompressionSweeps == 0};
}
namespace {
// PLACEINLAYER Generate ordered radii inside one bL-high centre-height band.
// Empty candidate calls do not alter the radius cursor; the caller closes a
// layer only after its prescribed empty-attempt rule is reached.
void placeInLayer(const Context &c, State &s,
                  const std::vector<double> &radii, std::size_t &next,
                  const Options &o, Random &random,
                  const GravityFrame *frame = nullptr, double layerBase = 0) {
  std::size_t heightAxis = 0;
  if (frame)
    for (std::size_t a = 1; a < 3; ++a)
      if (std::abs(frame->up[a]) > std::abs(frame->up[heightAxis]))
        heightAxis = a;
  while (next < radii.size()) {
    double r = radii[next];
    Vec3 lo = c.lower + Vec3{r, r, r}, hi = c.upper - Vec3{r, r, r};
    if (lo.x > hi.x || lo.y > hi.y || lo.z > hi.z)
      return;
    bool placed = false;
    for (std::size_t attempt = 0; attempt < o.maxAttempts; ++attempt) {
      Vec3 draw{random.uniform(), random.uniform(), random.uniform()}, p;
      for (std::size_t a = 0; a < 3; ++a)
        p[a] = lo[a] + draw[a] * (hi[a] - lo[a]);
      if (frame) {
        // Algorithm 1's bL-high candidate slab, measured opposite gravity.
        // Intersect the dominant-axis line with both the slab and AABB.
        double otherHeight = 0;
        for (std::size_t a = 0; a < 3; ++a)
          if (a != heightAxis)
            otherHeight += (p[a] - frame->origin[a]) * frame->up[a];
        double a = frame->origin[heightAxis] +
                   (layerBase - otherHeight) / frame->up[heightAxis];
        double b = frame->origin[heightAxis] +
                   (std::min(frame->height, layerBase + c.cellSize) - otherHeight) /
                       frame->up[heightAxis];
        double lower = std::max(lo[heightAxis], std::min(a, b));
        double upper = std::min(hi[heightAxis], std::max(a, b));
        if (lower > upper)
          continue;
        p[heightAxis] = lower + draw[heightAxis] * (upper - lower);
      }
      if (canPlace(c, s, p, r)) {
        addSphere(c, s, p, r);
        placed = true;
        break;
      }
    }
    if (!placed)
      return;
    ++next;
  }
}
} // namespace
// INITIALPLACEMENT Run Algorithm 1 through the current generation layer. The
// caller supplies the ordered radii and receives the first unresolved index.
void initialPlacement(const Context &c, State &s,
                      const std::vector<double> &radii, std::size_t &next,
                      const Options &o, Random &random) {
  placeInLayer(c, s, radii, next, o, random);
}
// REFILL Run Algorithm 4 over inward-facing surface facets. A successful
// traversal is settled immediately; only empty traversals count toward Mr.
void refill(const Context &c, State &s, const std::vector<double> &radii,
            std::size_t &next, const Options &o, Random &random) {
  if (!o.maxRefillPasses || next >= radii.size())
    return;
  auto f = gravityFrame(c, o.gravity);
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < c.inwardNormals.size(); ++i)
    if (dot(c.inwardNormals[i], f.direction) > 64 * eps)
      active.push_back(i);
  if (active.empty())
    return;
  std::size_t rejects = 0, round = 0;
  while (rejects < o.maxRefillPasses && next < radii.size()) {
    // N0/N1 span the complete surface traversal for this refill round. A
    // sphere accepted on an early face is therefore available as an obstacle
    // on every later face in the same round.
    auto start = s.radii.size();
    ++round;
    for (auto id : active) {
      const auto &t = c.mesh.triangles[id];
      while (next < radii.size()) {
        double r = radii[next];
        bool placed = false;
        for (std::size_t attempt = 0; attempt < o.maxAttempts; ++attempt) {
          double u = random.uniform(), v = random.uniform();
          // Reflect the unit-square sample into the reference triangle. This
          // keeps the barycentric candidate distribution uniform.
          if (u + v > 1) {
            u = 1 - u;
            v = 1 - v;
          }
          Vec3 p =
              t.a + (t.b - t.a) * u + (t.c - t.a) * v + c.inwardNormals[id] * r;
          if (canPlace(c, s, p, r)) {
            addSphere(c, s, p, r);
            placed = true;
            break;
          }
        }
        // Keep the same unresolved radius for the next face; a failed face does
        // not consume the prescribed sequence.
        if (!placed)
          break;
        ++next;
      }
    }
    if (!o.quiet)
      std::cout << "[Refilling] Round " << round << "; added "
                << s.radii.size() - start << "; total " << s.radii.size() << '/'
                << radii.size() << '\n';
    // An empty traversal consumes Mr; every nonempty traversal is compressed
    // immediately and leaves the empty-traversal counter unchanged.
    if (s.radii.size() == start)
      ++rejects;
    else {
      settleBatch(c, s, o, start, random);
    }
  }
}
// PACK Execute placement, A-range settling, surface refilling and final mass
// properties. Valid partial assemblies are returned with capacity_reached.
Result pack(const Context &c, const std::vector<double> &radii,
            const Options &input, Random &random) {
  Options o = input;
  if (o.buffer == 0 && !radii.empty())
    o.buffer = *std::max_element(radii.begin(), radii.end());
  validateOptions(o);
  auto f = gravityFrame(c, o.gravity);
  // Validate the entire prescribed sequence before inserting any sphere. The
  // grid was built from the same maximum radius, so no later candidate can
  // silently exceed the spatial-index invariant.
  for (double r : radii)
    if (!std::isfinite(r) || r <= 0 || 2 * r > c.cellSize * (1 + 64 * eps))
      throw std::invalid_argument(
          "Radii must be finite, positive and fit the context grid spacing");
  Result result;
  auto &s = result.state;
  s.centres.reserve(radii.size());
  s.radii.reserve(radii.size());
  s.cellIndices.reserve(radii.size());
  std::size_t next = 0, failures = 0;
  double layerBase = 0;
  while (next < radii.size() && layerBase < f.height) {
    // A generation layer owns all successful insertion calls until three empty
    // calls close it. Candidate A keeps that ownership during compression even
    // if a sphere subsequently settles below its original height band.
    auto start = s.radii.size();
    auto lastBatchStart = start;
    std::size_t rejects = 0;
    while (next < radii.size() && rejects <= 2) {
      auto before = s.radii.size();
      placeInLayer(c, s, radii, next, o, random, &f, layerBase);
      if (s.radii.size() == before) {
        ++rejects;
        continue;
      }
      // Keep the generation cursor separate from the selected moving range.
      // The same active range is used for motion, energy and cell updates.
      lastBatchStart = before;
      if (!o.quiet)
        std::cout << "[Initial packing] Batch " << s.radii.size() - before
                  << "; total " << s.radii.size() << '/' << radii.size() << '\n';
      auto activeStart = o.initialRelaxation == InitialRelaxation::All ? 0 :
                         o.initialRelaxation == InitialRelaxation::Layer ? start : before;
      settleBatch(c, s, o, activeStart, random);
    }
    failures += rejects;
    if (s.radii.size() == start) {
      // A complex/disconnected domain can have an empty horizontal slab.
      // Advance through it rather than declaring the whole object full.
      layerBase += c.cellSize;
      continue;
    }
    double top = -std::numeric_limits<double>::infinity();
    for (auto i = lastBatchStart; i < s.radii.size(); ++i)
      top = std::max(top, dot(s.centres[i] - f.origin, f.up) + s.radii[i]);
    if (top > f.height - c.cellSize)
      break;
    // The paper fixes slab height but leaves its anchor implicit. Visit
    // contiguous centre-height bands: advancing to a sphere's surface top
    // could skip a still-unvisited interval (up to one radius in height).
    layerBase += c.cellSize;
  }
  // Any radius not placed in the layered pass is offered to the active surface
  // refill routine before final mass properties are assembled.
  if (next < radii.size())
    refill(c, s, radii, next, o, random);
  auto &report = result.report;
  report.initialRelaxationScope = o.initialRelaxation == InitialRelaxation::All ? "all" :
                                  o.initialRelaxation == InitialRelaxation::Layer ? "layer" : "batch";
  report.requestedCount = radii.size();
  report.acceptedCount = s.radii.size();
  report.unplacedCount = radii.size() - s.radii.size();
  report.nextUnplacedRadiusIndex = next + 1;
  report.initialFailures = failures;
  report.refillPasses = o.maxRefillPasses;
  report.boundingBoxDimensions = c.upper - c.lower;
  report.stlVolume = meshSignedVolume(c.mesh);
  report.coordinateFrame = o.coordinateFrame;
  // Compute volume and mass in insertion order. This makes the accepted ID,
  // radius, output centre and mass vectors remain one-to-one.
  Vec3 weighted{};
  for (std::size_t i = 0; i < s.radii.size(); ++i) {
    double volume = (4.0 / 3) * std::acos(-1.0) * s.radii[i] * s.radii[i] *
                    s.radii[i],
           mass = o.density * volume;
    if (!std::isfinite(volume) || !std::isfinite(mass) || volume <= 0 ||
        mass <= 0)
      throw std::overflow_error(
          "Sphere mass or volume is outside the floating-point range");
    result.totalVolume += volume;
    result.masses.push_back(mass);
    report.totalMass += mass;
    weighted = weighted + s.centres[i] * mass;
  }
  // The selected frame is applied only after world-frame mass properties are
  // known; translating coordinates must not change the physical CoM itself.
  if (!s.radii.empty())
    report.centreOfMass = weighted / report.totalMass;
  if (!std::isfinite(report.totalMass) || !std::isfinite(result.totalVolume) ||
      !finite(report.centreOfMass))
    throw std::overflow_error(
        "Assembly mass properties are outside the floating-point range");
  report.sphereAssemblyVolume = result.totalVolume;
  report.coordinateShift =
      o.coordinateFrame == "world" ? Vec3{} : report.centreOfMass;
  Vec3 centredWeighted{};
  // Apply the parallel-axis theorem to the physical CoM. The residual pass
  // mirrors the MATLAB grouping so cancellation is reproducible.
  for (std::size_t i = 0; i < s.radii.size(); ++i) {
    result.outputCentres.push_back(s.centres[i] - report.coordinateShift);
    centredWeighted = centredWeighted +
                      (s.centres[i] - report.centreOfMass) * result.masses[i];
  }
  report.centreOfMassAfterShift =
      centredWeighted / std::max(report.totalMass, eps);
  // MATLAB recomputes the residual CoM before applying the parallel-axis
  // theorem.
  Vec3 residual = s.radii.empty() ? Vec3{} : centredWeighted / report.totalMass;
  for (std::size_t i = 0; i < s.radii.size(); ++i) {
    Vec3 d = s.centres[i] - report.centreOfMass - residual;
    double mass = result.masses[i],
           intrinsic = 0.4 * mass * s.radii[i] * s.radii[i];
    for (std::size_t a = 0; a < 3; ++a)
      for (std::size_t b = 0; b < 3; ++b)
        // Keep MATLAB's grouping: cancel the axial geometric contribution
        // before scaling it, and add intrinsic inertia independently.
        result.inertia[a * 3 + b] =
            (result.inertia[a * 3 + b] + (a == b ? intrinsic : 0)) +
            mass * ((a == b ? dot(d, d) : 0) - d[a] * d[b]);
  }
  for (double value : result.inertia)
    if (!std::isfinite(value))
      throw std::overflow_error(
          "Assembly inertia is outside the floating-point range");
  if (report.unplacedCount) {
    report.capacityWarning = true;
    report.stopReason = "capacity_reached";
  }
  return result;
}
} // namespace sp
