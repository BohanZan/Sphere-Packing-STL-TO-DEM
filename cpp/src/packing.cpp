#include "sp/packing.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace sp {
namespace {
constexpr double eps = std::numeric_limits<double>::epsilon();
bool finite(Vec3 p) {
  return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
template <class F> void neighbours(const Context &c, Cell index, F f) {
  for (auto x = std::max<std::int64_t>(0, index[0] - 1);
       x <= std::min(c.cellCount[0] - 1, index[0] + 1); ++x)
    for (auto y = std::max<std::int64_t>(0, index[1] - 1);
         y <= std::min(c.cellCount[1] - 1, index[1] + 1); ++y)
      for (auto z = std::max<std::int64_t>(0, index[2] - 1);
           z <= std::min(c.cellCount[2] - 1, index[2] + 1); ++z)
        f(Cell{x, y, z});
}
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
double contact(const Context &c, const State &s, std::size_t id, Vec3 u,
               Stamps &seen) {
  const Vec3 p = s.centres.at(id);
  const double r = s.radii.at(id);
  double distance = std::numeric_limits<double>::infinity();
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
    });
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
void validateOptions(const Options &o) {
  if (!o.maxAttempts || !std::isfinite(o.density) || o.density <= 0 ||
      !std::isfinite(o.compressionTolerance) || o.compressionTolerance < 0 ||
      !std::isfinite(o.buffer) || o.buffer < 0 || !std::isfinite(o.tolerance) ||
      o.tolerance < 0)
    throw std::invalid_argument("Invalid packing options");
  if (o.coordinateFrame != "world" && o.coordinateFrame != "center_of_mass")
    throw std::invalid_argument(
        "coordinateFrame must be world or center_of_mass");
}
} // namespace
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
double Random::uniform() {
  if (replay_ || !uniformReplay.empty()) {
    if (uniformCount >= uniformReplay.size())
      throw std::runtime_error("Uniform replay exhausted");
    return uniformReplay[uniformCount++];
  }
  ++uniformCount;
  return std::generate_canonical<double, 53>(engine_);
}
double Random::normal() {
  if (replay_ || !normalReplay.empty()) {
    if (normalCount >= normalReplay.size())
      throw std::runtime_error("Normal replay exhausted");
    return normalReplay[normalCount++];
  }
  ++normalCount;
  return normal_(engine_);
}
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
double firstContactDistance(const Context &c, const State &s, std::size_t id,
                            Vec3 direction) {
  if (!finite(direction) || std::abs(norm(direction) - 1) > 1e-10)
    throw std::invalid_argument("Movement direction must have unit length");
  Stamps seen(c.mesh.triangles.size());
  return contact(c, s, id, direction, seen);
}
void relax(const Context &c, State &s, const Options &o, bool globalPass,
           std::size_t start, Random &random) {
  if (start >= s.radii.size())
    return;
  auto frame = gravityFrame(c, o.gravity);
  auto energy = [&]() {
    double e = 0;
    for (auto i = start; i < s.radii.size(); ++i)
      e += dot(s.centres[i] - frame.origin, frame.up) * s.radii[i] *
           s.radii[i] * s.radii[i];
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
      double d = contact(c, s, i, u, seen);
      if (d > 0) {
        s.centres[i] = s.centres[i] + u * d;
        reindex(c, s, i);
      }
    }
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
      return;
    previous = e;
  }
}
void initialPlacement(const Context &c, State &s,
                      const std::vector<double> &radii, std::size_t &next,
                      const Options &o, Random &random) {
  while (next < radii.size()) {
    double r = radii[next];
    Vec3 lo = c.lower + Vec3{r, r, r}, hi = c.upper - Vec3{r, r, r};
    if (lo.x > hi.x || lo.y > hi.y || lo.z > hi.z)
      return;
    bool placed = false;
    for (std::size_t attempt = 0; attempt < o.maxAttempts; ++attempt) {
      Vec3 p;
      for (std::size_t a = 0; a < 3; ++a)
        p[a] = lo[a] + random.uniform() * (hi[a] - lo[a]);
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
    auto start = s.radii.size();
    ++round;
    for (auto id : active) {
      const auto &t = c.mesh.triangles[id];
      while (next < radii.size()) {
        double r = radii[next];
        bool placed = false;
        for (std::size_t attempt = 0; attempt < o.maxAttempts; ++attempt) {
          double u = random.uniform(), v = random.uniform();
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
        if (!placed)
          break;
        ++next;
      }
    }
    if (!o.quiet)
      std::cout << "[Refilling] Round " << round << "; added "
                << s.radii.size() - start << "; total " << s.radii.size() << '/'
                << radii.size() << '\n';
    if (s.radii.size() == start)
      ++rejects;
    else {
      relax(c, s, o, true, start, random);
      relax(c, s, o, false, start, random);
    }
  }
}
Result pack(const Context &c, const std::vector<double> &radii,
            const Options &o, Random &random) {
  validateOptions(o);
  auto f = gravityFrame(c, o.gravity);
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
  while (next < radii.size()) {
    auto start = s.radii.size();
    std::size_t rejects = 0;
    while (next < radii.size() && rejects <= 2) {
      auto before = s.radii.size();
      initialPlacement(c, s, radii, next, o, random);
      if (s.radii.size() == before)
        ++rejects;
    }
    failures += rejects;
    if (s.radii.size() == start)
      break;
    if (!o.quiet)
      std::cout << "[Initial packing] Batch " << s.radii.size() - start
                << "; total " << s.radii.size() << '/' << radii.size() << '\n';
    relax(c, s, o, true, start, random);
    relax(c, s, o, false, start, random);
    double top = -std::numeric_limits<double>::infinity();
    for (auto i = start; i < s.radii.size(); ++i)
      top = std::max(top, dot(s.centres[i] - f.origin, f.up) + s.radii[i]);
    if (top > f.height - c.cellSize)
      break;
  }
  if (next < radii.size())
    refill(c, s, radii, next, o, random);
  auto &report = result.report;
  report.requestedCount = radii.size();
  report.acceptedCount = s.radii.size();
  report.unplacedCount = radii.size() - s.radii.size();
  report.nextUnplacedRadiusIndex = next + 1;
  report.initialFailures = failures;
  report.refillPasses = o.maxRefillPasses;
  report.boundingBoxDimensions = c.upper - c.lower;
  report.stlVolume = meshSignedVolume(c.mesh);
  report.coordinateFrame = o.coordinateFrame;
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
