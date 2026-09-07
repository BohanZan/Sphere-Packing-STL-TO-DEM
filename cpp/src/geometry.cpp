#include "sp/geometry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace sp {
Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator-(Vec3 a) { return {-a.x, -a.y, -a.z}; }
Vec3 operator*(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
Vec3 operator*(double s, Vec3 a) { return a * s; }
Vec3 operator/(Vec3 a, double s) { return {a.x / s, a.y / s, a.z / s}; }
double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 a) { return std::hypot(a.x, a.y, a.z); }
Vec3 min(Vec3 a, Vec3 b) {
  return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}
Vec3 max(Vec3 a, Vec3 b) {
  return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}
std::size_t CellHash::operator()(const Cell &c) const noexcept {
  std::uint64_t h = 0x9e3779b97f4a7c15ULL;
  for (auto a : c) {
    std::uint64_t x = static_cast<std::uint64_t>(a);
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    h ^= (x ^ (x >> 31)) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  }
  return static_cast<std::size_t>(h);
}
const std::vector<std::size_t> &lookup(const Grid &g, const Cell &c) {
  static const std::vector<std::size_t> empty;
  auto i = g.find(c);
  return i == g.end() ? empty : i->second;
}
namespace {
constexpr double eps = std::numeric_limits<double>::epsilon(),
                 inf = std::numeric_limits<double>::infinity();
Vec3 abs3(Vec3 a) { return {std::abs(a.x), std::abs(a.y), std::abs(a.z)}; }
double largest(Vec3 a) { return std::max({a.x, a.y, a.z}); }
bool finite(Vec3 a) {
  return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
double spacing(double x) {
  x = std::abs(x);
  return std::nextafter(x, inf) - x;
}
std::int64_t index1(double x, std::int64_t n) {
  if (std::isnan(x))
    throw std::invalid_argument("nonfinite grid point");
  if (x <= 0)
    return 0;
  if (x >= static_cast<double>(n - 1))
    return n - 1;
  return static_cast<std::int64_t>(x);
}
Cell indices(Vec3 p, Vec3 lower, double h, Cell counts) {
  Cell c;
  for (int k = 0; k < 3; k++)
    c[k] = index1(std::floor((p[k] - lower[k]) / h), counts[k]);
  return c;
}
Cell countsFor(Vec3 extent, double h) {
  Cell n;
  for (int k = 0; k < 3; k++) {
    double d = std::max(1., std::ceil(extent[k] / h));
    if (!std::isfinite(d) ||
        d >= static_cast<double>(std::numeric_limits<std::int64_t>::max()))
      throw std::overflow_error("grid cell count overflow");
    n[k] = static_cast<std::int64_t>(d);
  }
  return n;
}
std::size_t offset(Cell c, Cell n) {
  return static_cast<std::size_t>(c[0]) +
         static_cast<std::size_t>(n[0]) *
             (static_cast<std::size_t>(c[1]) +
              static_cast<std::size_t>(n[1]) * static_cast<std::size_t>(c[2]));
}
} // namespace
Cell cellIndex(const Context &c, Vec3 p) {
  if (!finite(p))
    throw std::invalid_argument("nonfinite point");
  return indices(p, c.lower, c.cellSize, c.cellCount);
}
Vec3 closestPoint(Vec3 p, const Triangle &t) {
  if (!finite(p) || !finite(t.a) || !finite(t.b) || !finite(t.c))
    throw std::invalid_argument("nonfinite triangle distance input");
  Vec3 a = t.a, b = t.b, c = t.c, ab = b - a, ac = c - a, ap = p - a;
  double scale =
      std::max({largest(abs3(ab)), largest(abs3(ac)), largest(abs3(ap))});
  if (!std::isfinite(scale))
    throw std::overflow_error("triangle distance arithmetic overflow");
  if (scale > 0 && (scale < 1e-50 || scale > 1e50))
    return a + scale * closestPoint(ap / scale,
                                    Triangle{{}, ab / scale, ac / scale});
  if (norm(cross(ab, ac)) == 0) {
    Vec3 q = a;
    double best = inf;
    Vec3 v[] = {a, b, c};
    for (int i = 0; i < 3; i++) {
      Vec3 e = v[(i + 1) % 3] - v[i], x = v[i];
      double den = dot(e, e);
      if (den > 0)
        x = x + std::clamp(dot(p - x, e) / den, 0., 1.) * e;
      double d = dot(p - x, p - x);
      if (d < best) {
        q = x;
        best = d;
      }
    }
    return q;
  }
  double d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0 && d2 <= 0)
    return a;
  Vec3 bp = p - b;
  double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3)
    return b;
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0)
    return a + d1 / (d1 - d3) * ab;
  Vec3 cp = p - c;
  double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6)
    return c;
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0)
    return a + d2 / (d2 - d6) * ac;
  double va = d3 * d6 - d5 * d4;
  if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
    return b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b);
  return a + vb / (va + vb + vc) * ab + vc / (va + vb + vc) * ac;
}
double pointTriangleDistance(Vec3 p, const Triangle &t) {
  return norm(p - closestPoint(p, t));
}
namespace {
double entryRoot(Vec3 off, Vec3 vel, double r, double lo, double hi,
                 Vec3 oe = {}, Vec3 ve = {}) {
  oe = oe + std::abs(lo) * ve + 16 * eps * (abs3(off) + abs3(lo * vel));
  off = off + lo * vel;
  double A = dot(vel, vel), B = dot(off, vel), C = dot(off, off) - r * r;
  if (C <= 0)
    return lo;
  if (A == 0 || B >= 0)
    return inf;
  Vec3 proj = (B / A) * vel, perp = off - proj;
  double ps = dot(perp, perp), disc = r * r - ps;
  Vec3 ce = 16 * eps * (abs3(off) + abs3(proj));
  double extra = 4 * norm(oe) + 4 * norm(off) * norm(ve) / std::sqrt(A);
  ce = ce + Vec3{extra, extra, extra};
  double eb =
      dot(2 * abs3(perp), ce) + dot(ce, ce) + 32 * spacing(std::max(r * r, ps));
  if (disc < -eb)
    return inf;
  double root = C / (-B + std::sqrt(A) * std::sqrt(std::max(0., disc)));
  return root <= hi - lo ? lo + root : inf;
}
double edgeContact(Vec3 off, Vec3 u, Vec3 edge, double r, double lo,
                   double hi) {
  if (lo > hi)
    return inf;
  double l2 = dot(edge, edge);
  if (l2 == 0)
    return entryRoot(off, u, r, lo, hi);
  double s0 = dot(off, edge) / l2, s1 = dot(u, edge) / l2;
  std::vector<double> cuts{lo, hi};
  if (s1 != 0) {
    for (double e : {-s0 / s1, (1 - s0) / s1})
      if (e > lo && e < hi)
        cuts.push_back(e);
  }
  std::sort(cuts.begin(), cuts.end());
  for (size_t i = 0; i + 1 < cuts.size(); i++) {
    double start = cuts[i], finish = cuts[i + 1],
           s = s0 + s1 * (start + (finish - start) / 2), hit;
    if (s < 0)
      hit = entryRoot(off, u, r, start, finish);
    else if (s > 1)
      hit = entryRoot(off - edge, u, r, start, finish);
    else {
      Vec3 op = s0 * edge, vp = s1 * edge;
      hit = entryRoot(off - op, u - vp, r, start, finish,
                      16 * eps * (abs3(off) + abs3(op)),
                      16 * eps * (abs3(u) + abs3(vp)));
    }
    if (std::isfinite(hit))
      return hit;
  }
  return inf;
}
} // namespace
double triangleContactDistance(Vec3 p, double r, Vec3 u, const Triangle &source,
                               double limit, double tolerance) {
  if (!finite(p) || !finite(u) || !finite(source.a) || !finite(source.b) ||
      !finite(source.c) || !std::isfinite(r) || r < 0 || std::isnan(limit) ||
      !std::isfinite(tolerance) || tolerance < 0)
    throw std::invalid_argument("invalid triangle contact input");
  if (limit < 0)
    return inf;
  p = p - source.a;
  Triangle t{{}, source.b - source.a, source.c - source.a};
  double scale =
      std::max({r, largest(abs3(p)), largest(abs3(t.b)), largest(abs3(t.c))});
  if (!std::isfinite(scale))
    throw std::overflow_error("triangle contact arithmetic overflow");
  if (scale > 0 && (scale < 1e-50 || scale > 1e50))
    return scale *
           triangleContactDistance(p / scale, r / scale, u,
                                   Triangle{{}, t.b / scale, t.c / scale},
                                   limit / scale, tolerance / scale);
  Vec3 ab = t.b, ac = t.c, normal = cross(ab, ac);
  double nl = norm(normal), lo = 0, hi = limit, distance = inf;
  if (nl > 0) {
    Vec3 n = normal / nl;
    double h = dot(p, n), vel = dot(u, n);
    if (vel == 0) {
      if (std::abs(h) > r + tolerance)
        return inf;
    } else {
      double t0 = (-r - h) / vel, t1 = (r - h) / vel;
      lo = std::max(0., std::min(t0, t1));
      hi = std::min(limit, std::max(t0, t1));
      if (lo > hi)
        return inf;
    }
    if (std::abs(h) <= r + tolerance) {
      Vec3 delta = p - closestPoint(p, t);
      if (norm(delta) <= r + tolerance)
        return dot(delta, u) < -64 * eps * dot(abs3(delta), abs3(u)) ? 0 : inf;
    }
    if (lo > 0) {
      Vec3 at = p + lo * u;
      double alpha = dot(cross(at, ac), n) / nl,
             beta = dot(cross(ab, at), n) / nl;
      if (alpha >= 0 && beta >= 0 && alpha + beta <= 1)
        return lo;
    }
  } else {
    Vec3 delta = p - closestPoint(p, t);
    if (norm(delta) <= r + tolerance)
      return dot(delta, u) < -64 * eps * dot(abs3(delta), abs3(u)) ? 0 : inf;
  }
  Vec3 v[] = {t.a, t.b, t.c};
  for (int i = 0; i < 3; i++)
    distance =
        std::min(distance, edgeContact(p - v[i], u, v[(i + 1) % 3] - v[i], r,
                                       lo, std::min(hi, distance)));
  return distance;
}
namespace {
// Distance to a conservative projected finite-feature superset. Unlike an
// unguarded barycentric distance this retains numerically ambiguous skinny
// faces.
double finiteDistance(Vec3 point, const Triangle &t) {
  Vec3 rel[] = {Vec3{}, t.b - t.a, t.c - t.a}, q = point - t.a;
  double scale = std::max(largest(abs3(rel[1])), largest(abs3(rel[2])));
  Vec3 n{};
  if (scale > 0)
    n = cross(rel[1] / scale, rel[2] / scale);
  double nl = norm(n), slab = 0;
  if (nl > 0) {
    n = n / nl;
    double su[] = {0, dot(rel[1], n), dot(rel[2], n)}, height = dot(q, n);
    slab = std::max(0., std::max(std::min({su[0], su[1], su[2]}) - height,
                                 height - std::max({su[0], su[1], su[2]})));
    for (int i = 0; i < 3; i++)
      rel[i] = rel[i] - su[i] * n;
    q = q - height * n;
  }
  double projected = inf;
  for (int i = 0; i < 3; i++) {
    Vec3 e = rel[(i + 1) % 3] - rel[i], v = q - rel[i];
    double len = norm(e);
    if (len > 0) {
      Vec3 unit = e / len;
      v = v - std::clamp(dot(v, unit), 0., len) * unit;
    }
    projected = std::min(projected, norm(v));
  }
  if (nl > 0) {
    int axis = 0;
    if (std::abs(n.y) > std::abs(n[axis]))
      axis = 1;
    if (std::abs(n.z) > std::abs(n[axis]))
      axis = 2;
    int other[2], k = 0;
    for (int j = 0; j < 3; j++)
      if (j != axis)
        other[k++] = j;
    Vec3 e = rel[1] / scale, f = rel[2] / scale, v = q / scale;
    int x = other[0], y = other[1];
    double det = e[x] * f[y] - e[y] * f[x];
    if (det != 0) {
      double a = (v[x] * f[y] - v[y] * f[x]) / det,
             b = (e[x] * v[y] - e[y] * v[x]) / det,
             s = std::abs(e[x]) + std::abs(e[y]) + std::abs(f[x]) +
                 std::abs(f[y]);
      double error = 128 * eps * std::max(1., s * s / std::abs(det));
      if (a >= -error && b >= -error && a + b <= 1 + error)
        projected = 0;
    }
  }
  return std::hypot(projected, slab);
}
template <class F> void boxCells(Cell lo, Cell hi, F &&fn) {
  for (auto z = lo[2]; z <= hi[2]; z++)
    for (auto y = lo[1]; y <= hi[1]; y++)
      for (auto x = lo[0]; x <= hi[0]; x++)
        fn(Cell{x, y, z});
}
void raster(const Triangle &original, Vec3 lower, double h, Cell counts,
            double radius, double halo, Grid &grid, size_t id) {
  Triangle tri{original.a - lower, original.b - lower, original.c - lower};
  double magnitude =
      std::max({largest(abs3(original.a)), largest(abs3(original.b)),
                largest(abs3(original.c)), largest(abs3(lower)), h});
  double R = radius + 128 * spacing(magnitude) + halo;
  Vec3 mn = min(min(tri.a, tri.b), tri.c), mx = max(max(tri.a, tri.b), tri.c);
  Cell lo, hi;
  for (int k = 0; k < 3; k++) {
    lo[k] = index1(std::ceil((mn[k] - R) / h + .5) - 2, counts[k]);
    hi[k] = index1(std::floor((mx[k] + R) / h + .5), counts[k]);
  }
  auto put = [&](Cell c) {
    Vec3 p{(c[0] + .5) * h, (c[1] + .5) * h, (c[2] + .5) * h};
    if (finiteDistance(p, tri) <= R)
      grid[c].push_back(id);
  };
  Vec3 v[] = {tri.a, tri.b, tri.c},
       edges[] = {tri.b - tri.a, tri.c - tri.b, tri.a - tri.c};
  int longest = 0;
  for (int i = 1; i < 3; i++)
    if (norm(edges[i]) > norm(edges[longest]))
      longest = i;
  Vec3 P = v[longest], E = edges[longest];
  double length = norm(E);
  long double total = 1;
  for (int k = 0; k < 3; k++)
    total *= static_cast<long double>(hi[k] - lo[k] + 1);
  if (length == 0 || total <= 4096) {
    boxCells(lo, hi, put);
    return;
  }
  Vec3 third = v[(longest + 2) % 3];
  double t = std::clamp(dot(third - P, E / length) / length, 0., 1.),
         thickness = norm(third - P - t * E);
  int axis = 0;
  if (thickness <= h) {
    double S = R + thickness;
    for (int k = 1; k < 3; k++)
      if (std::abs(E[k]) > std::abs(E[axis]))
        axis = k;
    for (auto slice = lo[axis]; slice <= hi[axis]; slice++) {
      double centre = (slice + .5) * h, t0 = (centre - S - P[axis]) / E[axis],
             t1 = (centre + S - P[axis]) / E[axis];
      if (t0 > t1)
        std::swap(t0, t1);
      t0 = std::max(0., t0);
      t1 = std::min(1., t1);
      if (t0 > t1)
        continue;
      Vec3 a = P + t0 * E, b = P + t1 * E;
      Cell l, u;
      for (int k = 0; k < 3; k++) {
        l[k] = std::max(
            lo[k], index1(std::ceil((std::min(a[k], b[k]) - S) / h + .5) - 2,
                          counts[k]));
        u[k] = std::min(
            hi[k],
            index1(std::floor((std::max(a[k], b[k]) + S) / h + .5), counts[k]));
      }
      l[axis] = u[axis] = slice;
      boxCells(l, u, put);
    }
    return;
  }
  Vec3 ab = tri.b - tri.a, ac = tri.c - tri.a;
  double scale = std::max(largest(abs3(ab)), largest(abs3(ac)));
  Vec3 m = cross(ab / scale, ac / scale);
  for (int k = 1; k < 3; k++)
    if (std::abs(m[k]) > std::abs(m[axis]))
      axis = k;
  m = m / largest(abs3(m));
  double s0 = dot(tri.b - tri.a, m), s1 = dot(tri.c - tri.a, m),
         sl = std::min({0., s0, s1}) - R * norm(m),
         su = std::max({0., s0, s1}) + R * norm(m);
  int other[2], k = 0;
  for (int j = 0; j < 3; j++)
    if (j != axis)
      other[k++] = j;
  for (auto j = lo[other[1]]; j <= hi[other[1]]; j++)
    for (auto i = lo[other[0]]; i <= hi[other[0]]; i++) {
      double base = ((i + .5) * h - tri.a[other[0]]) * m[other[0]] +
                    ((j + .5) * h - tri.a[other[1]]) * m[other[1]];
      double a = (sl - base) / m[axis] + tri.a[axis],
             b = (su - base) / m[axis] + tri.a[axis];
      Cell l = lo, u = hi;
      l[other[0]] = u[other[0]] = i;
      l[other[1]] = u[other[1]] = j;
      l[axis] =
          std::max(lo[axis], index1(std::ceil(std::min(a, b) / h + .5) - 2,
                                    counts[axis]));
      u[axis] = std::min(
          hi[axis], index1(std::floor(std::max(a, b) / h + .5), counts[axis]));
      boxCells(l, u, put);
    }
}
std::pair<Vec3, Vec3> rayBounds(const Triangle &t, double tolerance,
                                double magnitude) {
  Vec3 ab = t.b - t.a, ac = t.c - t.a;
  double d = ab.x * ac.y - ab.y * ac.x,
         s = std::abs(ab.x) + std::abs(ab.y) + std::abs(ac.x) + std::abs(ac.y);
  double roundoff = 64 * spacing(magnitude) * std::max(1., s * s / std::abs(d));
  Vec3 pad =
      2 * tolerance * (abs3(ab) + abs3(ac)) + Vec3{roundoff, roundoff, 0};
  return {min(min(t.a, t.b), t.c) - pad, max(max(t.a, t.b), t.c) + pad};
}
} // namespace
bool exactPointInside(const Context &c, Vec3 p) {
  if (!finite(p))
    throw std::invalid_argument("nonfinite point");
  Cell cell = indices(p, c.lower, c.rayCellSize, c.rayCellCount);
  cell[2] = 0;
  std::vector<double> heights;

  for (auto id : lookup(c.rayCells, cell)) {
    auto &t = c.mesh.triangles[id];
    const auto &ray = c.rayTriangles[id];
    if (p.x < ray.lower.x || p.y < ray.lower.y || p.x > ray.upper.x ||
        p.y > ray.upper.y)
      continue;
    Vec3 ab = ray.ab, ac = ray.ac, ap = p - t.a;
    double alpha = (ap.x * ac.y - ap.y * ac.x) * ray.inverseDeterminant,
           beta = (ab.x * ap.y - ab.y * ap.x) * ray.inverseDeterminant;
    if (alpha >= -c.tolerance && beta >= -c.tolerance &&
        alpha + beta <= 1 + c.tolerance) {
      double h = t.a.z + alpha * ab.z + beta * ac.z;
      if (h < p.z - c.tolerance)
        heights.push_back(h);
    }
  }
  if (heights.empty())
    return false;
  std::sort(heights.begin(), heights.end());
  size_t hits = 1;
  for (size_t i = 1; i < heights.size(); i++)
    if (heights[i] - heights[i - 1] > c.tolerance)
      hits++;
  return hits % 2 == 1;
}
bool pointInside(const Context &c, Vec3 p) {
  if (!finite(p))
    throw std::invalid_argument("nonfinite point");
  if (!c.occupancy.empty()) {
    Cell ix = indices(p, c.lower, c.occupancyCellSize, c.occupancyCellCount);
    auto label = c.occupancy[offset(ix, c.occupancyCellCount)];
    if (label < 2)
      return label == 1;
  }
  return exactPointInside(c, p);
}
Context buildContext(Mesh mesh, double maxRadius, double buffer,
                     double relativeTolerance, const GeometryOptions &opt) {
  if (mesh.triangles.empty())
    throw std::invalid_argument("empty mesh");
  if (!std::isfinite(maxRadius) || maxRadius <= 0 || !std::isfinite(buffer) ||
      buffer < 0 || !std::isfinite(2 * maxRadius + buffer) ||
      !std::isfinite(relativeTolerance) || relativeTolerance < 0)
    throw std::invalid_argument("invalid geometry scale");
  if (!std::isfinite(opt.occupancyCellSize) || opt.occupancyCellSize < 0 ||
      opt.occupancyMaxCells == 0 ||
      opt.occupancyMaxCells > std::numeric_limits<std::uint32_t>::max() ||
      opt.threads < 1)
    throw std::invalid_argument("invalid geometry options");
  Context c;
  c.mesh = std::move(mesh);
  c.lower = {inf, inf, inf};
  c.upper = {-inf, -inf, -inf};
  for (auto &t : c.mesh.triangles)
    for (auto p : {t.a, t.b, t.c}) {
      if (!finite(p))
        throw std::invalid_argument("nonfinite mesh vertex");
      c.lower = min(c.lower, p);
      c.upper = max(c.upper, p);
    }
  Vec3 extent = c.upper - c.lower;
  if (!finite(extent))
    throw std::overflow_error("mesh extent overflow");
  double modelScale = largest(extent),
         magnitude = std::max(largest(abs3(c.lower)), largest(abs3(c.upper)));
  if (modelScale <= 0)
    throw std::invalid_argument("zero-size mesh");
  c.tolerance =
      std::max(relativeTolerance * modelScale, 64 * spacing(magnitude));
  if (!std::isfinite(c.tolerance))
    throw std::overflow_error("geometry tolerance overflow");
  c.cellSize = 2 * maxRadius + buffer;
  c.cellCount = countsFor(extent, c.cellSize);
  c.rayCellSize = c.cellSize;
  c.rayCellCount = {c.cellCount[0], c.cellCount[1], 1};
  c.rayTriangles.resize(c.mesh.triangles.size());
  for (size_t id = 0; id < c.mesh.triangles.size(); id++) {
    auto &t = c.mesh.triangles[id];
    Vec3 ab = t.b - t.a, ac = t.c - t.a;
    double det = ab.x * ac.y - ab.y * ac.x;
    if (!std::isfinite(det) || !finite(cross(ab, ac)) ||
        !std::isfinite(dot(ab, ab)) || !std::isfinite(dot(ac, ac)))
      throw std::overflow_error("mesh feature arithmetic overflow");
    raster(t, c.lower, c.cellSize, c.cellCount,
           (1 + std::sqrt(3.)) * c.cellSize / 2, c.tolerance, c.triangleCells,
           id);
    if (std::abs(det) > c.tolerance * c.tolerance) {
      auto [lo, hi] = rayBounds(t, c.tolerance, magnitude);
      c.rayTriangles[id] = {lo, hi, ab, ac, 1 / det};
      Cell l = indices(lo, c.lower, c.rayCellSize, c.rayCellCount),
           h = indices(hi, c.lower, c.rayCellSize, c.rayCellCount);
      l[2] = h[2] = 0;
      boxCells(l, h, [&](Cell a) { c.rayCells[a].push_back(id); });
    }
  }
  double probe = std::max(100 * c.tolerance, 1e-8 * c.cellSize);
  c.inwardNormals.reserve(c.mesh.triangles.size());
  for (auto &t : c.mesh.triangles) {
    Vec3 n = cross(t.b - t.a, t.c - t.a);
    double len = norm(n);
    if (len > 0)
      n = n / len;
    Vec3 centre = t.a + (t.b - t.a) / 3 + (t.c - t.a) / 3;
    if (!exactPointInside(c, centre + probe * n))
      n = -n;
    c.inwardNormals.push_back(n);
  }
  if (opt.occupancyAcceleration) {
    double width = opt.occupancyCellSize;
    if (width == 0)
      width = std::max(maxRadius / 2, 32 * c.tolerance);
    Cell n = countsFor(extent, width);
    auto totalOf = [](Cell a) {
      return static_cast<long double>(a[0]) * a[1] * a[2];
    };
    while (totalOf(n) > opt.occupancyMaxCells) {
      width *= std::ceil(std::cbrt(totalOf(n) / opt.occupancyMaxCells));
      if (!std::isfinite(width))
        throw std::overflow_error("occupancy spacing overflow");
      n = countsFor(extent, width);
    }
    size_t total = static_cast<size_t>(totalOf(n));
    c.occupancyCellSize = width;
    c.occupancyCellCount = n;
    c.occupancy.assign(total, 3);
    double halo = std::max(c.tolerance, 8 * spacing(magnitude));
    for (auto &t : c.mesh.triangles) {
      Vec3 pad{width + halo, width + halo, width + halo};
      Cell lo = indices(min(min(t.a, t.b), t.c) - pad, c.lower, width, n),
           hi = indices(max(max(t.a, t.b), t.c) + pad, c.lower, width, n);
      boxCells(lo, hi, [&](Cell a) { c.occupancy[offset(a, n)] = 2; });
    }
    std::vector<size_t> fifo;
    fifo.reserve(total);
    auto decode = [&](size_t k) {
      Cell a;
      a[0] = static_cast<std::int64_t>(k % static_cast<size_t>(n[0]));
      k /= static_cast<size_t>(n[0]);
      a[1] = static_cast<std::int64_t>(k % static_cast<size_t>(n[1]));
      a[2] = static_cast<std::int64_t>(k / static_cast<size_t>(n[1]));
      return a;
    };
    auto flood = [&](size_t seed) {
      fifo.clear();
      fifo.push_back(seed);
      c.occupancy[seed] = 0;
      for (size_t head = 0; head < fifo.size(); head++) {
        Cell a = decode(fifo[head]);
        for (int axis = 0; axis < 3; axis++)
          for (int direction : {-1, 1}) {
            Cell b = a;
            b[axis] += direction;
            if (b[axis] < 0 || b[axis] >= n[axis])
              continue;
            size_t j = offset(b, n);
            if (c.occupancy[j] == 3) {
              c.occupancy[j] = 0;
              fifo.push_back(j);
            }
          }
      }
    };
    for (size_t seed = 0; seed < total; seed++) {
      if (c.occupancy[seed] != 3)
        continue;
      Cell a = decode(seed);
      if (a[0] == 0 || a[1] == 0 || a[2] == 0 || a[0] == n[0] - 1 ||
          a[1] == n[1] - 1 || a[2] == n[2] - 1)
        flood(seed);
    }
    for (size_t seed = 0; seed < total; seed++) {
      if (c.occupancy[seed] != 3)
        continue;
      flood(seed);
      Cell a = decode(seed);
      Vec3 p = c.lower + Vec3{(a[0] + .5) * width, (a[1] + .5) * width,
                              (a[2] + .5) * width};
      if (exactPointInside(c, p))
        for (auto j : fifo)
          c.occupancy[j] = 1;
    }
  }
  return c;
}
double meshSignedVolume(const Mesh &m) {
  if (m.triangles.empty())
    return 0;
  Vec3 origin = m.triangles.front().a;
  long double sum = 0;
  for (auto &t : m.triangles) {
    Vec3 a = t.a - origin, b = t.b - origin, c = t.c - origin;
    sum +=
        static_cast<long double>(a.x) * (static_cast<long double>(b.y) * c.z -
                                         static_cast<long double>(b.z) * c.y) +
        static_cast<long double>(a.y) * (static_cast<long double>(b.z) * c.x -
                                         static_cast<long double>(b.x) * c.z) +
        static_cast<long double>(a.z) * (static_cast<long double>(b.x) * c.y -
                                         static_cast<long double>(b.y) * c.x);
  }
  double volume = static_cast<double>(sum / 6);
  if (!std::isfinite(volume))
    throw std::overflow_error("mesh volume overflow");
  return volume;
}
} // namespace sp
