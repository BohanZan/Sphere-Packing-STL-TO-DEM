#include "sp/geometry.hpp"
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace sp;
void check(bool x) {
  if (!x)
    throw std::runtime_error("geometry assertion");
}
Mesh box(double s) {
  Vec3 v[] = {{0, 0, 0}, {s, 0, 0}, {s, s, 0}, {0, s, 0},
              {0, 0, s}, {s, 0, s}, {s, s, s}, {0, s, s}};
  int f[][3] = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
                {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
  Mesh m;
  for (auto &a : f)
    m.triangles.push_back({v[a[0]], v[a[1]], v[a[2]]});
  return m;
}
void regressions() {
  Triangle tri{{0, 0, 0}, {4, 0, 0}, {0, 4, 0}};
  for (double s : {1e-100, 1e100}) {
    Triangle scaled{{}, s * tri.b, s * tri.c};
    check(std::abs(pointTriangleDistance(s * Vec3{1, 1, 2}, scaled) / s - 2) <
          1e-12);
    check(
        std::abs(triangleContactDistance(s * Vec3{1, 1, 2}, .5 * s, {0, 0, -1},
                                         scaled, 10 * s, 1e-10 * s) /
                     s -
                 1.5) < 1e-12);
  }
  Vec3 p[] = {{1, 1, 2},   {1, -1, 1}, {-1, -1, 1}, {2, -2, 0},
              {-1, 1, .5}, {2, -1, 2}, {1, -1, .1}};
  Vec3 u[] = {{0, 0, -1},
              {0, 1, -1},
              {1, 1, -1},
              {0, 1, 0},
              {1, 0, 0},
              {0, 0, -1},
              {0, std::sqrt(.99), .1}};
  double expected[] = {1.5, std::sqrt(2.) - .5, std::sqrt(3.) - .5, 1.5,
                       1,   INFINITY,           .5265114276384391};
  for (int k = 0; k < 7; k++)
    for (int winding = 0; winding < 2; winding++) {
      Triangle t = tri;
      if (winding)
        std::swap(t.b, t.c);
      double actual =
          triangleContactDistance(p[k], .5, u[k] / norm(u[k]), t, 10, 1e-10);
      check((std::isinf(actual) && std::isinf(expected[k])) ||
            std::abs(actual - expected[k]) < 1e-9);
    }
  check(std::isinf(triangleContactDistance({-1e8, -.5001, 0}, .5, {1, 0, 0},
                                           tri, 2e8, 1e-9)));
  Triangle longEdge{{0, 0, 0}, {50000, 120000, 0}, {-12, 5, 0}};
  check(std::abs(
            triangleContactDistance({25001.846153846152, 59999.230769230766, 4},
                                    2, {0, 0, -1}, longEdge, 8, 1e-9) -
            4) < 1e-5);
  for (double s : {1e-4, 1., 1e4}) {
    Vec3 shift = s * Vec3{1000, -2000, 3000};
    Triangle t{shift, s * tri.b + shift, s * tri.c + shift};
    check(std::abs(triangleContactDistance(s * Vec3{1, -1, 1} + shift, .5 * s,
                                           Vec3{0, 1, -1} / std::sqrt(2.), t,
                                           10 * s, 1e-10 * s) -
                   (std::sqrt(2.) - .5) * s) < 1e-8 * s);
  }
  std::mt19937 rng(872);
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> uniform;
  auto randomVec = [&]() {
    return Vec3{normal(rng), normal(rng), normal(rng)};
  };
  // Independent convex distance minimization followed by bisection verifies
  // first entry.
  for (int k = 0; k < 150; k++) {
    Triangle t{randomVec(), randomVec(), randomVec()};
    Vec3 p = 3 * randomVec();
    double r = .1 + .4 * uniform(rng);
    if (pointTriangleDistance(p, t) <= r + 1e-6)
      continue;
    Vec3 u = (t.a + t.b + t.c) / 3 - p + .8 * randomVec();
    u = u / norm(u);
    double actual = triangleContactDistance(p, r, u, t, 10, 1e-10);
    auto d = [&](double x) { return pointTriangleDistance(p + x * u, t); };
    double a = 0, b = 10;
    for (int j = 0; j < 120; j++) {
      double x = a + (b - a) / 3, y = b - (b - a) / 3;
      if (d(x) < d(y))
        b = y;
      else
        a = x;
    }
    double tm = (a + b) / 2, dm = d(tm);
    if (dm < r - 1e-7) {
      a = 0;
      b = tm;
      for (int j = 0; j < 60; j++) {
        double m = (a + b) / 2;
        if (d(m) > r)
          a = m;
        else
          b = m;
      }
      check(std::abs(actual - b) < 1e-7);
    } else if (dm > r + 1e-7)
      check(std::isinf(actual));
  }
  GeometryOptions opt;
  opt.occupancyAcceleration = false;
  Mesh thin;
  thin.triangles.push_back({{1, 1, 1}, {98, 98, 98}, {49, 49 + 1e-10, 49}});
  auto thinContext = buildContext(thin, .5, 0, 1e-12, opt);
  check(thinContext.triangleCells.size() < 20000);
  Mesh oblique;
  oblique.triangles.push_back({{1, 1, 1}, {78, 78, 3}, {3, 76, 78}});
  auto c = buildContext(oblique, .5, 0, 1e-12, opt);
  check(c.triangleCells.size() < 80000);
  for (int i = 0; i < 5000; i++) {
    Vec3 p =
        c.lower + Vec3{uniform(rng) * 77, uniform(rng) * 77, uniform(rng) * 77};
    auto ix = cellIndex(c, p);
    Vec3 centre =
        c.lower + Vec3{(ix[0] + .5) * c.cellSize, (ix[1] + .5) * c.cellSize,
                       (ix[2] + .5) * c.cellSize};
    if (pointTriangleDistance(centre, oblique.triangles[0]) <=
        (1 + std::sqrt(3.)) * c.cellSize / 2)
      check(!lookup(c.triangleCells, ix).empty());
  }
  Mesh cavity = box(4), inner = box(2);
  for (auto t : inner.triangles) {
    t.a = t.a + Vec3{1, 1, 1};
    t.b = t.b + Vec3{1, 1, 1};
    t.c = t.c + Vec3{1, 1, 1};
    std::swap(t.b, t.c);
    cavity.triangles.push_back(t);
  }
  opt.occupancyAcceleration = true;
  opt.occupancyCellSize = .125;
  auto hollow = buildContext(cavity, .25, 0, 1e-9, opt);
  check(!pointInside(hollow, {2, 2, 2}));
  check(pointInside(hollow, {.5, .5, .5}));
  Mesh rotated = box(4);
  double theta = .63;
  for (auto &t : rotated.triangles)
    for (Vec3 *p : {&t.a, &t.b, &t.c}) {
      Vec3 a = *p;
      *p = {std::cos(theta) * a.x + std::sin(theta) * a.z, a.y,
            -std::sin(theta) * a.x + std::cos(theta) * a.z};
    }
  auto rc = buildContext(rotated, .25, 0, 1e-9, opt);
  for (int k = 0; k < 1000; k++) {
    Vec3 e = rc.upper - rc.lower,
         p = rc.lower +
             Vec3{uniform(rng) * e.x, uniform(rng) * e.y, uniform(rng) * e.z};
    check(pointInside(rc, p) == exactPointInside(rc, p));
  }
}
int main() {
  try {
    regressions();
    Triangle t{{0, 0, 0}, {2, 0, 0}, {0, 2, 0}};
    check(std::abs(pointTriangleDistance({.5, .5, 3}, t) - 3) < 1e-12);
    check(std::abs(
              triangleContactDistance({.5, .5, 3}, 1, {0, 0, -1}, t, 9, 1e-10) -
              2) < 1e-12);
    check(std::isinf(
        triangleContactDistance({.5, .5, 1}, 1, {1, 0, 0}, t, 9, 1e-10)));
    check(triangleContactDistance({.5, .5, 1}, 1, {0, 0, -1}, t, 9, 1e-10) ==
          0);
    check(std::abs(
              triangleContactDistance({3, 0, 0}, .5, {-1, 0, 0}, t, 9, 1e-10) -
              .5) < 1e-12);
    for (double s : {1e-6, 1., 1e6}) {
      GeometryOptions opt;
      opt.occupancyMaxCells = 50000;
      auto c = buildContext(box(s), s / 20, 0, 1e-10, opt);
      check(exactPointInside(c, {s / 2, s / 2, s / 2}));
      check(!exactPointInside(c, {s / 2, s / 2, 2 * s}));
      check(!pointInside(c, {-s, s / 2, s / 2}));
      for (int x = 0; x < 13; x++)
        for (int y = 0; y < 13; y++)
          for (int z = 0; z < 13; z++) {
            Vec3 p{s * (x + .13) / 13, s * (y + .27) / 13, s * (z + .39) / 13};
            check(pointInside(c, p) == exactPointInside(c, p));
          }
      check(!lookup(c.triangleCells, cellIndex(c, {s * .5, s * .5, s * .01}))
                 .empty());
      for (size_t i = 0; i < c.mesh.triangles.size(); i++) {
        auto a = c.mesh.triangles[i];
        Vec3 p = (a.a + a.b + a.c) / 3 + c.inwardNormals[i] * (s * 1e-6);
        check(exactPointInside(c, p));
      }
    }
    bool thrown = false;
    try {
      buildContext(box(1), 0, 0, 1e-10);
    } catch (const std::exception &) {
      thrown = true;
    }
    check(thrown);
    std::cout << "geometry tests passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
