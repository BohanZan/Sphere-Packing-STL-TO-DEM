#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
namespace sp {
struct Vec3 {
  double x = 0, y = 0, z = 0;
  double &operator[](std::size_t i) { return i == 0 ? x : i == 1 ? y : z; }
  double operator[](std::size_t i) const { return i == 0 ? x : i == 1 ? y : z; }
};
Vec3 operator+(Vec3, Vec3);
Vec3 operator-(Vec3, Vec3);
Vec3 operator-(Vec3);
Vec3 operator*(Vec3, double);
Vec3 operator*(double, Vec3);
Vec3 operator/(Vec3, double);
double dot(Vec3, Vec3);
Vec3 cross(Vec3, Vec3);
double norm(Vec3);
Vec3 min(Vec3, Vec3);
Vec3 max(Vec3, Vec3);
struct Triangle {
  Vec3 a, b, c;
};
struct Mesh {
  std::vector<Triangle> triangles;
};
using Cell = std::array<std::int64_t, 3>;
struct CellHash {
  std::size_t operator()(const Cell &) const noexcept;
};
using Grid = std::unordered_map<Cell, std::vector<std::size_t>, CellHash>;
struct GeometryOptions {
  bool occupancyAcceleration = true;
  double occupancyCellSize = 0;
  std::size_t occupancyMaxCells = 2000000;
  int threads = 1;
};
struct RayTriangle {
  Vec3 lower, upper, ab, ac;
  double inverseDeterminant = 0;
};
struct Context {
  Mesh mesh;
  Vec3 lower, upper;
  double cellSize = 0, tolerance = 0;
  Cell cellCount{};
  Grid triangleCells;
  std::vector<Vec3> inwardNormals;
  Grid rayCells;
  std::vector<RayTriangle> rayTriangles;
  double rayCellSize = 0;
  Cell rayCellCount{};
  double occupancyCellSize = 0;
  Cell occupancyCellCount{};
  std::vector<std::int8_t> occupancy;
};
Context buildContext(Mesh mesh, double maxRadius, double buffer,
                     double relativeTolerance,
                     const GeometryOptions &options = {});
Cell cellIndex(const Context &, Vec3);
bool exactPointInside(const Context &, Vec3);
bool pointInside(const Context &, Vec3);
Vec3 closestPoint(Vec3, const Triangle &);
double pointTriangleDistance(Vec3, const Triangle &);
double triangleContactDistance(Vec3 p, double r, Vec3 u, const Triangle &,
                               double limit, double tolerance);
const std::vector<std::size_t> &lookup(const Grid &, const Cell &);
double meshSignedVolume(const Mesh &);
} // namespace sp
