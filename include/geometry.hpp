#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
namespace sp {
// Vec3 is the small value type used by all geometric queries. Keeping it
// trivially movable makes the hot query loops cheap while retaining readable
// component-wise expressions in the implementation.
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
  // STL facets are kept as finite triangles, including degenerate facets.
  // Their original order is preserved for deterministic parity and output.
  Vec3 a, b, c;
};
struct Mesh {
  std::vector<Triangle> triangles;
};
using Cell = std::array<std::int64_t, 3>;
struct CellHash {
  // Hash sparse integer cell coordinates without allocating the full grid.
  std::size_t operator()(const Cell &) const noexcept;
};
using Grid = std::unordered_map<Cell, std::vector<std::size_t>, CellHash>;
struct GeometryOptions {
  // Occupancy is only a conservative shortcut; exact ray parity remains the
  // authority for margin cells and for contexts where the shortcut is off.
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
  // Immutable geometry and spatial indices shared by placement and relaxation.
  // Triangle and sphere grids use the same origin and contact cell spacing;
  // the ray grid stores the projected 2-D candidates separately.
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
// BUILD CONTEXT Read-only preprocessing for all sphere, triangle and ray
// queries. The main cell edge is at least the largest sphere diameter plus the
// movement buffer, so a local neighbourhood is complete for contact tests.
Context buildContext(Mesh mesh, double maxRadius, double buffer,
                     double relativeTolerance,
                     const GeometryOptions &options = {});
// Map a physical point to a clamped zero-based sparse-grid cell.
Cell cellIndex(const Context &, Vec3);
// Apply the exact downward-ray parity test to a point in the closed mesh.
bool exactPointInside(const Context &, Vec3);
// Use a conservative occupancy label when safe, otherwise call exact parity.
bool pointInside(const Context &, Vec3);
// Return the closest point on a finite triangle, including edge and vertex
// regions; repeated or collinear vertices are handled without division by 0.
Vec3 closestPoint(Vec3, const Triangle &);
// Euclidean distance to a finite triangle, not to its infinite supporting plane.
double pointTriangleDistance(Vec3, const Triangle &);
// First translation distance at which a moving sphere touches a finite triangle.
// The result includes face, edge and vertex contacts along the complete path.
double triangleContactDistance(Vec3 p, double r, Vec3 u, const Triangle &,
                               double limit, double tolerance);
// Return an empty static vector for an unoccupied sparse cell.
const std::vector<std::size_t> &lookup(const Grid &, const Cell &);
// Calculate the oriented STL volume using triangle tetrahedra about one anchor.
double meshSignedVolume(const Mesh &);
} // namespace sp
