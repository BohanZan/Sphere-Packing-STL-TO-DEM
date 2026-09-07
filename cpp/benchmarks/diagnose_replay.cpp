#include "sp/io.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace sp;
void runRelax(const Context &c, State &s, const Options &o, bool global,
              Random &random, std::ofstream &out) {
  auto f = gravityFrame(c, o.gravity);
  auto energy = [&]() {
    double e = 0;
    for (std::size_t i = 0; i < s.radii.size(); ++i)
      e += dot(s.centres[i] - f.origin, f.up) * s.radii[i] * s.radii[i] *
           s.radii[i];
    return e;
  };
  double previous = energy();
  for (std::size_t sweep = 1; sweep <= o.maxCompressionSweeps; ++sweep) {
    auto shakes = global ? 0 : o.shakeSweeps;
    for (std::size_t shake = 1; shake <= shakes + 1; ++shake)
      for (std::size_t id = 0; id < s.radii.size(); ++id) {
        Vec3 u = f.direction;
        if (shake <= shakes) {
          u = {random.normal(), random.normal(), random.normal()};
          u = u - f.direction * dot(u, f.direction);
          u = u / norm(u);
        }
        Vec3 p = s.centres[id];
        double d = firstContactDistance(c, s, id, u);
        if (global && sweep == 1 && id == 6) {
          std::ofstream detail("tmp/cpp_first_contact.csv");
          detail << std::setprecision(17);
          for (std::size_t face = 0; face < c.mesh.triangles.size(); ++face) {
            const auto &tri = c.mesh.triangles[face];
            double hit = triangleContactDistance(p, s.radii[id], u, tri,
                                                 d + 1e-8, c.tolerance);
            if (std::isfinite(hit) && std::abs(hit - d) < 1e-8) {
              Vec3 normal = cross(tri.b - tri.a, tri.c - tri.a);
              double len = norm(normal);
              Vec3 n = normal / len;
              detail << face + 1 << ',' << hit << ',' << normal.x << ','
                     << normal.y << ',' << normal.z << ',' << len << ',' << n.x
                     << ',' << n.y << ',' << n.z << ',' << dot(p - tri.a, n)
                     << ',' << dot(u, n) << '\n';
            }
          }
        }
        out << global << ',' << sweep << ',' << shake << ',' << id + 1 << ','
            << p.x << ',' << p.y << ',' << p.z << ',' << u.x << ',' << u.y
            << ',' << u.z << ',' << d << '\n';
        if (d > 0) {
          s.centres[id] = p + u * d;
          reindex(c, s, id);
        }
      }
    double e = energy();
    if (std::abs(e - previous) <=
        o.compressionTolerance *
            std::max(std::abs(previous), std::numeric_limits<double>::min()))
      return;
    previous = e;
  }
}
void validateSameState(const Context &c) {
  std::ifstream input("tmp/matlab_moves.csv");
  std::vector<std::array<double, 11>> rows;
  std::string line;
  while (std::getline(input, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream values(line);
    std::array<double, 11> row{};
    for (auto &value : row)
      if (!(values >> value))
        throw std::runtime_error("Invalid MATLAB move record");
    rows.push_back(row);
  }
  if (rows.size() < 50)
    throw std::runtime_error("Missing MATLAB trace");
  std::vector<Vec3> after(rows.size());
  std::array<Vec3, 50> nextPosition{};
  std::array<bool, 50> known{};
  for (std::size_t k = rows.size(); k-- > 0;) {
    const auto &row = rows[k];
    auto id = static_cast<std::size_t>(row[3]) - 1;
    Vec3 p{row[4], row[5], row[6]}, u{row[7], row[8], row[9]};
    after[k] = known[id] ? nextPosition[id] : p + u * row[10];
    nextPosition[id] = p;
    known[id] = true;
  }
  State state;
  for (auto p : nextPosition)
    addSphere(c, state, p, .625);
  double worst = 0;
  std::size_t worstIndex = 0, above = 0;
  std::ofstream errors("tmp/same_state_errors.csv");
  errors << std::setprecision(17);
  for (std::size_t k = 0; k < rows.size(); ++k) {
    const auto &row = rows[k];
    auto id = static_cast<std::size_t>(row[3]) - 1;
    state.centres[id] = {row[4], row[5], row[6]};
    reindex(c, state, id);
    Vec3 u{row[7], row[8], row[9]};
    double actual = firstContactDistance(c, state, id, u),
           error = std::abs(actual - row[10]);
    if (error > worst) {
      worst = error;
      worstIndex = k;
    }
    if (error > 1e-10)
      ++above;
    if (error > 1e-12)
      errors << k + 1 << ',' << row[10] << ',' << actual << ',' << error
             << '\n';
    state.centres[id] = after[k];
    reindex(c, state, id);
  }
  std::ofstream summary("tmp/same_state_summary.txt");
  summary << std::setprecision(17) << "queries=" << rows.size()
          << " maxAbsoluteError=" << worst << " row1based=" << worstIndex + 1
          << " above1e-10=" << above << '\n';
  std::cout << std::setprecision(17) << "Same-state contacts: " << rows.size()
            << " max error " << worst << " at row " << worstIndex + 1
            << " above1e-10 " << above << '\n';
  if (above)
    throw std::runtime_error("Same-state contact discrepancy exceeds 1e-10");
}
int main(int argc, char **) {
  auto mesh = readStl("inputs/greatBudda/greatBudda.stl");
  auto c = buildContext(std::move(mesh), .625, 0, 1e-9);
  if (argc > 1) {
    validateSameState(c);
    return 0;
  }
  Options o;
  o.maxAttempts = 60;
  o.maxCompressionSweeps = 200;
  o.shakeSweeps = 5;
  o.compressionTolerance = 1e-7;
  o.quiet = true;
  Random r;
  r.loadReplay(readNumbers("cpp/benchmarks/output/buddha/uniforms.txt"),
               readNumbers("cpp/benchmarks/output/buddha/normals.txt"));
  State s;
  std::size_t next = 0;
  initialPlacement(c, s, std::vector<double>(50, .625), next, o, r);
  std::ofstream out("tmp/cpp_moves.csv");
  out << std::setprecision(17);
  runRelax(c, s, o, true, r, out);
  runRelax(c, s, o, false, r, out);
}
