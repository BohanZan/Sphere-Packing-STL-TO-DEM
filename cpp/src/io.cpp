#include "sp/io.hpp"
#include <algorithm>
#include <bit>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace sp {
namespace {
double number(std::string_view token) {
  if (!token.empty() && token.front() == '+')
    token.remove_prefix(1);
  double result = 0;
  const auto parsed =
      std::from_chars(token.data(), token.data() + token.size(), result);
  if (token.empty() || parsed.ec != std::errc{} ||
      parsed.ptr != token.data() + token.size() || !std::isfinite(result))
    throw std::runtime_error("Invalid finite number: " + std::string(token));
  return result;
}
std::vector<std::string_view> tokens(std::string_view line) {
  std::vector<std::string_view> result;
  while (!line.empty()) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == line.npos)
      break;
    line.remove_prefix(first);
    const auto end = line.find_first_of(" \t\r\n");
    result.push_back(line.substr(0, end));
    if (end == line.npos)
      break;
    line.remove_prefix(end);
  }
  return result;
}
std::uint32_t little32(const unsigned char *b) {
  return std::uint32_t(b[0]) | (std::uint32_t(b[1]) << 8) |
         (std::uint32_t(b[2]) << 16) | (std::uint32_t(b[3]) << 24);
}
std::ofstream output(const std::filesystem::path &path) {
  std::ofstream stream(path);
  stream.exceptions(std::ios::failbit | std::ios::badbit);
  stream.imbue(std::locale::classic());
  stream << std::setprecision(17);
  return stream;
}
std::string label(const Cell &c) {
  return std::to_string(c[0] + 1) + "," + std::to_string(c[1] + 1) + "," +
         std::to_string(c[2] + 1);
}
std::string jsonString(const std::string &value) {
  std::ostringstream out;
  out << '"';
  for (unsigned char c : value) {
    if (c == '"' || c == '\\')
      out << '\\' << c;
    else if (c < 32)
      out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c)
          << std::dec;
    else
      out << c;
  }
  out << '"';
  return out.str();
}
void vectorJson(std::ostream &out, Vec3 p) {
  out << '[' << p.x << ',' << p.y << ',' << p.z << ']';
}
} // namespace

Mesh readStl(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("Cannot open STL: " + path.string());
  const auto bytes = std::filesystem::file_size(path);
  unsigned char header[84]{};
  in.read(reinterpret_cast<char *>(header), 84);
  const auto count = little32(header + 80);
  const bool binary = bytes >= 84 && bytes == 84ULL + 50ULL * count;
  in.clear();
  in.seekg(binary ? 84 : 0);
  Mesh mesh;
  if (binary) {
    mesh.triangles.reserve(count);
    for (std::uint32_t id = 0; id < count; ++id) {
      unsigned char record[50];
      in.read(reinterpret_cast<char *>(record), 50);
      if (!in)
        throw std::runtime_error("Truncated binary STL triangle");
      double values[12];
      for (int i = 0; i < 12; ++i) {
        values[i] = std::bit_cast<float>(little32(record + 4 * i));
        if (!std::isfinite(values[i]))
          throw std::runtime_error("Nonfinite binary STL number");
      }
      mesh.triangles.push_back({{values[3], values[4], values[5]},
                                {values[6], values[7], values[8]},
                                {values[9], values[10], values[11]}});
    }
  } else {
    // Strict line grammar avoids silently accepting truncated faces or
    // atof("bad")==0.
    int stage = 0;
    bool inSolid = false, ended = false;
    std::string line;
    Triangle tri;
    std::size_t lineId = 0;
    while (std::getline(in, line)) {
      ++lineId;
      const auto t = tokens(line);
      if (t.empty())
        continue;
      auto syntax = [&]() {
        throw std::runtime_error("Malformed/truncated ASCII STL at line " +
                                 std::to_string(lineId));
      };
      if (t[0] == "solid") {
        if (stage != 0 || inSolid)
          syntax();
        inSolid = true;
        ended = false;
        continue;
      }
      if (t[0] == "endsolid") {
        if (stage != 0 || !inSolid)
          syntax();
        inSolid = false;
        ended = true;
        continue;
      }
      if (!inSolid)
        syntax();
      if (stage == 0) {
        if (t.size() != 5 || t[0] != "facet" || t[1] != "normal")
          syntax();
        for (int i = 2; i < 5; ++i)
          (void)number(t[i]);
        stage = 1;
      } else if (stage == 1) {
        if (t.size() != 2 || t[0] != "outer" || t[1] != "loop")
          syntax();
        stage = 2;
      } else if (stage >= 2 && stage <= 4) {
        if (t.size() != 4 || t[0] != "vertex")
          syntax();
        Vec3 p{number(t[1]), number(t[2]), number(t[3])};
        if (stage == 2)
          tri.a = p;
        else if (stage == 3)
          tri.b = p;
        else
          tri.c = p;
        ++stage;
      } else if (stage == 5) {
        if (t.size() != 1 || t[0] != "endloop")
          syntax();
        stage = 6;
      } else {
        if (t.size() != 1 || t[0] != "endfacet")
          syntax();
        mesh.triangles.push_back(tri);
        stage = 0;
      }
    }
    if (in.bad() || stage != 0 || inSolid || !ended)
      throw std::runtime_error("Truncated or invalid STL file: " +
                               path.string());
  }
  if (mesh.triangles.empty())
    throw std::runtime_error("STL contains no triangles");
  return mesh;
}

std::vector<double> readNumbers(const std::filesystem::path &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open numeric file: " + path.string());
  std::vector<double> values;
  std::string line;
  while (std::getline(in, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    for (auto t : tokens(line))
      values.push_back(number(t));
  }
  if (in.bad())
    throw std::runtime_error("Error reading numeric file");
  return values;
}

void writeOutputs(const std::filesystem::path &directory,
                  const std::string &prefix, const Context &c, Result &result) {
  if (prefix.empty() || prefix == "." || prefix == ".." ||
      prefix.find_first_of("/\\:\r\n") != prefix.npos)
    throw std::invalid_argument(
        "Output prefix must be a filename without a directory");
  std::filesystem::create_directories(directory);
  const auto spheres = directory / (prefix + "_spheres.csv"),
             summary = directory / (prefix + "_summary.csv"),
             points = directory / (prefix + "_grid_points.csv"),
             cells = directory / (prefix + "_grid_hexahedra.csv");
  auto sf = output(spheres);
  sf << "id,x,y,z,radius,diameter,mass\n";
  for (std::size_t i = 0; i < result.state.radii.size(); ++i) {
    auto p = result.outputCentres.at(i);
    double r = result.state.radii[i];
    sf << i + 1 << ',' << p.x << ',' << p.y << ',' << p.z << ',' << r << ','
       << 2 * r << ',' << result.masses.at(i) << '\n';
  }
  sf.close();
  const auto &r = result.report;
  auto sm = output(summary);
  sm << "requested_count,accepted_count,unplaced_count,stop_reason,capacity_"
        "warning,total_volume,inertia_xx,inertia_yy,inertia_zz\n"
     << r.requestedCount << ',' << r.acceptedCount << ',' << r.unplacedCount
     << ',' << r.stopReason << ',' << r.capacityWarning << ','
     << result.totalVolume << ',' << result.inertia[0] << ','
     << result.inertia[4] << ',' << result.inertia[8] << '\n';
  sm.close();
  std::map<std::string, Cell> occupied;
  for (const auto &[cell, ids] : c.triangleCells)
    if (!ids.empty())
      occupied.emplace(label(cell), cell);
  for (const auto &[cell, ids] : result.state.sphereCells)
    if (!ids.empty())
      occupied.emplace(label(cell), cell);
  auto pf = output(points), cf = output(cells);
  pf << "point_id,x,y,z\n";
  cf << "cell_id,p1,p2,p3,p4,p5,p6,p7,p8,sphere_count,triangle_count\n";
  std::size_t id = 0;
  for (const auto &[name, index] : occupied) {
    (void)name;
    Vec3 lo =
        c.lower +
        Vec3{double(index[0]), double(index[1]), double(index[2])} * c.cellSize;
    Vec3 hi = min(lo + Vec3{c.cellSize, c.cellSize, c.cellSize}, c.upper);
    lo = lo - r.coordinateShift;
    hi = hi - r.coordinateShift;
    cf << id + 1;
    for (int corner = 0; corner < 8; ++corner) {
      const auto pid = 8 * id + corner + 1;
      pf << pid << ',' << ((corner & 1) ? hi.x : lo.x) << ','
         << ((corner & 2) ? hi.y : lo.y) << ',' << ((corner & 4) ? hi.z : lo.z)
         << '\n';
      cf << ',' << pid;
    }
    cf << ',' << lookup(result.state.sphereCells, index).size() << ','
       << lookup(c.triangleCells, index).size() << '\n';
    ++id;
  }
  pf.close();
  cf.close();
  result.report.outputFiles = {spheres.string(), summary.string(),
                               points.string(), cells.string()};
}

void printSummary(std::ostream &out, const Context &c, const Result &result) {
  const auto &r = result.report;
  auto d = r.boundingBoxDimensions;
  out << std::setprecision(8) << "Bounding Box Dimensions Lx=" << d.x
      << "; Ly=" << d.y << "; Lz=" << d.z
      << "\nSpatial Grid Discretisation\nGrid Cell Edge Length = " << c.cellSize
      << "\nGrid Counts Nx=" << c.cellCount[0] << "; Ny=" << c.cellCount[1]
      << "; Nz=" << c.cellCount[2] << "\nSTL Volume: " << r.stlVolume
      << ", Sphere Assembly Volume: " << result.totalVolume
      << "\nSum of sphere masses: " << r.totalMass
      << ", total mass: " << r.totalMass
      << "\nCoM of the cluster: " << r.centreOfMass.x << ' ' << r.centreOfMass.y
      << ' ' << r.centreOfMass.z
      << "\nCoM of the cluster after shifting: " << r.centreOfMassAfterShift.x
      << ' ' << r.centreOfMassAfterShift.y << ' ' << r.centreOfMassAfterShift.z
      << "\n========================================\nFinished\nNumber of "
         "final spheres = "
      << r.acceptedCount << "\nAssembly volume = " << result.totalVolume
      << "\nMoment of inertia tensor:\n";
  for (int row = 0; row < 3; ++row)
    out << result.inertia[3 * row] << ' ' << result.inertia[3 * row + 1] << ' '
        << result.inertia[3 * row + 2] << '\n';
  out << "Status: " << r.stopReason << "; unplaced=" << r.unplacedCount
      << "\n========================================\nCSV result files:\n";
  for (const auto &f : r.outputFiles)
    out << "  " << f << '\n';
}

void writeReportJson(const std::filesystem::path &path, const Result &result,
                     double importSeconds, double preprocessingSeconds,
                     double packingSeconds, double outputSeconds,
                     const Random &random) {
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
  auto out = output(path);
  const auto &r = result.report;
  out << "{\n\"requestedCount\":" << r.requestedCount
      << ",\n\"acceptedCount\":" << r.acceptedCount
      << ",\n\"unplacedCount\":" << r.unplacedCount
      << ",\n\"nextUnplacedRadiusIndex\":" << r.nextUnplacedRadiusIndex
      << ",\n\"initialFailures\":" << r.initialFailures
      << ",\n\"refillPasses\":" << r.refillPasses
      << ",\n\"stopReason\":" << jsonString(r.stopReason)
      << ",\n\"capacityWarning\":" << (r.capacityWarning ? "true" : "false")
      << ",\n\"coordinateFrame\":" << jsonString(r.coordinateFrame)
      << ",\n\"totalVolume\":" << result.totalVolume
      << ",\n\"totalMass\":" << r.totalMass
      << ",\n\"stlVolume\":" << r.stlVolume << ",\n\"boundingBoxDimensions\":";
  vectorJson(out, r.boundingBoxDimensions);
  out << ",\n\"sphereAssemblyVolume\":" << r.sphereAssemblyVolume;
  out << ",\n\"centreOfMass\":";
  vectorJson(out, r.centreOfMass);
  out << ",\n\"coordinateShift\":";
  vectorJson(out, r.coordinateShift);
  out << ",\n\"centreOfMassAfterShift\":";
  vectorJson(out, r.centreOfMassAfterShift);
  out << ",\n\"inertia\":[";
  for (std::size_t i = 0; i < 9; ++i) {
    if (i)
      out << ',';
    out << result.inertia[i];
  }
  out << "],\n\"outputFiles\":[";
  for (std::size_t i = 0; i < r.outputFiles.size(); ++i) {
    if (i)
      out << ',';
    out << jsonString(r.outputFiles[i]);
  }
  out << "],\n\"uniformDraws\":" << random.uniformCount
      << ",\n\"normalDraws\":" << random.normalCount
      << ",\n\"importSeconds\":" << importSeconds
      << ",\n\"preprocessingSeconds\":" << preprocessingSeconds
      << ",\n\"packingSeconds\":" << packingSeconds
      << ",\n\"outputSeconds\":" << outputSeconds
      << ",\n\"computeSeconds\":" << preprocessingSeconds + packingSeconds
      << ",\n\"totalSeconds\":"
      << importSeconds + preprocessingSeconds + packingSeconds + outputSeconds
      << "\n}\n";
  out.close();
}
} // namespace sp
