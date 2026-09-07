#include "sp/io.hpp"
#include <bit>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace sp;
void require(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
template <class F> void rejected(F f) {
  bool yes = false;
  try {
    f();
  } catch (const std::exception &) {
    yes = true;
  }
  require(yes, "Malformed input must be rejected");
}
int main() {
  try {
    auto dir =
        std::filesystem::temp_directory_path() /
        ("sp_io_" +
         std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    struct Cleanup {
      std::filesystem::path p;
      ~Cleanup() {
        std::error_code ec;
        std::filesystem::remove_all(p, ec);
      }
    } cleanup{dir};
    const auto text = dir / "test.stl";
    {
      std::ofstream f(text);
      f << "solid s\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 "
           "0\nvertex 0 1 0\nendloop\nendfacet\nendsolid s\n";
    }
    auto mesh = readStl(text);
    require(mesh.triangles.size() == 1, "ASCII triangle count");
    require(mesh.triangles[0].b.x == 1 && mesh.triangles[0].c.y == 1,
            "ASCII order/precision");
    {
      std::ofstream f(text);
      f << "solid s\nfacet normal 0 0 1\nouter loop\nvertex bad 0 0\nvertex 1 "
           "0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid s\n";
    }
    rejected([&] { readStl(text); });
    {
      std::ofstream f(text);
      f << "solid s\nfacet normal 0 0 1\nouter loop\nvertex nan 0 0\nvertex 1 "
           "0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid s\n";
    }
    rejected([&] { readStl(text); });
    {
      std::ofstream f(text);
      f << "solid s\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\n";
    }
    rejected([&] { readStl(text); });
    const auto binary = dir / "solid-header.stl";
    {
      std::ofstream f(binary, std::ios::binary);
      std::string head(80, ' ');
      head.replace(0, 5, "solid");
      f.write(head.data(), 80);
      unsigned char count[4] = {1, 0, 0, 0};
      f.write(reinterpret_cast<char *>(count), 4);
      float data[12] = {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0};
      f.write(reinterpret_cast<char *>(data), 48);
      f.put(0);
      f.put(0);
    }
    mesh = readStl(binary);
    require(mesh.triangles.size() == 1 && mesh.triangles[0].b.x == 1,
            "Binary solid-header disambiguation");
    std::filesystem::resize_file(binary, 133);
    rejected([&] { readStl(binary); });
    {
      std::ofstream f(text);
      f << "solid s\nfacet normal 0 0 0\nouter loop\nvertex 0 0 0\nvertex 0 0 "
           "0\nvertex 1 0 0\nendloop\nendfacet\nendsolid s\n";
    }
    require(readStl(text).triangles.size() == 1,
            "Keep degenerate facets instead of silently editing mesh");
    auto numbers = dir / "radii.txt";
    {
      std::ofstream f(numbers);
      f << "0.1,0.2\n0.3\n";
    }
    require(readNumbers(numbers).size() == 3, "Ordered radius sequence");
    {
      std::ofstream f(numbers);
      f << "0.2 nope";
    }
    rejected([&] { readNumbers(numbers); });
    Context c;
    c.lower = {0, 0, 0};
    c.upper = {2, 2, 2};
    c.cellSize = 1;
    c.cellCount = {2, 2, 2};
    c.triangleCells[{0, 0, 0}] = {0, 1};
    Result r;
    r.state.centres = {{.5, .5, .5}};
    r.state.radii = {.1};
    r.state.sphereCells[{0, 0, 0}] = {0};
    r.outputCentres = {{0, 0, 0}};
    r.masses = {1};
    r.report.coordinateShift = {.5, .5, .5};
    r.report.acceptedCount = 1;
    writeOutputs(dir, "check", c, r);
    require(r.report.outputFiles.size() == 4, "Exactly four CSV outputs");
    std::ifstream sf(dir / "check_spheres.csv");
    std::string line;
    std::getline(sf, line);
    require(line == "id,x,y,z,radius,diameter,mass", "Sphere CSV schema");
    std::ifstream gf(dir / "check_grid_points.csv");
    std::getline(gf, line);
    std::getline(gf, line);
    require(line == "1,-0.5,-0.5,-0.5", "Grid coordinate shift");
    rejected([&] { writeOutputs(dir, "../escape", c, r); });
    r.report.sphereAssemblyVolume = 1.25;
    Random random;
    writeReportJson(dir / "report.json", r, 0, 0, 0, 0, random);
    std::ifstream report(dir / "report.json");
    std::string json((std::istreambuf_iterator<char>(report)), {});
    require(json.find("\"sphereAssemblyVolume\":1.25") != json.npos,
            "Full MATLAB report includes sphereAssemblyVolume");
    std::cout << "IO regression tests passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
