# C++ port progress

Plan: CPP_PLAN.md
Baseline: main-cpp, user modifications already present; MATLAB .m files retained intact.
Environment detected: Intel oneAPI 2025.2, MATLAB R2025a, CMake, Ninja, VS 2022 and VS 18 Insiders.
Interface ownership: geometry agent owns geometry.hpp/cpp and geometry tests; packing agent owns packing.hpp/cpp and packing tests; controller owns IO, build, benchmark and docs. Geometry interface is agreed before integration.
Ruling: use the user-created main-cpp directory directly, as requested. Existing MATLAB changes are the source of truth.
Ruling: automatic execution authorization overrides skill requests for repeated design/phase approval.

Geometry: complete, independent review clean. Sparse raster/parity/occupancy and continuous finite-feature tests pass.
Packing: complete. Independent review found and fixed intrinsic inertia cancellation; failing then passing large-separation regression.
IO/build: complete. Strict ASCII/binary STL reader preserves triangles. Full report sphereAssemblyVolume fixed after review. Four CSV schemas and ordering match.
Environment fix: ASAN loader screenshot reproduced; aggregate oneAPI directory contained zero-length redirection entries. Copy real compiler/latest DLLs beside executables. CMake --fresh + full compiler path preserves Release/ASAN flags. Reduced-PATH launches pass.
MATLAB tests: 80/80 passed with explicit project absolute path.
C++ final tests: 7/7 Release and 7/7 AddressSanitizer, confirmed actual instrumentation flags.
Differential: box/refill/buddha_gravity strict numeric equality checks pass; all four cell CSV connectivity and counts identical. Long-shake Buddha trajectories differ .0291599 due initial dot rounding4.44e-16; all13400 same-state contacts within1e-10 (max7.916e-11). Report openly records numeric non-equivalence rather than claiming identical outputs.
Benchmark: 1 warmup +3 runs, both recorded replay and native RNG; final Buddha native1.527190s vs MATLAB29.408646s median=19.26x. Full22k requested C++74.862575s, accepts20108 capacity_reached; no full-scale MATLAB speedup claimed.
Physical validation: all20108 centres inside, pairgap>=-3.109e-15, wallgap>=-1.610e-15 using independent SciPy/NumPy calculations.
Delivery: README.md commands and BENCHMARK_REPORT.md evidence. Keep existing main-cpp workspace changes, no merge/push/commit.
