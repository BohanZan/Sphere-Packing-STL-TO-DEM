# Windows oneAPI C++ port

Goal: translate the current working MATLAB implementation, including its uncommitted fixes, to an independent Windows C++20 executable and library. Preserve ordered radii, batching (three empty initial attempts), continuous finite-feature contacts, batch-only compression/shaking, cumulative empty refill traversals, arbitrary gravity, relative tolerances, physical properties and all four CSV schemas.

Architecture: geometry.hpp/cpp own vectors, triangles, sparse surface/ray grids, exact parity, conservative occupancy, normals and finite-feature contact. packing.hpp/cpp own options, state, random source, incremental particle hash, DDA movement and packing. io.hpp/cpp and main.cpp own checked STL ingestion, CLI, CSV and reports. PowerShell initializes x64 Visual Studio and Intel oneAPI and builds with CMake/Ninja. RAII owns all allocations. Precise floating point is required; unsafe fast math is disabled.

STL decision: assess a small established header-only reader, pin source/license locally, wrap with input validation. Avoid a heavyweight rendering/geometry dependency. Preserve face order. Validate both ASCII and binary, truncation, nonfinite input and degenerate triangles.

Execution: use subagent-driven-development for independent geometry and packing units; controller implements IO/build/benchmarks. Existing user-created main-cpp checkout is the intended workspace; preserve pre-existing changes. No commit, merge or push is needed for local delivery.

## Tasks and verification
- [x] Geometry: tests first for finite triangle distances/contacts, parity seams, ray grid, occupied surface bands, occupancy and scale. Implement geometry translation against those and existing MATLAB tests.
- [x] Packing: tests first for swept sphere contact, zero-gap motion, batch isolation, gravity normalization, initial/refill rejection semantics, mass and inertia. Implement the state machine using geometry interfaces.
- [x] IO/build: add oneAPI x64 build script, checked STL reader, complete CLI option coverage, four CSV files and machine-readable report. Test invalid input and file schemas.
- [x] Differential validation: run MATLAB tests; export deterministic baseline fixtures and random draws for numerical replay; compare assembly/masses/inertia/report/grid CSVs and independently check physical gaps.
- [x] Benchmark: use matching mesh/radii/options, exclude application startup, warm up then repeated unprofiled elapsed measurements; report all timings and exact workload. Include real STL and a modest deterministic regression workload. Do not imply the full 22,000-sphere baseline was timed if it was not.
- [x] Safety/review: sanitizer build if supported on this Windows toolchain; memory/resource lifetime review; malformed-input tests and full CTest; independent final review.
- [x] Delivery: README with literal PowerShell build/run/benchmark commands and local measured timing report.

Randomness: native C++ normal distribution need not reproduce MATLAB's proprietary normal stream. Provide recorded MATLAB random-draw replay for rigorous like-for-like numerical validation; disclose seeded native stream differences. Timing runs do not include tracing overhead.

Progress and measured results will be recorded in CPP_PROGRESS.md and benchmark outputs.
