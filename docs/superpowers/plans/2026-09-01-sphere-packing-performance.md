# SpherePacking Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce SpherePacking runtime without changing the packing decision, random-number consumption, public API, or persisted CSV schema.

**Architecture:** Keep the current random-placement, relaxation, and refill algorithms. Replace only equivalent implementations on their hot paths: reject known colliding spheres before geometry work; cache static face data; evaluate ray and triangle queries in batches; classify provably safe occupancy cells before exact ray casting; then remove avoidable hash-query allocation. Every change is isolated behind the existing function interfaces and checked against a fixed-seed reference fixture.

**Tech Stack:** MATLAB R2025a, `matlab.unittest`, ASCII/binary STL read through the existing `stlRead`, `containers.Map`.

## Global Constraints

- Preserve the signatures and return layouts of `spawnSpheres`, `spCanPlace`, `spPointInside`, `spSphereHitsTriangles`, and every CSV writer.
- For the same `model`, `radii`, `maxAttempts`, `buffer`, `options`, and `options.randomSeed`, preserve every placement decision and make no additional calls to `rand`, `randn`, or `rng`.
- Preserve the current strict collision rule: `distance < radius - context.tolerance` for mesh tests and `distance < r1 + r2 - context.tolerance` for sphere tests.
- Preserve `spPointInside`'s existing tolerance, duplicate-height removal, and odd/even ray-parity rule. The occupancy grid in Task 5 is a conservative triage layer only: every `margin` cell must call the exact parity implementation, and a cell may be labelled `inside` or `outside` only when its closed volume is separated from every STL triangle by the rasterised margin halo.
- Do not use an approximate signed-distance field, BVH, GPU code, MEX file, or parallel placement in this plan.
- Preserve the four CSV schemas and the `report` fields. Benchmark-only output must go to a temporary directory and be deleted by test cleanup.
- Use `AbsTol = 1e-12` for existing scale-dependent regression comparisons, except when the existing test already specifies a different tolerance.
- Compare wall-clock runtime with `tic`/`toc`; use `profile` only to locate residual hot paths because the full profile records tens of millions of calls.

---

## Profiler-derived order of work

The supplied 937.921 s profile concentrates 98% of self time in five leaves: `verticalTriangleHit` (309.978 s), `spPointInside` (249.264 s), `spSphereHitsTriangles` (173.907 s), `spPointTriangleDistance` (94.835 s), and `spHashNeighbours` (89.311 s). `spReindex` (1.649 s), `spBuildContext` (0.555 s), and `spHashInsert` (1.231 s) are explicitly deferred.

The counts also expose an immediately safe fast path: `spHashNeighbours` ran 124,935 times and `spSphereHitsTriangles` ran 36,366 times. Since `spCanPlace` calls the sphere-neighbour query once and the triangle-neighbour query only after sphere clearance, about `124935 - 36366 = 88569` candidates reached the sphere check and `52203` of them were rejected for sphere overlap. Those candidates currently call `spPointInside` first. Moving the sphere test ahead of the ray query therefore avoids an expensive geometry query for roughly 59% of sphere-checked candidates without changing a Boolean result.

## File structure

| Path | Responsibility |
|---|---|
| `tests/helpers/spTestCubeMesh.m` | Shared closed cube mesh for performance-regression tests. |
| `tests/fixtures/seed53_cube_reference.mat` | Fixed-seed pre-optimization reference outputs; contains no output paths or timings. |
| `tests/benchmarkSpherePacking.m` | Three-run, fixed-seed wall-clock benchmark with automatic temporary-output cleanup. |
| `tests/testSpawnSpheres.m` | Existing physical invariants plus all new equivalence tests. |
| `tests/helpers/spReferencePointInside.m` | Scalar pre-vectorization oracle, copied verbatim from the current point-in-solid algorithm. |
| `tests/helpers/spReferenceSphereHitsTriangles.m` | Scalar pre-vectorization oracle for finite triangle collision. |
| `tests/helpers/spTestConcavePrism.m` | Watertight concave L-prism for occupancy-classification tests. |
| `tests/helpers/spTestHollowCubeMesh.m` | Watertight outer cube plus reversed-winding inner cube for cavity tests. |
| `tests/profileOccupancyAcceleration.m` | Fixed-seed comparison of exact-ray call counts with occupancy enabled and disabled. |
| `spCanPlace.m` | Reorder independent rejection tests; later consume separated sphere/triangle neighbour queries. |
| `spBuildContext.m` | Cache static ray-query and triangle-query data, then cache inward normals after the spatial hashes exist. |
| `spExactPointInside.m` | New exact, batched ray-parity implementation with no occupancy shortcut. |
| `spPointInside.m` | Public wrapper that returns a conservative occupancy result or calls `spExactPointInside`. |
| `spBuildOccupancyGrid.m` | Rasterise conservative margin cells, label six-connected free-space components, and store the triage grid. |
| `spOccupancyCellIndex.m` | Map a point to the independent occupancy-grid cell using the current clamp convention. |
| `spSphereHitsTriangles.m` | Batch finite-triangle distance calculation and stop on the same strict threshold. |
| `spHashNeighbours.m` | Triangle-neighbour query only; retains deduplication. |
| `spSphereNeighbours.m` | New sphere-only neighbour query; no deduplication because each sphere occupies one cell. |
| `spRefill.m` | Use cached inward normals and skip active-face work when no refill pass can run. |
| `spEmptyState.m`, `spAddSphere.m`, `spReindex.m`, `spawnSpheres.m` | Deferred preallocation and incremental-indexing work after the measured geometry fixes land. |

### Task 1: Freeze correctness and wall-clock baselines

**Files:**
- Create: `tests/helpers/spTestCubeMesh.m`
- Create: `tests/fixtures/seed53_cube_reference.mat`
- Create: `tests/benchmarkSpherePacking.m`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Produces `spTestCubeMesh(sideLength)`, returning `struct('vertices', double(8,3), 'faces', double(12,3))` with outward-consistent triangular faces.
- Produces `benchmarkSpherePacking(runCount)`, returning `struct('seconds', double(1,runCount), 'medianSeconds', double, 'acceptedCount', double, 'report', struct)`.
- Creates the fixture variable `expected` with fields `assembly`, `masses`, `totalVolume`, `inertia`, `acceptedCount`, `unplacedCount`, `centreOfMass`, and `coordinateFrame`.

- [ ] **Step 1: Add the shared deterministic cube helper**

Create `tests/helpers/spTestCubeMesh.m` with the same 8 vertices and 12 consistently oriented faces currently used by the local `cubeMesh` function. Use this public test helper body:

```matlab
function mesh = spTestCubeMesh(sideLength)
v = [0 0 0; sideLength 0 0; sideLength sideLength 0; 0 sideLength 0; ...
     0 0 sideLength; sideLength 0 sideLength; sideLength sideLength sideLength; 0 sideLength sideLength];
f = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
     2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8];
mesh = struct('vertices', v, 'faces', f);
end
```

- [ ] **Step 2: Generate the pre-change reference fixture once**

Before changing implementation files, run the following MATLAB expression from the SpherePacking root. It intentionally uses a temporary nonempty output directory because an empty directory currently triggers output-directory inference for STL paths.

```matlab
addpath(pwd); addpath('tests/helpers');
out = tempname; cleanup = onCleanup(@() rmdir(out, 's'));
options = struct('randomSeed', 53, 'outputDirectory', out, ...
    'outputPrefix', 'reference', 'coordinateFrame', 'world', ...
    'maxCompressionSweeps', 8, 'shakeSweeps', 1, 'maxRefillPasses', 1);
[assembly, masses, totalVolume, inertia, report] = ...
    spawnSpheres(spTestCubeMesh(20), [0.5; 0.75; 1.0; 1.25], 300, 0.01, options);
expected = struct('assembly', assembly, 'masses', masses, ...
    'totalVolume', totalVolume, 'inertia', inertia, ...
    'acceptedCount', report.acceptedCount, 'unplacedCount', report.unplacedCount, ...
    'centreOfMass', report.centreOfMass, 'coordinateFrame', report.coordinateFrame);
save('tests/fixtures/seed53_cube_reference.mat', 'expected');
clear cleanup
```

- [ ] **Step 3: Add a fixed-seed output-regression test**

Add this test to `tests/testSpawnSpheres.m`; omit `report.outputFiles`, which contains a temporary directory name rather than model semantics.

```matlab
function testFixedSeedPackingMatchesReferenceFixture(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
loaded = load(fullfile(root, 'tests', 'fixtures', 'seed53_cube_reference.mat'), 'expected');
out = tempname; cleanup = onCleanup(@() rmdir(out, 's'));
options = struct('randomSeed', 53, 'outputDirectory', out, ...
    'outputPrefix', 'reference', 'coordinateFrame', 'world', ...
    'maxCompressionSweeps', 8, 'shakeSweeps', 1, 'maxRefillPasses', 1);
[assembly, masses, totalVolume, inertia, report] = ...
    spawnSpheres(spTestCubeMesh(20), [0.5; 0.75; 1.0; 1.25], 300, 0.01, options);
verifyEqual(testCase, assembly, loaded.expected.assembly, 'AbsTol', 1e-12);
verifyEqual(testCase, masses, loaded.expected.masses, 'AbsTol', 1e-12);
verifyEqual(testCase, totalVolume, loaded.expected.totalVolume, 'AbsTol', 1e-12);
verifyEqual(testCase, inertia, loaded.expected.inertia, 'AbsTol', 1e-12);
verifyEqual(testCase, report.acceptedCount, loaded.expected.acceptedCount);
verifyEqual(testCase, report.unplacedCount, loaded.expected.unplacedCount);
verifyEqual(testCase, report.centreOfMass, loaded.expected.centreOfMass, 'AbsTol', 1e-12);
verifyEqual(testCase, report.coordinateFrame, loaded.expected.coordinateFrame);
clear cleanup
end
```

- [ ] **Step 4: Add the three-run wall-clock benchmark**

Create `tests/benchmarkSpherePacking.m` with this benchmark body. It writes outside the input directory and deletes only its own explicitly created temporary directory.

```matlab
function result = benchmarkSpherePacking(runCount)
if nargin < 1, runCount = 3; end
root = fileparts(fileparts(mfilename('fullpath')));
model = fullfile(root, 'inputs', 'ironParticle', 'ironParticle.stl');
radii = repmat(3.0e-6, 100, 1);
result.seconds = zeros(1, runCount);
for runId = 1:runCount
    out = tempname; cleanup = onCleanup(@() rmdir(out, 's'));
    options = struct('randomSeed', 53, 'outputDirectory', out, ...
        'outputPrefix', 'benchmark', 'coordinateFrame', 'world', ...
        'maxCompressionSweeps', 100, 'shakeSweeps', 2, 'maxRefillPasses', 3);
    started = tic;
    [~, ~, ~, ~, report] = spawnSpheres(model, radii, 100, 0, options);
    result.seconds(runId) = toc(started);
    result.acceptedCount(runId) = report.acceptedCount;
    clear cleanup
end
result.medianSeconds = median(result.seconds);
end
```

- [ ] **Step 5: Run the regression suite and establish the unprofiled baseline**

Run:

```powershell
& 'D:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(pwd); addpath('tests/helpers'); results=runtests('tests'); assertSuccess(results); benchmark=benchmarkSpherePacking(3); disp(benchmark)"
```

Expected: all tests pass; record the three seconds values, their median, and the three identical accepted counts in the implementation PR description.

- [ ] **Step 6: Commit the baseline guard**

```powershell
git add tests/helpers/spTestCubeMesh.m tests/fixtures/seed53_cube_reference.mat tests/benchmarkSpherePacking.m tests/testSpawnSpheres.m
git commit -m "test: add sphere packing performance regression baseline"
```

### Task 2: Reject known sphere overlaps before ray casting

**Files:**
- Modify: `spCanPlace.m:7-25`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Consumes the existing `spHashNeighbours(context, state.sphereCells, idx)` and unchanged `spPointInside(context, centre)`.
- Produces the same scalar logical result from `spCanPlace(context, state, centre, radius, ignoreId)`.

- [ ] **Step 1: Add direct equivalence cases for `spCanPlace`**

Add a test that creates a cube context, inserts one sphere at `[10 10 10]`, and checks one overlapping interior candidate, one non-overlapping interior candidate, and one exterior candidate.

```matlab
function testCanPlacePreservesCollisionAndGeometryDecisions(testCase)
context = spBuildContext(spTestCubeMesh(20), 1.0, 0.01, 1e-9);
state = spEmptyState(context); state = spAddSphere(context, state, [10 10 10], 1.0);
verifyFalse(testCase, spCanPlace(context, state, [10.5 10 10], 1.0));
verifyTrue(testCase, spCanPlace(context, state, [15 15 15], 1.0));
verifyFalse(testCase, spCanPlace(context, state, [21 10 10], 1.0));
end
```

- [ ] **Step 2: Run the new test against the current ordering**

Run:

```powershell
& 'D:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(pwd); addpath('tests/helpers'); results=runtests('tests/testSpawnSpheres.m'); assertSuccess(results)"
```

Expected: PASS before the refactor; this establishes that the test checks semantic outcomes rather than implementation order.

- [ ] **Step 3: Reorder only independent rejection checks**

In `spCanPlace.m`, retain the bounds test first. Move the current sphere-neighbour retrieval and overlap loop immediately after it; leave the point-inside test and triangle-collision test unchanged after that loop.

```matlab
if any(centre-radius<context.lower) || any(centre+radius>context.upper), return; end
idx = spCellIndex(context, centre);
ids = spHashNeighbours(context, state.sphereCells, idx);
for id = ids
    if id == ignoreId, continue; end
    if norm(centre-state.centres(id,:)) < radius+state.radii(id)-context.tolerance, return; end
end
if ~spPointInside(context, centre), return; end
triIds = spHashNeighbours(context, context.triangleCells, idx);
if spSphereHitsTriangles(context, centre, radius, triIds), return; end
yes = true;
```

- [ ] **Step 4: Verify reference equality and measure the isolated gain**

Run the fixture test and `benchmarkSpherePacking(3)`. Require all three accepted counts to match Task 1 and the fixture test to pass. Require the median unprofiled time to be lower than Task 1; record both medians.

- [ ] **Step 5: Commit the ordering fix**

```powershell
git add spCanPlace.m tests/testSpawnSpheres.m
git commit -m "perf: reject overlapping sphere candidates before ray tests"
```

### Task 3: Cache inward face normals and avoid disabled-refill work

**Files:**
- Modify: `spBuildContext.m:62-65`
- Modify: `spRefill.m:3-56`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Adds `context.faceCentres` as `double(F,3)` and `context.inwardNormals` as `double(F,3)`, where `F = size(context.faces,1)`.
- Keeps `spRefill(context, state, radii, nextRadius, options)` unchanged.

- [ ] **Step 1: Add tests for cached normal orientation and disabled refill**

Add this cached-normal test using the cube. It checks the same inward-probe rule currently used by `inwardNormal`.

```matlab
function testCachedInwardNormalsPointIntoTheClosedMesh(testCase)
context = spBuildContext(spTestCubeMesh(20), 1.0, 0.01, 1e-9);
probeDistance = max(context.tolerance*100, 1e-8*context.cellSize);
for id = 1:size(context.faces, 1)
    probe = context.faceCentres(id,:) + probeDistance*context.inwardNormals(id,:);
    verifyTrue(testCase, spPointInside(context, probe));
end
end
```

Also add a refill-disabled case with `maxRefillPasses = 0`, a state containing one sphere, and `nextRadius = 1`; verify `state.centres`, `state.radii`, `state.count`, and `nextRadius` are unchanged.

- [ ] **Step 2: Build face centres and normals after the ray hash exists**

After `xyCells` is built, construct a temporary context with all fields needed by `spPointInside`. For each face, compute exactly the old normal/probe decision once; append `faceCentres` and `inwardNormals` to the final context.

```matlab
faceVertices = reshape(vertices(faces.', :), 3, size(faces,1), 3);
faceCentres = squeeze(mean(faceVertices, 1));
rawNormals = cross(vertices(faces(:,2),:) - vertices(faces(:,1),:), ...
                   vertices(faces(:,3),:) - vertices(faces(:,1),:), 2);
rawNormals = rawNormals ./ vecnorm(rawNormals, 2, 2);
```

For each row of `rawNormals`, use the existing probe rule and `spPointInside` to select either `rawNormals(id,:)` or its negation. Do not use signed-volume orientation in this task because it could alter behaviour for inconsistent face winding.

- [ ] **Step 3: Consume cached normals in refill**

At the start of `spRefill`, return before computing `active` when `options.maxRefillPasses < 1` or `nextRadius > numel(radii)`. Replace active-face selection and sampled-face normal calculation with cached values.

```matlab
if options.maxRefillPasses < 1 || nextRadius > numel(radii), return; end
active = find(context.inwardNormals * options.gravity.' > 0).';
if isempty(active), active = 1:size(context.faces,1); end
% Inside the candidate loop:
n = context.inwardNormals(id,:);
```

- [ ] **Step 4: Run all correctness checks and the benchmark**

Run the full suite, fixed-seed fixture, and three-run benchmark. Require fixture equality and a lower median time than the Task 2 baseline. Inspect a fresh profile: `spRefill>inwardNormal` must no longer appear as a repeated hot function.

- [ ] **Step 5: Commit the cached geometry data**

```powershell
git add spBuildContext.m spRefill.m tests/testSpawnSpheres.m
git commit -m "perf: cache inward face normals for refill"
```

### Task 4: Batch vertical-ray intersection without changing parity semantics

**Files:**
- Create: `tests/helpers/spReferencePointInside.m`
- Create: `spExactPointInside.m`
- Modify: `spBuildContext.m`
- Modify: `spPointInside.m`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Adds `context.ray.a`, `context.ray.ab`, `context.ray.ac`, `context.ray.determinant`, `context.ray.inverseDeterminant`, and `context.ray.zDelta` as face-indexed numeric arrays.
- Produces `inside = spExactPointInside(context, point)`, the scalar logical exact-parity query used by later occupancy construction.
- Keeps `inside = spPointInside(context, point)` scalar and logical; until Task 5 it delegates directly to `spExactPointInside`.

- [ ] **Step 1: Preserve the scalar algorithm as a test oracle**

Create `tests/helpers/spReferencePointInside.m` by copying the pre-vectorization bodies of `spPointInside` and `verticalTriangleHit`, renaming only the outer function. It must retain `sort`, `diff(heights) > context.tolerance`, and the final odd/even expression exactly.

```matlab
function inside = spReferencePointInside(context, point)
% Copy the current scalar spPointInside body here before changing it.
% Keep its local verticalTriangleHit helper unchanged.
end
```

- [ ] **Step 2: Add scalar-versus-batched point tests**

Use all cube-cell centres plus 100 deterministic random points from the cube bounding box. Test both ordinary-scale and `1e-6` scale meshes.

```matlab
rng(71, 'twister');
points = [10 10 10; 0 0 0; 20 20 20; 20.1 10 10; 20*rand(100,3)];
for k = 1:size(points,1)
    verifyEqual(testCase, spPointInside(context, points(k,:)), ...
        spReferencePointInside(context, points(k,:)));
end
```

- [ ] **Step 3: Precompute ray coefficients in `spBuildContext`**

Store arrays indexed by face ID so `spPointInside` never reconstructs triangle rows inside its hot loop.

```matlab
a = vertices(faces(:,1),:);
b = vertices(faces(:,2),:);
c = vertices(faces(:,3),:);
ab = b(:,1:2) - a(:,1:2);
ac = c(:,1:2) - a(:,1:2);
determinant = ab(:,1).*ac(:,2) - ab(:,2).*ac(:,1);
context.ray = struct('a', a, 'ab', ab, 'ac', ac, ...
    'determinant', determinant, 'inverseDeterminant', 1./determinant, ...
    'zDelta', [b(:,3)-a(:,3), c(:,3)-a(:,3)]);
```

- [ ] **Step 4: Move the exact per-face loop into a batched exact-query function**

For `ids = context.xyCells(key)`, compute `alpha`, `beta`, `hit`, and `zHit` for every selected face at once. Mask degenerate projections before using the inverse determinant, apply the same z tolerance, sort the remaining heights, and retain the current duplicate-removal expression.

```matlab
ray = context.ray; ids = ids(:);
ap = point(1:2) - ray.a(ids,1:2);
alpha = (ap(:,1).*ray.ac(ids,2) - ap(:,2).*ray.ac(ids,1)) .* ray.inverseDeterminant(ids);
beta = (ray.ab(ids,1).*ap(:,2) - ray.ab(ids,2).*ap(:,1)) .* ray.inverseDeterminant(ids);
valid = abs(ray.determinant(ids)) > context.tolerance^2 & ...
    alpha >= -context.tolerance & beta >= -context.tolerance & ...
    alpha + beta <= 1 + context.tolerance;
zHit = ray.a(ids,3) + alpha.*ray.zDelta(ids,1) + beta.*ray.zDelta(ids,2);
heights = sort(zHit(valid & zHit < point(3)-context.tolerance));
inside = ~isempty(heights) && mod(1 + sum(diff(heights) > context.tolerance), 2) == 1;
```

Put this body in `spExactPointInside.m`. Replace `spPointInside.m` with the unchanged public delegation below; Task 5 will add the occupancy shortcut around this exact implementation.

```matlab
function inside = spPointInside(context, point)
inside = spExactPointInside(context, point);
end
```

- [ ] **Step 5: Verify geometry equivalence before measuring speed**

Run all scalar-oracle tests, the fixture test, then the full suite. Only after all pass, run `benchmarkSpherePacking(3)` and inspect a profile. Require the scalar `verticalTriangleHit` leaf to be absent from the hot-path table.

- [ ] **Step 6: Commit the batched exact-ray query**

```powershell
git add spBuildContext.m spExactPointInside.m spPointInside.m tests/helpers/spReferencePointInside.m tests/testSpawnSpheres.m
git commit -m "perf: batch ray intersections in point-in-solid tests"
```

### Task 5: Add a conservative occupancy grid before exact ray casting

**Files:**
- Create: `spBuildOccupancyGrid.m`
- Create: `spOccupancyCellIndex.m`
- Create: `tests/helpers/spTestConcavePrism.m`
- Create: `tests/helpers/spTestHollowCubeMesh.m`
- Create: `tests/profileOccupancyAcceleration.m`
- Modify: `spawnSpheres.m:11,80-118`
- Modify: `spBuildContext.m:1-7,62-65`
- Modify: `spPointInside.m`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Extends the internal call to `spBuildContext(model, maxRadius, buffer, tolerance, occupancyOptions)`; its fifth argument is optional, so every existing four-argument test call remains valid and builds no occupancy shortcut.
- Adds `options.occupancyAcceleration` (`true`), `options.occupancyCellSize` (`[]`), and `options.occupancyMaxCells` (`2e6`) in `spDefaultOptions`.
- Adds `context.occupancy = struct('enabled', logical, 'lower', double(1,3), 'cellSize', double, 'cellCount', double(1,3), 'labels', uint8(Nx,Ny,Nz))`.
- Uses label values `0 = outside`, `1 = inside`, `2 = margin`, `3 = unlabelled`. `unlabelled` must never remain after construction.
- `spPointInside(context, point)` remains scalar/logical and is still exact for every `margin` cell because it delegates to `spExactPointInside` there.

- [ ] **Step 1: Add convex, concave, and cavity mesh helpers plus failing classification tests**

Create `tests/helpers/spTestConcavePrism.m` as this watertight L-prism. Its top and bottom use a four-triangle ear decomposition; every side quad is split consistently.

```matlab
function mesh = spTestConcavePrism()
v = [0 0 0; 3 0 0; 3 1 0; 1 1 0; 1 3 0; 0 3 0; ...
     0 0 3; 3 0 3; 3 1 3; 1 1 3; 1 3 3; 0 3 3];
bottom = [1 4 2; 2 4 3; 1 6 4; 4 6 5];
top = [7 8 10; 8 9 10; 7 10 12; 10 11 12];
sides = [1 2 8; 1 8 7; 2 3 9; 2 9 8; 3 4 10; 3 10 9; ...
         4 5 11; 4 11 10; 5 6 12; 5 12 11; 6 1 7; 6 7 12];
mesh = struct('vertices', v, 'faces', [bottom; top; sides]);
end
```

Create `tests/helpers/spTestHollowCubeMesh.m` by adding a reversed-winding inner cube to the existing outward-oriented cube.

```matlab
function mesh = spTestHollowCubeMesh()
outer = spTestCubeMesh(20);
innerVertices = [5 5 5; 15 5 5; 15 15 5; 5 15 5; ...
                 5 5 15; 15 5 15; 15 15 15; 5 15 15];
innerFaces = outer.faces(:, [1 3 2]) + 8;
mesh = struct('vertices', [outer.vertices; innerVertices], ...
    'faces', [outer.faces; innerFaces]);
end
```

Add tests for the points below with occupancy explicitly enabled. They must all agree with `spExactPointInside`; the listed truth value is an independent geometric assertion.

```matlab
function testOccupancyHandlesConvexConcaveAndCavityMeshes(testCase)
cfg = struct('enabled', true, 'cellSize', 0.25, 'maxCells', 2e6);
convex = spBuildContext(spTestCubeMesh(20), 1.0, 0, 1e-9, cfg);
concave = spBuildContext(spTestConcavePrism(), 0.25, 0, 1e-9, cfg);
hollow = spBuildContext(spTestHollowCubeMesh(), 0.5, 0, 1e-9, cfg);
cases = {convex, [10 10 10], true; convex, [21 10 10], false; ...
         concave, [0.5 2.0 1.5], true; concave, [2.0 2.0 1.5], false; ...
         hollow, [2 10 10], true; hollow, [10 10 10], false; hollow, [21 10 10], false};
for k = 1:size(cases,1)
    context = cases{k,1}; point = cases{k,2}; expected = cases{k,3};
    verifyEqual(testCase, spExactPointInside(context, point), expected);
    verifyEqual(testCase, spPointInside(context, point), expected);
end
end
```

- [ ] **Step 2: Thread a bounded, optional occupancy configuration through context construction**

In `spDefaultOptions`, add the three scalar defaults and validate them. In `spawnSpheres`, pass a struct with exact field names to the optional fifth `spBuildContext` argument.

```matlab
defaults = struct(..., 'occupancyAcceleration', true, ...
    'occupancyCellSize', [], 'occupancyMaxCells', 2e6);
occupancyOptions = struct('enabled', logical(options.occupancyAcceleration), ...
    'cellSize', options.occupancyCellSize, 'maxCells', options.occupancyMaxCells);
context = spBuildContext(model, max(radii), options.buffer, ...
    options.tolerance, occupancyOptions);
```

For the four-argument `spBuildContext` form, set `occupancyOptions = struct('enabled', false, 'cellSize', [], 'maxCells', 2e6)`. This preserves the exact current behaviour for direct helper and legacy test callers.

- [ ] **Step 3: Build a conservative independent grid**

Create `spBuildOccupancyGrid(context, maxRadius, occupancyOptions)`. Select the requested cell width as `max(maxRadius/2, 32*context.tolerance)` when no explicit width was supplied. If `prod(cellCount) > occupancyOptions.maxCells`, multiply the width by `ceil((prod(cellCount)/occupancyOptions.maxCells)^(1/3))`, recompute counts, and continue. A coarser grid may produce more `margin` cells but must never produce a nonconservative label.

Rasterise every face's **expanded axis-aligned bounding box** into `margin`. The expansion is one full occupancy cell plus `max(context.tolerance, 8*eps(max(abs(context.vertices(:)))))`; the AABB approximation intentionally marks extra cells rather than risking a surface hole.

```matlab
labels = 3 * ones(cellCount, 'uint8');
margin = false(cellCount);
halo = max(context.tolerance, 8*eps(max(abs(context.vertices(:)))));
for id = 1:size(context.faces,1)
    triangle = context.vertices(context.faces(id,:),:);
    lo = min(triangle, [], 1) - occupancyCellSize - halo;
    hi = max(triangle, [], 1) + occupancyCellSize + halo;
    loIndex = spOccupancyCellIndex(occupancy, lo);
    hiIndex = spOccupancyCellIndex(occupancy, hi);
    margin(loIndex(1):hiIndex(1), loIndex(2):hiIndex(2), loIndex(3):hiIndex(3)) = true;
end
labels(margin) = uint8(2);
```

`spOccupancyCellIndex(occupancy, point)` must use the existing lower-bound, `floor`, one-based indexing, and clamp convention:

```matlab
index = floor((point - occupancy.lower) ./ occupancy.cellSize) + 1;
index = min(max(index, 1), occupancy.cellCount);
```

- [ ] **Step 4: Label only non-margin components, using exact parity once per component**

Flood-fill all six-connected, non-margin boundary components as `outside`. Then flood-fill every remaining six-connected non-margin component; use its first cell centre as a representative and call `spExactPointInside` once to label the whole component `inside` or `outside`. Store no label `3` when returning.

```matlab
% For each unlabelled component returned by a six-neighbour FIFO flood fill:
[ix, iy, iz] = ind2sub(occupancy.cellCount, component(1));
representative = occupancy.lower + ([ix iy iz] - 0.5) * occupancy.cellSize;
if spExactPointInside(context, representative)
    labels(component) = uint8(1);
else
    labels(component) = uint8(0);
end
```

Use a `uint32` FIFO and six axis-neighbours only. The one-cell AABB expansion plus six-connectivity is deliberate: it over-classifies uncertainty as `margin` instead of allowing diagonal surface leaks to create a false `inside` cell.

- [ ] **Step 5: Add the public triage wrapper with exact margin fallback**

Replace the temporary delegation from Task 4 with this wrapper. It is the only production consumer of occupancy labels; `spBuildOccupancyGrid` and every validation test call `spExactPointInside` directly.

```matlab
function inside = spPointInside(context, point)
if isfield(context, 'occupancy') && context.occupancy.enabled
    index = spOccupancyCellIndex(context.occupancy, point);
    label = context.occupancy.labels(index(1), index(2), index(3));
    if label == 0, inside = false; return; end
    if label == 1, inside = true; return; end
end
inside = spExactPointInside(context, point);
end
```

Do not change `spSphereHitsTriangles`: a classified-inside centre still requires the existing local triangle-distance test to prove the entire sphere lies inside the STL surface.

- [ ] **Step 6: Add conservative-cell and full-packing equivalence tests**

For every non-margin cell in each of the three helper contexts, compare `spExactPointInside` at the cell centre and at offsets `[-0.49, 0, 0.49] * cellSize` along every axis with the label. Add a fixed-seed packing test that runs the same cube mesh once with `occupancyAcceleration = false` and once with it `true`; require equal assembly, masses, volume, inertia, accepted count, and sphere-volume fraction.

```matlab
fractionExact = exactReport.sphereAssemblyVolume / exactReport.stlVolume;
fractionFast = fastReport.sphereAssemblyVolume / fastReport.stlVolume;
verifyEqual(testCase, fastAssembly, exactAssembly, 'AbsTol', 1e-12);
verifyEqual(testCase, fractionFast, fractionExact, 'AbsTol', 1e-12);
verifyEqual(testCase, fastReport.acceptedCount, exactReport.acceptedCount);
```

- [ ] **Step 7: Measure exact-ray avoidance and benchmark only after equivalence passes**

Create `tests/profileOccupancyAcceleration.m`. For the iron-particle benchmark input, run once with `occupancyAcceleration = false` and once with it `true`, each with seed 53 and profiling enabled. Extract the `NumCalls` entry whose `FunctionName` ends in `spExactPointInside`; also return accepted count and `sphereAssemblyVolume/stlVolume` for both runs.

```matlab
entry = info.FunctionTable(endsWith({info.FunctionTable.FunctionName}, 'spExactPointInside'));
rayCalls = sum([entry.NumCalls]);
```

Require identical accepted count and sphere-volume fraction, fewer exact-ray calls with occupancy enabled, and a lower median from `benchmarkSpherePacking(3)`. If the conservative margin fills nearly all cells and ray calls do not fall for a particular STL/radius combination, retain the implemented feature and record its measured bypass ratio; do not weaken its conservative classification rules or alter packing logic to manufacture a speedup.

- [ ] **Step 8: Commit occupancy acceleration separately**

```powershell
git add spawnSpheres.m spBuildContext.m spBuildOccupancyGrid.m spOccupancyCellIndex.m spPointInside.m tests/helpers/spTestConcavePrism.m tests/helpers/spTestHollowCubeMesh.m tests/profileOccupancyAcceleration.m tests/testSpawnSpheres.m
git commit -m "perf: bypass exact rays in conservatively classified cells"
```

### Task 6: Batch finite-triangle distance tests

**Files:**
- Create: `tests/helpers/spReferenceSphereHitsTriangles.m`
- Modify: `spBuildContext.m`
- Modify: `spSphereHitsTriangles.m`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Adds `context.triangles.a`, `context.triangles.b`, and `context.triangles.c`, each `double(F,3)`.
- Keeps `intersects = spSphereHitsTriangles(context, centre, radius, triangleIds)` scalar and logical.

- [ ] **Step 1: Preserve scalar mesh-hit logic as an oracle**

Create `tests/helpers/spReferenceSphereHitsTriangles.m` from the current implementation and its scalar `spPointTriangleDistance` calculation. The helper must use the same strict test:

```matlab
if spPointTriangleDistance(centre, triangle) < radius - context.tolerance
    intersects = true;
    return;
end
```

- [ ] **Step 2: Add exact boundary and random oracle cases**

For a cube context, compare the optimized and reference result for the six face-centre directions at radii `0.99`, `1.00`, and `1.01`, plus 100 deterministic random centres/radii and the complete face-ID list.

```matlab
ids = 1:size(context.faces,1);
rng(73, 'twister');
for k = 1:100
    centre = 20*rand(1,3); radius = 0.05 + 2.0*rand;
    verifyEqual(testCase, spSphereHitsTriangles(context, centre, radius, ids), ...
        spReferenceSphereHitsTriangles(context, centre, radius, ids));
end
```

- [ ] **Step 3: Cache triangle vertices by face ID**

Build the numeric triangle cache once in `spBuildContext`.

```matlab
context.triangles = struct('a', vertices(faces(:,1),:), ...
    'b', vertices(faces(:,2),:), 'c', vertices(faces(:,3),:));
```

- [ ] **Step 4: Implement a batched squared-distance kernel**

Replace the outer `for id = triangleIds` loop with vector operations that implement the same seven Voronoi-region cases as `spPointTriangleDistance`. Compare against `limitSquared = (radius-context.tolerance)^2` only after guarding `radius <= context.tolerance`.

```matlab
limit = radius - context.tolerance;
if isempty(triangleIds) || limit <= 0, intersects = false; return; end
tri = context.triangles; ids = triangleIds(:);
% Compute d1..d6, va/vb/vc, closest points, and d2 for all ids.
% Assign each row exactly once in the scalar function's branch order.
intersects = any(d2 < limit^2);
```

- [ ] **Step 5: Verify the oracle and baseline fixture**

Run the distance oracle, all unit tests, and the fixture test. Then run the benchmark. A profile must show neither millions of calls to `spPointTriangleDistance` nor a dominant self-time in `spSphereHitsTriangles`.

- [ ] **Step 6: Commit the batched surface distance path**

```powershell
git add spBuildContext.m spSphereHitsTriangles.m tests/helpers/spReferenceSphereHitsTriangles.m tests/testSpawnSpheres.m
git commit -m "perf: batch finite triangle collision checks"
```

### Task 7: Remove avoidable hash-neighbour allocation

**Files:**
- Create: `spSphereNeighbours.m`
- Modify: `spHashNeighbours.m`
- Modify: `spCanPlace.m`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- `spSphereNeighbours(context, hash, index)` returns a row vector of unique sphere IDs without calling `unique`.
- `spHashNeighbours(context, hash, index)` remains the triangle query and returns unique triangle IDs.

- [ ] **Step 1: Add a neighbour-equivalence test**

Build a cube context, insert spheres whose centres occupy the central cell and every valid adjacent cell, then compare sorted IDs from the new sphere query with the old generic query.

```matlab
actual = sort(spSphereNeighbours(context, state.sphereCells, [2 2 2]));
expected = sort(spHashNeighbours(context, state.sphereCells, [2 2 2]));
verifyEqual(testCase, actual, expected);
```

- [ ] **Step 2: Implement a preallocated cell-list collector for sphere IDs**

Use a 27-element cell array, append nonempty map values into its next slot, and concatenate once. Do not call `unique`, because `spAddSphere` indexes every accepted sphere in exactly one spatial cell.

```matlab
parts = cell(27,1); partCount = 0;
for ix = max(1,index(1)-1):min(context.cellCount(1),index(1)+1)
    for iy = max(1,index(2)-1):min(context.cellCount(2),index(2)+1)
        for iz = max(1,index(3)-1):min(context.cellCount(3),index(3)+1)
            key = sprintf('%d,%d,%d', ix, iy, iz);
            if isKey(hash, key)
                partCount = partCount + 1;
                parts{partCount} = hash(key);
            end
        end
    end
end
ids = [parts{1:partCount}];
```

- [ ] **Step 3: Vectorize the sphere overlap test**

Replace the sphere-ID loop in `spCanPlace` with squared distances, while keeping `ignoreId` semantics.

```matlab
ids = spSphereNeighbours(context, state.sphereCells, idx);
ids(ids == ignoreId) = [];
if ~isempty(ids)
    delta = state.centres(ids,:) - centre;
    overlap = sum(delta.^2, 2) < (radius + state.radii(ids) - context.tolerance).^2;
    if any(overlap), return; end
end
```

- [ ] **Step 4: Keep triangle deduplication but remove incremental concatenation**

Apply the same `parts` collector to `spHashNeighbours`, then keep exactly one `unique([parts{1:partCount}])` call. This retains correctness for triangles that belong to multiple overlapped cells.

- [ ] **Step 5: Verify and profile the narrowed hash target**

Run the neighbour-equivalence test, full suite, fixture, and benchmark. In a fresh profile, require lower self time for `spHashNeighbours`; `unique` may remain for triangle candidates.

- [ ] **Step 6: Commit hash-query allocation changes**

```powershell
git add spSphereNeighbours.m spHashNeighbours.m spCanPlace.m tests/testSpawnSpheres.m
git commit -m "perf: avoid duplicate sphere neighbour allocation"
```

### Task 8: Defer large-N memory improvements until geometry no longer dominates

**Files:**
- Modify: `spEmptyState.m`, `spAddSphere.m`, `spReindex.m`, `spawnSpheres.m`
- Modify: `tests/testSpawnSpheres.m`

**Interfaces:**
- Changes `spEmptyState(context)` to `spEmptyState(context, capacity)` where `capacity = numel(radii)`.
- Adds `state.cellIndices` as `double(capacity,3)` and keeps `state.count` as the count of valid rows.

- [ ] **Step 1: Add a capacity and slicing regression test**

Run a request for four spheres that accepts fewer than four and verify no preallocated tail appears in the returned assembly or CSV table.

```matlab
[assembly, ~, ~, ~, report] = spawnSpheres(spTestCubeMesh(2.10), [1;1;1;1], 100, 0.01);
verifyEqual(testCase, size(assembly,2), report.acceptedCount);
verifyEqual(testCase, size(assembly,1), 4);
verifyLessThan(testCase, report.acceptedCount, 4);
```

- [ ] **Step 2: Preallocate state arrays and slice at every public boundary**

Allocate `centres`, `radii`, and `cellIndices` to `capacity` rows. Pass `numel(radii)` from `spawnSpheres`. Before mass properties, assembly creation, CSV output, and any loops that consume all spheres, use `1:state.count` explicitly.

```matlab
state = struct('centres', zeros(capacity,3), 'radii', zeros(capacity,1), ...
    'cellIndices', zeros(capacity,3), 'count', 0, ...
    'sphereCells', containers.Map('KeyType','char','ValueType','any'), ...
    'nextProgressPercent', 25);
```

- [ ] **Step 3: Make reindex incremental only after preallocation passes**

For each moved sphere, compare its new `spCellIndex` with `state.cellIndices(id,:)`. If equal, do not touch the hash. If different, remove that ID from the old bucket, insert it into the new bucket, and update only `state.cellIndices(id,:)`. Keep a full `spReindex` helper for test setup and assert that its rebuilt map agrees with the incremental map after every relaxation sweep.

- [ ] **Step 4: Verify this deferred work does not trade correctness for a tiny measured gain**

Run the complete suite, fixture, and benchmark. Retain the change only if all correctness checks pass and the median is no slower than after Task 7. If it is slower, keep preallocation but revert only the incremental-map substep in its own commit.

- [ ] **Step 5: Commit the large-N hardening separately**

```powershell
git add spEmptyState.m spAddSphere.m spReindex.m spawnSpheres.m tests/testSpawnSpheres.m
git commit -m "perf: preallocate packing state for larger assemblies"
```

### Task 9: Final evidence and documentation

**Files:**
- Modify: `README.md`
- Modify: `docs/superpowers/plans/2026-09-01-sphere-packing-performance.md`

**Interfaces:**
- Documents the unchanged input/output contract and the exact benchmark command.

- [ ] **Step 1: Run a clean three-run final benchmark**

Run `benchmarkSpherePacking(3)` with MATLAB profiling disabled. Report the three raw times, median, fixed random seed, accepted counts, MATLAB version, and machine CPU model in the README performance section.

- [ ] **Step 2: Run the complete regression suite**

Run:

```powershell
& 'D:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(pwd); addpath('tests/helpers'); results=runtests('tests'); assertSuccess(results)"
```

Expected: every test passes, including CSV schema, physical non-overlap, scale-independent inside checks, fixed-seed fixture equality, ray oracle equivalence, triangle oracle equivalence, and neighbour-query equivalence.

- [ ] **Step 3: Run one final profile only to find residual work**

Profile a short, fixed-seed case first. If a flame graph is required, use a reduced sphere count that completes with fewer than five million recorded calls instead of raising `historysize` on the full 100-sphere job. Confirm that the old leaf functions no longer dominate the table.

- [ ] **Step 4: Document the result without claiming algorithmic changes**

Add a README note stating that the optimizer preserves the existing placement/relaxation/refill algorithm and accelerates only its equivalent geometry and hash-query implementation. Include the before/after median wall-clock times and the command used to obtain them.

- [ ] **Step 5: Commit the evidence**

```powershell
git add README.md docs/superpowers/plans/2026-09-01-sphere-packing-performance.md
git commit -m "docs: record sphere packing performance verification"
```

## Self-review

**Spec coverage:** Tasks 2–7 address every profiler-confirmed hotspot. Task 5 makes conservative occupancy classification a required part of the implementation, not an optional future redesign. Tasks 1 and 9 enforce unchanged semantics and measure unprofiled wall time. Task 8 retains dynamic-array improvements but correctly postpones them because their measured cost is below 0.2% of the supplied profile.

**Placeholder scan:** The plan names every modified/created file, function interface, test fixture, and benchmark command. The only intentionally deferred element is an approximate/BVH/parallel redesign, which is excluded by the global constraints because it would expand the behavioural scope.

**Type consistency:** `context.ray`, `context.triangles`, `context.faceCentres`, and `context.inwardNormals` are created in `spBuildContext` before their consumers. `spSphereNeighbours` is introduced before `spCanPlace` consumes it. `spEmptyState` capacity is supplied by `spawnSpheres` before state allocation is changed.

## Execution handoff

Plan complete and saved to `docs/superpowers/plans/2026-09-01-sphere-packing-performance.md`.

1. **Subagent-driven execution (recommended):** implement one task at a time, with an independent review after each task.
2. **Inline execution:** implement Tasks 1–5 in order, checkpoint against the benchmark, then continue only if the fixed-seed regression remains exact.
