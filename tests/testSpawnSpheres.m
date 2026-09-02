function tests = testSpawnSpheres
% Regression tests for the random non-overlapping sphere insertion core.
%Return local function tests so every physical invariant is run independently.
tests = functiontests(localfunctions);
end

function testFixedSeedPackingMatchesReferenceFixture(testCase)
%FIXEDSEEDPACKINGMATCHESREFERENCEFIXTURE Guard future optimisations exactly.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
fixture = fullfile(root, 'tests', 'fixtures', 'seed53_cube_reference.mat');
if ~isfile(fixture)
    verifyTrue(testCase, false, ...
        'The pre-optimisation fixed-seed fixture must exist before this test runs.');
    return;
end

loaded = load(fixture, 'expected');
out = tempname;
cleanup = onCleanup(@() removeOutputDirectory(out));
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

function testBenchmarkReturnsOneDeterministicSmokeSample(testCase)
%BENCHMARKRETURNSONEDETERMINISTICSMOKESAMPLE Exercise the benchmark harness.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
if isempty(which('benchmarkSpherePacking'))
    verifyNotEmpty(testCase, which('benchmarkSpherePacking'), ...
        'Task 1 must provide the benchmarkSpherePacking smoke-test harness.');
    return;
end

config = struct('model', spTestCubeMesh(20), 'radii', [0.5; 0.75], ...
    'maxAttempts', 300, 'buffer', 0.01, ...
    'options', struct('randomSeed', 59, 'maxCompressionSweeps', 3, ...
    'shakeSweeps', 0, 'maxRefillPasses', 0, 'coordinateFrame', 'world'));
result = benchmarkSpherePacking(1, config);
expectedVolume = (4/3) * pi * sum(config.radii.^3);

verifySize(testCase, result.seconds, [1 1]);
verifyTrue(testCase, isfinite(result.seconds));
verifyEqual(testCase, result.acceptedCount, 2);
verifyEqual(testCase, result.sphereVolumeFraction, expectedVolume / 20^3, ...
    'AbsTol', 1e-12);
verifyEqual(testCase, result.medianSeconds, result.seconds, 'AbsTol', 1e-12);
end

function testRequestedRadiiArePreservedAndDoNotOverlap(testCase)
% Catches accepting a sphere with a changed radius or an intersecting pair.
rng(31, 'twister');
mesh = cubeMesh(20);
radii = [0.50; 0.75; 1.00; 1.25];

[assembly, ~, ~, ~, report] = spawnSpheres(mesh, radii, 300, 0.01);

verifyEqual(testCase, report.stopReason, 'completed');
verifyEqual(testCase, assembly(4, :).', radii, 'AbsTol', 1e-12);
for i = 1:size(assembly, 2)
    for j = i + 1:size(assembly, 2)
        separation = norm(assembly(1:3, i) - assembly(1:3, j));
        verifyGreaterThanOrEqual(testCase, separation, ...
            assembly(4, i) + assembly(4, j) - 1e-10);
    end
end
end

function testStopsAndReportsWhenRequestedRadiiDoNotFit(testCase)
% Catches silently returning a partial packing without a capacity warning.
rng(7, 'twister');
mesh = cubeMesh(2.10);
radii = [1.0; 1.0];

[assembly, ~, ~, ~, report] = spawnSpheres(mesh, radii, 1000, 0.01);

verifyLessThan(testCase, size(assembly, 2), numel(radii));
verifyEqual(testCase, report.stopReason, 'capacity_reached');
verifyTrue(testCase, report.capacityWarning);
verifyEqual(testCase, report.requestedCount, 2);
verifyEqual(testCase, report.acceptedCount, size(assembly, 2));
end

function testWritesHeaderedCommaSeparatedCsvOutputs(testCase)
% Catches regressions to headerless or non-CSV result output.
rng(19, 'twister');
outputDirectory = tempname;
cleanup = onCleanup(@() removeOutputDirectory(outputDirectory));
options = struct('outputDirectory', outputDirectory, 'outputPrefix', 'cube');

[assembly, ~, ~, ~, report] = spawnSpheres(cubeMesh(20), [0.5; 0.75], 300, 0.01, options);

verifyEqual(testCase, numel(report.outputFiles), 4);
verifyTrue(testCase, all(endsWith(report.outputFiles, '.csv')));
sphereTable = readtable(report.outputFiles{1});
summaryTable = readtable(report.outputFiles{2});
gridPointTable = readtable(report.outputFiles{3});
gridCellTable = readtable(report.outputFiles{4});
verifyEqual(testCase, sphereTable.Properties.VariableNames, ...
    {'id', 'x', 'y', 'z', 'radius', 'diameter', 'mass'});
verifyEqual(testCase, height(sphereTable), size(assembly, 2));
verifyTrue(testCase, ismember('requested_count', summaryTable.Properties.VariableNames));
verifyEqual(testCase, gridPointTable.Properties.VariableNames, {'point_id', 'x', 'y', 'z'});
verifyEqual(testCase, gridCellTable.Properties.VariableNames, ...
    {'cell_id', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'sphere_count', 'triangle_count'});
verifyGreaterThan(testCase, height(gridCellTable), 0);
verifyEqual(testCase, height(gridPointTable), 8 * height(gridCellTable));
clear cleanup
end

function testReportProvidesCentredAssemblyProperties(testCase)
% Catches returning world-coordinate output without the reported COM data.
rng(23, 'twister');
[assembly, masses, volume, inertia, report] = ...
    spawnSpheres(cubeMesh(20), [0.5; 0.75; 1.0], 300, 0.01);

verifyEqual(testCase, report.boundingBoxDimensions, [20 20 20], 'AbsTol', 1e-12);
verifyEqual(testCase, report.sphereAssemblyVolume, volume, 'AbsTol', 1e-12);
verifyEqual(testCase, report.totalMass, sum(masses), 'AbsTol', 1e-12);
verifyLessThan(testCase, norm(report.centreOfMassAfterShift), 1e-10);
verifyEqual(testCase, assembly(1:3, :) * masses.' / sum(masses), zeros(3,1), 'AbsTol', 1e-10);
verifySize(testCase, inertia, [3 3]);
end

function testPrintsPreprocessingAndQuarterFillProgress(testCase)
% Catches delaying domain/grid status until the final summary or printing duplicate milestones.
outputDirectory = tempname;
cleanup = onCleanup(@() removeOutputDirectory(outputDirectory));
options = struct('outputDirectory', outputDirectory, 'outputPrefix', 'progress', ...
    'randomSeed', 53, 'coordinateFrame', 'world');

text = evalc('spawnSpheres(cubeMesh(20), [0.5; 0.5; 0.5; 0.5], 300, 0.01, options);');

verifyEqual(testCase, count(text, 'Bounding Box Dimensions'), 2);
verifyEqual(testCase, count(text, 'Spatial Grid Discretisation'), 1);
verifyEqual(testCase, count(text, sprintf('\nGrid Cell Edge Length =')), 1);
verifyEqual(testCase, count(text, 'Total Spatial Cells'), 1);
verifyEqual(testCase, count(text, 'Filling Progress: 25%'), 1);
verifyEqual(testCase, count(text, 'Filling Progress: 50%'), 1);
verifyEqual(testCase, count(text, 'Filling Progress: 75%'), 1);
verifyEqual(testCase, count(text, 'Filling Progress: 100%'), 1);

clear cleanup
end

function testCoordinateFrameKeepsSphereAndGridCsvCoordinatesAligned(testCase)
% Catches exporting centred spheres alongside world-coordinate grid points.
worldDirectory = tempname;
centredDirectory = tempname;
cleanup = onCleanup(@() removeOutputDirectory(worldDirectory));
cleanup2 = onCleanup(@() removeOutputDirectory(centredDirectory));
radii = [0.5; 0.75; 1.0];

worldOptions = struct('outputDirectory', worldDirectory, 'outputPrefix', 'cube', ...
    'randomSeed', 47, 'coordinateFrame', 'world');
[worldAssembly, worldMasses, ~, ~, worldReport] = ...
    spawnSpheres(cubeMesh(20), radii, 300, 0.01, worldOptions);

centredOptions = struct('outputDirectory', centredDirectory, 'outputPrefix', 'cube', ...
    'randomSeed', 47, 'coordinateFrame', 'center_of_mass');
[centredAssembly, centredMasses, ~, ~, centredReport] = ...
    spawnSpheres(cubeMesh(20), radii, 300, 0.01, centredOptions);

verifyEqual(testCase, worldAssembly(1:3, :) * worldMasses.' / sum(worldMasses), ...
    worldReport.centreOfMass, 'AbsTol', 1e-10);
verifyEqual(testCase, centredAssembly(1:3, :) * centredMasses.' / sum(centredMasses), ...
    zeros(3, 1), 'AbsTol', 1e-10);

worldSphereTable = readtable(worldReport.outputFiles{1});
centredSphereTable = readtable(centredReport.outputFiles{1});
worldGridTable = readtable(worldReport.outputFiles{3});
centredGridTable = readtable(centredReport.outputFiles{3});
shift = worldReport.centreOfMass.';
verifyEqual(testCase, centredSphereTable{:, {'x', 'y', 'z'}}, ...
    worldSphereTable{:, {'x', 'y', 'z'}} - shift, 'AbsTol', 1e-10);
verifyEqual(testCase, centredGridTable{:, {'x', 'y', 'z'}}, ...
    worldGridTable{:, {'x', 'y', 'z'}} - shift, 'AbsTol', 1e-10);
verifyEqual(testCase, worldReport.coordinateFrame, 'world');
verifyEqual(testCase, centredReport.coordinateFrame, 'center_of_mass');

clear cleanup cleanup2
end

function testSpatialIndexDoesNotAllocateEmptyCells(testCase)
% Catches a dense cell(nx,ny,nz) allocation for large STL-to-radius ratios.
context = spBuildContext(cubeMesh(1e6), 1e-3, 0, 1e-9);

verifyClass(testCase, context.triangleCells, 'containers.Map');
verifyClass(testCase, context.xyCells, 'containers.Map');
verifyLessThan(testCase, context.triangleCells.Count, 1000);
end

function testSpatialGridUsesScalarWidthForNonCubicFaces(testCase)
% Catches concatenating a scalar sphere diameter with a 1-by-3 face extent.
mesh = struct('vertices', [0 0 0; 6 0 0; 0 2 1; 0 0 4], ...
    'faces', [1 2 3; 1 4 2; 1 3 4; 2 4 3]);
context = spBuildContext(mesh, 0.5, 0, 1e-9);

verifySize(testCase, context.cellSize, [1 1]);
verifyEqual(testCase, context.cellSize, 6, 'AbsTol', 1e-12);
end

function testInsidePredicateIsScaleIndependent(testCase)
% Catches treating a length tolerance as an area tolerance for tiny STL meshes.
vertices = 1e-6 * [0 0 0; 1 0 0; 0 1 0; 0 0 1];
faces = [1 3 2; 1 2 4; 1 4 3; 2 3 4];
context = spBuildContext(struct('vertices', vertices, 'faces', faces), 0.05e-6, 0, 1e-9);

verifyTrue(testCase, spPointInside(context, [0.1 0.1 0.1]*1e-6));
verifyFalse(testCase, spPointInside(context, [0.8 0.8 0.8]*1e-6));
end

function testExactPointInsideMatchesScalarReferenceAcrossScales(testCase)
% Catches parity changes while replacing scalar vertical-ray intersections.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
for scale = [1, 1e-6]
    context = spBuildContext(spTestCubeMesh(20*scale), 1.0*scale, 0.01*scale, 1e-9);
    axes = cell(1, 3);
    for dimension = 1:3
        axes{dimension} = context.lower(dimension) + ...
            ((0:context.cellCount(dimension)-1) + 0.5) * context.cellSize;
    end
    [x, y, z] = ndgrid(axes{1}, axes{2}, axes{3});
    cellCentres = [x(:), y(:), z(:)];
    rng(71, 'twister');
    points = scale * [10 10 10; 0 0 0; 20 20 20; 20.1 10 10; 20*rand(100,3)];
    points = [cellCentres; points];
    for id = 1:size(points, 1)
        expected = spReferencePointInside(context, points(id,:));
        verifyEqual(testCase, spExactPointInside(context, points(id,:)), expected);
        verifyEqual(testCase, spPointInside(context, points(id,:)), expected);
    end
end
end

function testTriangleCacheAndBatchedHitsMatchScalarReference(testCase)
% Catches cache or batched closest-feature changes to strict collision logic.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
context = spBuildContext(spTestCubeMesh(20), 1.0, 0.01, 1e-9);
ids = 1:size(context.faces, 1);
for fieldName = {'a', 'b', 'c'}
    value = context.triangles.(fieldName{1});
    verifyClass(testCase, value, 'double');
    verifySize(testCase, value, [size(context.faces,1), 3]);
end
verifyEqual(testCase, context.triangles.a, double(context.vertices(context.faces(:,1),:)));
verifyEqual(testCase, context.triangles.b, double(context.vertices(context.faces(:,2),:)));
verifyEqual(testCase, context.triangles.c, double(context.vertices(context.faces(:,3),:)));
shapeProbe = [10 10 0.1];
expected = spSphereHitsTriangles(context, shapeProbe, 0.2, [1 2]);
verifyEqual(testCase, spSphereHitsTriangles(context, shapeProbe, 0.2, [1; 2]), expected);
verifyEqual(testCase, spSphereHitsTriangles(context, shapeProbe, 0.2, [1 1 2]), expected);
verifyFalse(testCase, spSphereHitsTriangles(context, shapeProbe, 0.2, []));

for centre = [1 10 10; 19 10 10; 10 1 10; 10 19 10; 10 10 1; 10 10 19].'
    for radius = [0.99 1.00 1.01]
        verifyEqual(testCase, spSphereHitsTriangles(context, centre.', radius, ids), ...
            spReferenceSphereHitsTriangles(context, centre.', radius, ids));
    end
end
rng(73, 'twister');
for id = 1:100
    centre = 20*rand(1,3); radius = 0.05 + 2.0*rand;
    verifyEqual(testCase, spSphereHitsTriangles(context, centre, radius, ids), ...
        spReferenceSphereHitsTriangles(context, centre, radius, ids));
end
end

function testBatchedTriangleHitsCoverEveryClosestFeature(testCase)
% Catches changed branch order across vertex, edge, and face regions.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
context = spBuildContext(spTestCubeMesh(20), 1.0, 0.01, 1e-9);
points = [-1 -1 0; 19 21 0; 21 -1 0; 5 5 1; 5 0 1; 20 5 1; 15 5 1];
distances = [sqrt(2), sqrt(2), sqrt(2), 1, 1, 1, 1];
for id = 1:size(points, 1)
    for delta = [-0.1 0.1]
        radius = distances(id) + context.tolerance + delta;
        verifyEqual(testCase, spSphereHitsTriangles(context, points(id,:), radius, 1), ...
            spReferenceSphereHitsTriangles(context, points(id,:), radius, 1));
    end
end
verifyFalse(testCase, spSphereHitsTriangles(context, [10 10 10], context.tolerance, 1));
verifyFalse(testCase, spSphereHitsTriangles(context, [10 10 10], context.tolerance/2, 1));
end

function testOccupancyHandlesConvexConcaveAndCavityMeshes(testCase)
% Catches nonconservative classification of solid, concave, and cavity cells.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
cfg = struct('enabled', true, 'cellSize', 0.5, 'maxCells', 2e6);
convex = spBuildContext(spTestCubeMesh(20), 1.0, 0, 1e-9, cfg);
concave = spBuildContext(spTestConcavePrism(), 0.25, 0, 1e-9, cfg);
hollow = spBuildContext(spTestHollowCubeMesh(), 0.5, 0, 1e-9, cfg);
cases = {convex, [10 10 10], true; convex, [21 10 10], false; ...
         concave, [0.5 2.0 1.5], true; concave, [2.0 2.0 1.5], false; ...
         hollow, [2 10 10], true; hollow, [10 10 10], false; hollow, [21 10 10], false};
for id = 1:size(cases, 1)
    context = cases{id,1}; point = cases{id,2}; expected = cases{id,3};
    verifyEqual(testCase, spExactPointInside(context, point), expected);
    verifyEqual(testCase, spPointInside(context, point), expected);
end
end

function testDisabledOccupancyAndClampedCellIndices(testCase)
% Catches changing the existing four-argument context construction path.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
disabled = spBuildContext(spTestCubeMesh(4), 0.5, 0, 1e-9);
verifyFalse(testCase, disabled.occupancy.enabled);
verifyEmpty(testCase, disabled.occupancy.labels);
point = [2 2 2];
verifyEqual(testCase, spPointInside(disabled, point), ...
    spExactPointInside(disabled, point));

cfg = struct('enabled', true, 'cellSize', 0.5, 'maxCells', 2e6);
enabled = spBuildContext(spTestCubeMesh(4), 0.5, 0, 1e-9, cfg);
occupancy = enabled.occupancy;
verifyEqual(testCase, spOccupancyCellIndex(occupancy, [-10 -10 -10]), [1 1 1]);
verifyEqual(testCase, spOccupancyCellIndex(occupancy, [10 10 10]), occupancy.cellCount);
verifyEqual(testCase, spOccupancyCellIndex(occupancy, occupancy.lower), [1 1 1]);
verifyEqual(testCase, spOccupancyCellIndex(occupancy, enabled.upper), occupancy.cellCount);
end

function testOccupancyGridCapsDenseStorageAndCompletesLabels(testCase)
% Catches allocating the requested fine grid before applying the hard cap.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
cfg = struct('enabled', true, 'cellSize', 1e-4, 'maxCells', 64);
context = spBuildContext(spTestCubeMesh(20), 1.0, 0, 1e-9, cfg);

verifyLessThanOrEqual(testCase, prod(double(context.occupancy.cellCount)), cfg.maxCells);
verifyClass(testCase, context.occupancy.labels, 'uint8');
verifyFalse(testCase, any(context.occupancy.labels(:) == 3));
end

function testEveryNonMarginOccupancyCellMatchesExactParity(testCase)
% Catches a 0/1 shortcut that is not valid throughout its entire cell.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
cfg = struct('enabled', true, 'cellSize', 0.5, 'maxCells', 2e6);
contexts = {spBuildContext(spTestCubeMesh(4), 0.5, 0, 1e-9, cfg), ...
    spBuildContext(spTestConcavePrism(), 0.25, 0, 1e-9, cfg), ...
    spBuildContext(spTestHollowCubeMesh(), 0.5, 0, 1e-9, ...
        struct('enabled', true, 'cellSize', 2.5, 'maxCells', 2e6))};
offsets = [-0.49 0 0.49];
for contextId = 1:numel(contexts)
    context = contexts{contextId};
    occupancy = context.occupancy;
    nonMargin = find(occupancy.labels ~= 2);
    for linearIndex = nonMargin.'
        [ix, iy, iz] = ind2sub(occupancy.cellCount, linearIndex);
        centre = occupancy.lower + ([ix iy iz] - 0.5) * occupancy.cellSize;
        expected = occupancy.labels(linearIndex) == 1;
        verifyEqual(testCase, spExactPointInside(context, centre), expected);
        for dimension = 1:3
            for offset = offsets
                point = centre;
                point(dimension) = point(dimension) + offset * occupancy.cellSize;
                verifyEqual(testCase, spExactPointInside(context, point), expected);
                verifyEqual(testCase, spPointInside(context, point), expected);
            end
        end
    end
end
end

function testSeededPackingMatchesWithOccupancyDisabledAndEnabled(testCase)
% Catches occupancy shortcuts changing the accepted seeded sphere assembly.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'helpers'));
disabledOutput = tempname;
enabledOutput = tempname;
cleanup = onCleanup(@() removeOutputDirectory(disabledOutput));
cleanup2 = onCleanup(@() removeOutputDirectory(enabledOutput));
baseOptions = struct('randomSeed', 53, 'coordinateFrame', 'world', ...
    'maxCompressionSweeps', 8, 'shakeSweeps', 1, 'maxRefillPasses', 1, ...
    'occupancyCellSize', 0.5, 'occupancyMaxCells', 2e6);
disabledOptions = baseOptions;
disabledOptions.occupancyAcceleration = false;
disabledOptions.outputDirectory = disabledOutput;
disabledOptions.outputPrefix = 'disabled';
enabledOptions = baseOptions;
enabledOptions.occupancyAcceleration = true;
enabledOptions.outputDirectory = enabledOutput;
enabledOptions.outputPrefix = 'enabled';

[disabledAssembly, disabledMasses, disabledVolume, disabledInertia, disabledReport] = ...
    spawnSpheres(spTestCubeMesh(20), [0.5; 0.75; 1.0; 1.25], 300, 0.01, disabledOptions);
[enabledAssembly, enabledMasses, enabledVolume, enabledInertia, enabledReport] = ...
    spawnSpheres(spTestCubeMesh(20), [0.5; 0.75; 1.0; 1.25], 300, 0.01, enabledOptions);

verifyEqual(testCase, enabledAssembly, disabledAssembly, 'AbsTol', 1e-12);
verifyEqual(testCase, enabledMasses, disabledMasses, 'AbsTol', 1e-12);
verifyEqual(testCase, enabledVolume, disabledVolume, 'AbsTol', 1e-12);
verifyEqual(testCase, enabledInertia, disabledInertia, 'AbsTol', 1e-12);
verifyEqual(testCase, enabledReport.acceptedCount, disabledReport.acceptedCount);
verifyEqual(testCase, enabledReport.sphereAssemblyVolume / enabledReport.stlVolume, ...
    disabledReport.sphereAssemblyVolume / disabledReport.stlVolume, 'AbsTol', 1e-12);
clear cleanup cleanup2
end

function testCachedInwardNormalsPointIntoTheClosedMesh(testCase)
% Catches cached face normals that no longer use the existing ray-probe rule.
context = spBuildContext(spTestCubeMesh(20), 1.0, 0.01, 1e-9);
probeDistance = max(context.tolerance*100, 1e-8*context.cellSize);
for id = 1:size(context.faces, 1)
    probe = context.faceCentres(id,:) + probeDistance*context.inwardNormals(id,:);
    verifyTrue(testCase, spPointInside(context, probe));
end
end

function testDisabledRefillLeavesStateAndNextRadiusUnchanged(testCase)
% Catches a disabled refill pass that mutates packing state.
context = spBuildContext(spTestCubeMesh(20), 1.0, 0.01, 1e-9);
state = spAddSphere(context, spEmptyState(context), [10 10 10], 1.0);
before = state;
options = struct('maxRefillPasses', 0, 'gravity', [0 0 -1]);

[state, nextRadius] = spRefill(context, state, 1.0, 1, options);

verifyEqual(testCase, state.centres, before.centres);
verifyEqual(testCase, state.radii, before.radii);
verifyEqual(testCase, state.count, before.count);
verifyEqual(testCase, nextRadius, 1);
end

function mesh = cubeMesh(sideLength)
%CUBEMESH Construct a consistently oriented closed cube for deterministic tests.
%Use two triangular faces per side to match the STL mesh representation.
v = [0 0 0; sideLength 0 0; sideLength sideLength 0; 0 sideLength 0; ...
     0 0 sideLength; sideLength 0 sideLength; sideLength sideLength sideLength; 0 sideLength sideLength];
f = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
     2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8];
mesh = struct('vertices', v, 'faces', f);
end

function removeOutputDirectory(pathName)
%REMOVEOUTPUTDIRECTORY Delete one test-only temporary output directory.
if isfolder(pathName)
    rmdir(pathName, 's');
end
end
