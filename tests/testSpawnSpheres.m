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
