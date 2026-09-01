function tests = testSpawnSpheres
% Regression tests for the random non-overlapping sphere insertion core.
tests = functiontests(localfunctions);
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

verifyEqual(testCase, numel(report.outputFiles), 2);
verifyTrue(testCase, all(endsWith(report.outputFiles, '.csv')));
sphereTable = readtable(report.outputFiles{1});
summaryTable = readtable(report.outputFiles{2});
verifyEqual(testCase, sphereTable.Properties.VariableNames, ...
    {'id', 'x', 'y', 'z', 'radius', 'diameter', 'mass'});
verifyEqual(testCase, height(sphereTable), size(assembly, 2));
verifyTrue(testCase, ismember('requested_count', summaryTable.Properties.VariableNames));
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

function mesh = cubeMesh(sideLength)
v = [0 0 0; sideLength 0 0; sideLength sideLength 0; 0 sideLength 0; ...
     0 0 sideLength; sideLength 0 sideLength; sideLength sideLength sideLength; 0 sideLength sideLength];
f = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
     2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8];
mesh = struct('vertices', v, 'faces', f);
end

function removeOutputDirectory(pathName)
if isfolder(pathName)
    rmdir(pathName, 's');
end
end
