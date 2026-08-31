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

function mesh = cubeMesh(sideLength)
v = [0 0 0; sideLength 0 0; sideLength sideLength 0; 0 sideLength 0; ...
     0 0 sideLength; sideLength 0 sideLength; sideLength sideLength sideLength; 0 sideLength sideLength];
f = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
     2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8];
mesh = struct('vertices', v, 'faces', f);
end
