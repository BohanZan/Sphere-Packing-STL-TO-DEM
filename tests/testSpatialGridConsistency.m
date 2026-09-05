function tests = testSpatialGridConsistency
%TESTSPATIALGRIDCONSISTENCY Regression tests for shared spatial indexing.
tests = functiontests(localfunctions);
end

function testRayGridUsesMainGridSpacingAndCounts(testCase)
% A tall box makes the old independently fitted XY grid use a different
% spacing from the 3-D triangle/sphere grid.
model = rectangularBoxMesh([0 0 0], [10 7 100]);
context = spBuildContext(model, 0.25, 0, 1e-9);

verifyEqual(testCase, context.xySize, context.cellSize);
verifyEqual(testCase, context.xyCount, context.cellCount(1:2));
end

function testNormalOrientationAutomaticPathMatchesSerialPath(testCase)
% The parallel-capable normal-orientation path must produce exactly the
% same result as the portable serial implementation.
model = rectangularBoxMesh([0 0 0], [4 3 2]);
context = spBuildContext(model, 0.25, 0, 1e-9);
faces = model.faces;
vertices = model.vertices;
rawNormals = cross(vertices(faces(:,2),:) - vertices(faces(:,1),:), ...
    vertices(faces(:,3),:) - vertices(faces(:,1),:), 2);
rawNormals = rawNormals ./ vecnorm(rawNormals, 2, 2);
probeContext = struct('vertices', vertices, 'faces', faces, ...
    'lower', context.lower, 'xySize', context.xySize, ...
    'xyCount', context.xyCount, 'xyCells', {context.xyCells}, ...
    'tolerance', context.tolerance, 'ray', context.ray);
probeDistance = max(context.tolerance * 100, 1e-8 * context.cellSize);

[serialNormals, usedParallel] = spOrientInwardNormals( ...
    context.faceCentres, rawNormals, probeDistance, probeContext, false);
[automaticNormals, ~] = spOrientInwardNormals( ...
    context.faceCentres, rawNormals, probeDistance, probeContext, true);

verifyFalse(testCase, usedParallel);
verifyEqual(testCase, automaticNormals, serialNormals, 'AbsTol', 0);
end

function model = rectangularBoxMesh(lower, upper)
x = [lower(1), upper(1)]; y = [lower(2), upper(2)]; z = [lower(3), upper(3)];
model.vertices = [x(1) y(1) z(1); x(2) y(1) z(1); x(2) y(2) z(1); x(1) y(2) z(1); ...
                  x(1) y(1) z(2); x(2) y(1) z(2); x(2) y(2) z(2); x(1) y(2) z(2)];
model.faces = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
               2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8];
end
