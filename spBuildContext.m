function context = spBuildContext(model, maxRadius, buffer, tolerance, occupancyOptions)
%SPBUILDCONTEXT Preprocess STL triangles into 2-D and 3-D uniform grids.
%The resulting sparse hashes accelerate sphere, triangle and ray queries.

% Keep direct callers on the existing sparse-only construction path.
if nargin < 5 || isempty(occupancyOptions)
    occupancyOptions = struct('enabled', false, 'cellSize', [], 'maxCells', 2e6);
end
occupancyOptions = spNormaliseOccupancyOptions(occupancyOptions);

%Read and bound the closed surface before selecting numerical scales.
[vertices, faces] = readMesh(model);
lower = min(vertices, [], 1); upper = max(vertices, [], 1);
modelScale = max(upper - lower);
dimensions = upper - lower;
fprintf('\n========================================\n');
fprintf('Geometry Preprocessing\n');
fprintf('Bounding Box Dimensions Lx=%.8g; Ly=%.8g; Lz=%.8g\n', ...
    dimensions(1), dimensions(2), dimensions(3));

% The public tolerance is relative to model size. Store a length tolerance
% internally so the same options work for micrometre and metre STL files.
tolerance = max(tolerance * modelScale, 64 * eps(max(abs(vertices(:)))));

%Measure face extents to keep each 3-D hash cell physically meaningful.
triangleExtent = zeros(size(faces,1), 3);
for id = 1:size(faces,1)
    triangle = vertices(faces(id,:), :);
    triangleExtent(id,:) = max(triangle, [], 1) - min(triangle, [], 1);
end
% Prefer a moderately finer grid when a coarse STL contains isolated large
% faces. Each face is still indexed in every cell it spans, so this changes
% candidate-list size only, not geometric coverage.
largeFaceCellFloor = 0.5 * max(triangleExtent(:));
cellSize = max(2 * maxRadius + buffer, largeFaceCellFloor);
count = max(1, ceil((upper - lower) / cellSize));

% A sparse hash is essential: the geometric grid can have billions of
% possible cells while only cells touched by STL faces need storage.
triCells = containers.Map('KeyType', 'char', 'ValueType', 'any');
for id = 1:size(faces,1)
    tri = vertices(faces(id,:),:);
    lo = toCell(min(tri,[],1)); hi = toCell(max(tri,[],1));
    for ix = lo(1):hi(1)
        for iy = lo(2):hi(2)
            for iz = lo(3):hi(3)
                key = cellKey([ix iy iz]);
                if isKey(triCells, key), ids = triCells(key); else, ids = []; end
                triCells(key) = [ids id];
            end
        end
    end
end

%Store the downward-ray projection in the XY layers of the same main grid.
%It remains a 2-D hash to avoid duplicating faces along Z, but shares the
%main grid origin, edge length and XY cell partition exactly.
xySize = cellSize;
xyCount = count(1:2);
xyCells = containers.Map('KeyType', 'char', 'ValueType', 'any');
for id=1:size(faces,1)
    tri=vertices(faces(id,:),1:2); lo=toXY(min(tri,[],1)); hi=toXY(max(tri,[],1));
    for ix = lo(1):hi(1)
        for iy = lo(2):hi(2)
            key = xyKey([ix iy]);
            if isKey(xyCells, key), ids = xyCells(key); else, ids = []; end
            xyCells(key) = [ids id];
        end
    end
end

%Precompute the face-indexed coefficients used by exact vertical ray tests.
a = vertices(faces(:,1),:);
b = vertices(faces(:,2),:);
c = vertices(faces(:,3),:);
triangles = struct('a', double(a), 'b', double(b), 'c', double(c));
ab = b(:,1:2) - a(:,1:2);
ac = c(:,1:2) - a(:,1:2);
determinant = ab(:,1).*ac(:,2) - ab(:,2).*ac(:,1);
ray = struct('a', a, 'ab', ab, 'ac', ac, 'determinant', determinant, ...
    'inverseDeterminant', 1./determinant, ...
    'zDelta', [b(:,3)-a(:,3), c(:,3)-a(:,3)]);

%Cache the face geometry and orient each normal by the existing ray probe.
faceVertices = reshape(vertices(faces.', :), 3, size(faces,1), 3);
faceCentres = reshape(mean(faceVertices, 1), size(faces,1), 3);
rawNormals = cross(vertices(faces(:,2),:) - vertices(faces(:,1),:), ...
    vertices(faces(:,3),:) - vertices(faces(:,1),:), 2);
rawNormals = rawNormals ./ vecnorm(rawNormals, 2, 2);
probeContext = struct('vertices',vertices,'faces',faces,'lower',lower, ...
    'xySize',xySize,'xyCount',xyCount,'xyCells',{xyCells}, ...
    'tolerance',tolerance,'ray',ray);
probeDistance = max(tolerance*100, 1e-8*cellSize);
inwardNormals = spOrientInwardNormals( ...
    faceCentres, rawNormals, probeDistance, probeContext);

%Collect all preprocessed geometry and spatial indexing information.
context=struct('vertices',vertices,'faces',faces,'lower',lower,'upper',upper,...
    'cellSize',cellSize,'cellCount',count,'triangleCells',{triCells},...
    'xySize',xySize,'xyCount',xyCount,'xyCells',{xyCells},'tolerance',tolerance, ...
    'faceCentres',faceCentres,'inwardNormals',inwardNormals,'ray',ray, ...
    'triangles',triangles);
if occupancyOptions.enabled
    context.occupancy = spBuildOccupancyGrid(context, maxRadius, occupancyOptions);
else
    context.occupancy = struct('enabled', false, 'lower', lower, 'cellSize', NaN, ...
        'cellCount', zeros(1,3), 'labels', zeros(0,0,0, 'uint8'));
end
fprintf('Spatial Grid Discretisation\n');
fprintf('Grid Cell Edge Length = %.8g\n', cellSize);
fprintf('Grid Counts Nx=%d; Ny=%d; Nz=%d\n', count(1), count(2), count(3));
fprintf('Total Spatial Cells = %.0f\n', prod(double(count)));
fprintf('Ray Grid Cell Edge Length = %.8g\n', xySize);
fprintf('Ray Grid Counts Nx=%d; Ny=%d; Total=%d\n', ...
    xyCount(1), xyCount(2), prod(double(xyCount)));
fprintf('========================================\n');

%Map a 3-D coordinate to its clamped sparse-hash cell index.
    function index=toCell(p)
        index=min(max(floor((p-lower)/cellSize)+1,1),count);
    end

%Map an XY coordinate to its clamped downward-ray cell index.
    function index=toXY(p)
        index=min(max(floor((p-lower(1:2))/xySize)+1,1),xyCount);
    end

%Encode 3-D and XY integer indices as containers.Map keys.
    function key = cellKey(index)
        key = sprintf('%d,%d,%d', index(1), index(2), index(3));
    end
    function key = xyKey(index)
        key = sprintf('%d,%d', index(1), index(2));
    end
end

function options = spNormaliseOccupancyOptions(options)
%SPNORMALISEOCCUPANCYOPTIONS Validate direct occupancy-grid construction options.
required = {'enabled', 'cellSize', 'maxCells'};
for id = 1:numel(required)
    if ~isfield(options, required{id})
        error('SpherePacking:InvalidOccupancyOptions', ...
            'occupancyOptions.%s is required.', required{id});
    end
end
if ~(isscalar(options.enabled) && (islogical(options.enabled) || ...
        (isnumeric(options.enabled) && isfinite(options.enabled) && ...
        any(options.enabled == [0 1]))))
    error('SpherePacking:InvalidOccupancyAcceleration', ...
        'occupancy acceleration must be a scalar logical value.');
end
options.enabled = logical(options.enabled);
if ~isempty(options.cellSize)
    validateattributes(options.cellSize, {'numeric'}, ...
        {'real','finite','scalar','positive'});
end
validateattributes(options.maxCells, {'numeric'}, ...
    {'real','finite','scalar','integer','positive','<=',double(intmax('uint32'))});
end

function [vertices, faces] = readMesh(model)
%READMESH Accept either an in-memory mesh structure or an STL filename.
%Validate the common vertex and triangular-face representation.
if isstruct(model)
    vertices=model.vertices; faces=model.faces;
else
    [vertices,faces,~,~]=stlRead(char(model));
end
validateattributes(vertices,{'numeric'},{'2d','ncols',3,'finite','real'});
validateattributes(faces,{'numeric'},{'2d','ncols',3,'positive','integer'});
end
