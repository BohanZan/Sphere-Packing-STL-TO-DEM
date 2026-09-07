function context = spBuildContext(model, maxRadius, buffer, tolerance, occupancyOptions)
%SPBUILDCONTEXT Preprocess STL triangles into 2-D and 3-D uniform grids.
%The resulting sparse hashes accelerate sphere, triangle and ray queries.

% Keep direct callers on the existing sparse-only construction path.
if nargin < 5 || isempty(occupancyOptions)
    occupancyOptions = struct('enabled', false, 'cellSize', [], 'maxCells', 2e6);
end
occupancyOptions = spNormaliseOccupancyOptions(occupancyOptions);
validateattributes(maxRadius,{'numeric'},{'real','finite','scalar','positive'});
if ~(isnumeric(buffer) && isreal(buffer) && isscalar(buffer) && ...
        isfinite(buffer) && buffer>=0 && isfinite(2*maxRadius+buffer))
    error('SpherePacking:InvalidGridBuffer', ...
        'buffer must be finite and nonnegative, with finite grid spacing.');
end

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

% Cell spacing >= the largest sphere diameter guarantees complete local
% sphere queries. Triangles have no size restriction; a conservative finite
% surface distance band covers every possible sphere-centre contact cell.
cellSize = 2 * maxRadius + buffer;
count = max(1, ceil((upper - lower) / cellSize));
[~,cellKeySpec]=spCellKeys(ones(1,3),count);
[~,xyKeySpec]=spCellKeys(ones(1,2),count(1:2));

% A sparse hash is essential: the geometric grid can have billions of
% possible cells while only cells touched by STL faces need storage.
triangleKeyParts=cell(size(faces,1),1); triangleIdParts=triangleKeyParts;
grid=struct('lower',lower,'cellSize',cellSize,'cellCount',count);
rasterStats=struct('enumerated',0,'references',0);
for id = 1:size(faces,1)
    tri = vertices(faces(id,:),:);
    [indices,~,stats]=spTriangleCellCandidates(tri,grid,(1+sqrt(3))*cellSize/2,tolerance);
    rasterStats.enumerated=rasterStats.enumerated+stats.enumerated;
    rasterStats.references=rasterStats.references+size(indices,1);
    triangleKeyParts{id}=spCellKeys(indices,cellKeySpec);
    triangleIdParts{id}=repmat(id,size(indices,1),1);
end
triCells=spBuildStaticGrid(vertcat(triangleKeyParts{:}),vertcat(triangleIdParts{:}),cellKeySpec.keyType);
% Eq.10 covers every sphere-centre cell (r <= h/2), not just surface cells.
% This is an index invariant, not a user-selectable approximation.
triCells.centreCoverage=true;
clear triangleKeyParts triangleIdParts

%Store the downward-ray projection in the XY layers of the same main grid.
%It remains a 2-D hash to avoid duplicating faces along Z, but shares the
%main grid origin, edge length and XY cell partition exactly.
xySize = cellSize;
xyCount = count(1:2);
xyKeyParts=cell(size(faces,1),1); xyIdParts=xyKeyParts;

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
% The parity test accepts barycentric coordinates slightly outside a face.
% Expand projected bounds for that same tolerance and cancellation in thin
% projections. These bounds only reject impossible hits; parity stays exact.
projectionCondition = (sum(abs(ab),2)+sum(abs(ac),2)).^2 ./ abs(determinant);
roundoff = 64 * eps(max(abs(vertices(:)))) * max(1, projectionCondition);
padding = 2*tolerance*(abs(ab)+abs(ac)) + roundoff;
ray.xyLower = min(min(a(:,1:2),b(:,1:2)),c(:,1:2)) - padding;
ray.xyUpper = max(max(a(:,1:2),b(:,1:2)),c(:,1:2)) + padding;
% The grid must cover exactly the same expanded projection as the ray test.
% Invalid projections cannot contribute a hit and need no XY references.
for id=find(abs(determinant)>tolerance^2).'
    lo=toXY(ray.xyLower(id,:)); hi=toXY(ray.xyUpper(id,:));
    [ix,iy]=ndgrid(lo(1):hi(1),lo(2):hi(2));
    xyKeyParts{id}=spCellKeys([ix(:),iy(:)],xyKeySpec);
    xyIdParts{id}=repmat(id,numel(ix),1);
end
xyCells=spBuildStaticGrid(vertcat(xyKeyParts{:}),vertcat(xyIdParts{:}),xyKeySpec.keyType);
clear xyKeyParts xyIdParts

%Cache the face geometry and orient each normal by the existing ray probe.
faceVertices = reshape(vertices(faces.', :), 3, size(faces,1), 3);
faceCentres = reshape(mean(faceVertices, 1), size(faces,1), 3);
rawNormals = cross(vertices(faces(:,2),:) - vertices(faces(:,1),:), ...
    vertices(faces(:,3),:) - vertices(faces(:,1),:), 2);
rawNormals = rawNormals ./ vecnorm(rawNormals, 2, 2);
probeContext = struct('vertices',vertices,'faces',faces,'lower',lower, ...
    'xySize',xySize,'xyCount',xyCount,'xyCells',{xyCells}, ...
    'tolerance',tolerance,'ray',ray,'xyKeySpec',xyKeySpec);
probeDistance = max(tolerance*100, 1e-8*cellSize);
inwardNormals = spOrientInwardNormals( ...
    faceCentres, rawNormals, probeDistance, probeContext);

%Collect all preprocessed geometry and spatial indexing information.
context=struct('vertices',vertices,'faces',faces,'lower',lower,'upper',upper,...
    'cellSize',cellSize,'cellCount',count,'triangleCells',{triCells},...
    'xySize',xySize,'xyCount',xyCount,'xyCells',{xyCells},'tolerance',tolerance, ...
    'faceCentres',faceCentres,'inwardNormals',inwardNormals,'ray',ray, ...
    'triangles',triangles,'cellKeySpec',cellKeySpec,'xyKeySpec',xyKeySpec, ...
    'rasterStats',rasterStats);
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
        key = spCellKeys(index,cellKeySpec);
        if iscell(key), key=key{1}; end
    end
    function key = xyKey(index)
        key = spCellKeys(index,xyKeySpec);
        if iscell(key), key=key{1}; end
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
