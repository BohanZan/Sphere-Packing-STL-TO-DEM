function context = spBuildContext(model, maxRadius, buffer, tolerance)
%SPBUILDCONTEXT Preprocess STL triangles into 2-D and 3-D uniform grids.
[vertices, faces] = readMesh(model);
lower = min(vertices, [], 1); upper = max(vertices, [], 1);
triangleExtent = zeros(size(faces,1), 3);
for id = 1:size(faces,1)
    triangle = vertices(faces(id,:), :);
    triangleExtent(id,:) = max(triangle, [], 1) - min(triangle, [], 1);
end
% Keep the grid sparse in both storage and construction. A cell may hold
% many spheres; making it at least as wide as the largest face prevents a
% large STL facet from being expanded into billions of indexed cells.
cellSize = max(2 * maxRadius + buffer, max(triangleExtent(:)));
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
span = max(upper(1:2)-lower(1:2));
xySize = max(maxRadius, span / max(1, ceil(span / cellSize)));
xyCount = max(1, ceil((upper(1:2)-lower(1:2)) / xySize));
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
context=struct('vertices',vertices,'faces',faces,'lower',lower,'upper',upper,...
    'cellSize',cellSize,'cellCount',count,'triangleCells',{triCells},...
    'xySize',xySize,'xyCount',xyCount,'xyCells',{xyCells},'tolerance',tolerance);
    function index=toCell(p)
        index=min(max(floor((p-lower)/cellSize)+1,1),count);
    end
    function index=toXY(p)
        index=min(max(floor((p-lower(1:2))/xySize)+1,1),xyCount);
    end
    function key = cellKey(index)
        key = sprintf('%d,%d,%d', index(1), index(2), index(3));
    end
    function key = xyKey(index)
        key = sprintf('%d,%d', index(1), index(2));
    end
end

function [vertices, faces] = readMesh(model)
if isstruct(model)
    vertices=model.vertices; faces=model.faces;
else
    [vertices,faces,~,~]=stlRead(char(model));
end
validateattributes(vertices,{'numeric'},{'2d','ncols',3,'finite','real'});
validateattributes(faces,{'numeric'},{'2d','ncols',3,'positive','integer'});
end
