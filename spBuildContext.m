function context = spBuildContext(model, maxRadius, buffer, tolerance)
%SPBUILDCONTEXT Preprocess STL triangles into 2-D and 3-D uniform grids.
[vertices, faces] = readMesh(model);
lower = min(vertices, [], 1); upper = max(vertices, [], 1);
cellSize = 2 * maxRadius + buffer;
count = max(1, ceil((upper - lower) / cellSize));
triCells = cell(count(1), count(2), count(3));
for id = 1:size(faces,1)
    tri = vertices(faces(id,:),:);
    lo = toCell(min(tri,[],1)); hi = toCell(max(tri,[],1));
    for ix = lo(1):hi(1)
        for iy = lo(2):hi(2)
            for iz = lo(3):hi(3)
                triCells{ix,iy,iz}(end+1) = id;
            end
        end
    end
end
span = max(upper(1:2)-lower(1:2));
xySize = max(maxRadius, span / max(1, ceil(span / cellSize)));
xyCount = max(1, ceil((upper(1:2)-lower(1:2)) / xySize));
xyCells = cell(xyCount(1),xyCount(2));
for id=1:size(faces,1)
    tri=vertices(faces(id,:),1:2); lo=toXY(min(tri,[],1)); hi=toXY(max(tri,[],1));
    for ix = lo(1):hi(1)
        for iy = lo(2):hi(2)
            xyCells{ix,iy}(end+1)=id;
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
