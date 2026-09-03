function ids = spSphereNeighbours(context, hash, index)
%SPSPHERENEIGHBOURS Return sphere IDs from INDEX and its 26 neighbours.
%Each sphere is indexed in one spatial cell, so this query needs no unique.
%Batch the map membership and value lookups; calling containers.Map once per
%neighbour cell otherwise dominates this very small (at most 27-cell) query.
cellCount = context.cellCount;
lo = max(index - 1, 1);
hi = min(index + 1, cellCount);
x = lo(1):hi(1);
y = lo(2):hi(2);
z = lo(3):hi(3);
xKeys = cellstr(string(x));
yKeys = cellstr(string(y));
zKeys = cellstr(string(z));
keys = strcat(repelem(xKeys, numel(y) * numel(z)), ',', ...
    repmat(repelem(yKeys, numel(z)), 1, numel(x)), ',', ...
    repmat(zKeys, 1, numel(x) * numel(y)));
present = isKey(hash, keys);
if ~any(present)
    ids = [];
else
    parts = values(hash, keys(present));
    ids = [parts{:}];
end
end
