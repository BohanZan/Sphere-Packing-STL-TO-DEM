function ids = spSphereNeighbours(context, hash, index)
%SPSPHERENEIGHBOURS Return sphere IDs from INDEX and its 26 neighbours.
%Each sphere is indexed in one spatial cell, so this query needs no unique.
%Batch the map membership and value lookups; calling containers.Map once per
%neighbour cell otherwise dominates this very small (at most 27-cell) query.
cellCount = context.cellCount;
lo = max(index - 1, 1);
hi = min(index + 1, cellCount);
keys = spBoxCellKeys(lo,hi,context.cellKeySpec);
present = isKey(hash, keys);
if ~any(present)
    ids = [];
else
    parts = values(hash, keys(present));
    ids = [parts{:}];
end
end
