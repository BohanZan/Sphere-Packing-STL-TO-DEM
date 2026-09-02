function ids = spSphereNeighbours(context, hash, index)
%SPSPHERENEIGHBOURS Return sphere IDs from INDEX and its 26 neighbours.
%Each sphere is indexed in one spatial cell, so this query needs no unique.
parts = cell(27, 1);
partCount = 0;
for ix = max(1, index(1)-1):min(context.cellCount(1), index(1)+1)
    for iy = max(1, index(2)-1):min(context.cellCount(2), index(2)+1)
        for iz = max(1, index(3)-1):min(context.cellCount(3), index(3)+1)
            key = sprintf('%d,%d,%d', ix, iy, iz);
            if isKey(hash, key)
                partCount = partCount + 1;
                parts{partCount} = hash(key);
            end
        end
    end
end
if partCount == 0
    ids = [];
else
    ids = [parts{1:partCount}];
end
end
