function ids = spHashNeighbours(context, hash, index)
%SPHASHNEIGHBOURS Return unique IDs from INDEX and its 26 neighbours.
ids = [];
for ix = max(1, index(1)-1):min(context.cellCount(1), index(1)+1)
    for iy = max(1, index(2)-1):min(context.cellCount(2), index(2)+1)
        for iz = max(1, index(3)-1):min(context.cellCount(3), index(3)+1)
            key = sprintf('%d,%d,%d', ix, iy, iz);
            if isKey(hash, key)
                ids = [ids hash(key)]; %#ok<AGROW>
            end
        end
    end
end
ids = unique(ids);
end
