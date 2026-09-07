function ids = spHashNeighbours(context, hash, index)
%SPHASHNEIGHBOURS Return unique IDs from INDEX and its 26 neighbours.
%Search the local 3-by-3-by-3 block, clipped at domain boundaries.
%Batch map membership and value retrieval to avoid up to 54 scalar map calls.
cellCount = context.cellCount;
lo = max(index - 1, 1);
hi = min(index + 1, cellCount);
keys = spBoxCellKeys(lo,hi,context.cellKeySpec);
ids=spGridLookup(hash,keys);
if ~isempty(ids)
    %Triangles can cover adjacent cells, so keep the generic query deduplicated.
    ids = sort(ids);
    ids = ids([true, diff(ids) ~= 0]);
end
end
