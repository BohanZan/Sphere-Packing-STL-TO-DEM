function ids = spHashNeighbours(context, hash, index)
%SPHASHNEIGHBOURS Return unique IDs from INDEX and its 26 neighbours.
%Search the local 3-by-3-by-3 block, clipped at domain boundaries.
%Batch map membership and value retrieval to avoid up to 54 scalar map calls.
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
    %Triangles can cover adjacent cells, so keep the generic query deduplicated.
    ids = sort([parts{:}]);
    ids = ids([true, diff(ids) ~= 0]);
end
end
