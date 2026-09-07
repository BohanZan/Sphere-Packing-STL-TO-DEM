function spHashInsert(hash, index, id, spec)
%SPHASHINSERT Append ID to the sparse 3-D cell INDEX.
% The grid encoding is immutable; sphere buckets remain independent ID vectors.
key = spCellKeys(index,spec);
if iscell(key), key=key{1}; end

%Preserve prior members when the cell is already occupied.
if isKey(hash, key)
    ids = hash(key);
else
    ids = [];
end
hash(key) = [ids id];
end
