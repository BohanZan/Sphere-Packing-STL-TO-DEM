function spHashInsert(hash, index, id)
%SPHASHINSERT Append ID to the sparse 3-D cell INDEX.
key = sprintf('%d,%d,%d', index(1), index(2), index(3));
if isKey(hash, key)
    ids = hash(key);
else
    ids = [];
end
hash(key) = [ids id];
end
