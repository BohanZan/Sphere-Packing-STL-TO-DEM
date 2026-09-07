function indices=spCellIndices(keys,counts)
%SPCELLINDICES Decode occupied keys without converting wide integers to double.
if iscell(keys) && (isempty(keys) || ischar(keys{1}))
    indices=zeros(numel(keys),numel(counts));
    format=[repmat('%d,',1,numel(counts)-1) '%d'];
    for k=1:numel(keys), indices(k,:)=sscanf(keys{k},format).'; end
    return
end
if iscell(keys), keys=[keys{:}]; end
remaining=uint64(keys(:))-uint64(1);
indices=zeros(numel(keys),numel(counts));
for axis=1:numel(counts)
    indices(:,axis)=double(rem(remaining,uint64(counts(axis))))+1;
    remaining=idivide(remaining,uint64(counts(axis)),'floor');
end
end
