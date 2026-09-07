function grid=spBuildStaticGrid(cellKeys,faceIds,keyType)
%SPBUILDSTATICGRID Compact immutable face lists, proportional to occupied cells.
% Sort once by key/face ID instead of repeatedly expanding Map value vectors.
faceIds=faceIds(:);
if strcmp(keyType,'char')
    [sorted,order]=sort(cellKeys(:)); ids=double(faceIds(order));
    if isempty(sorted), starts=[]; else
        starts=[1;find(~strcmp(sorted(2:end),sorted(1:end-1)))+1];
    end
    occupied=sorted(starts); parts=cell(numel(starts),1);
    ends=[starts(2:end)-1;numel(sorted)];
    for k=1:numel(starts), parts{k}=unique(ids(starts(k):ends(k))); end
    ids=vertcat(parts{:}); lengths=cellfun(@numel,parts);
else
    pairs=sortrows([uint64(cellKeys(:)),uint64(faceIds(:))],[1 2]);
    if isempty(pairs), keep=false(0,1); else
        keep=[true;any(pairs(2:end,:)~=pairs(1:end-1,:),2)];
    end
    pairs=pairs(keep,:);
    if isempty(pairs), starts=[]; else
        starts=[1;find(pairs(2:end,1)~=pairs(1:end-1,1))+1];
    end
    occupied=pairs(starts,1); ids=pairs(:,2);
    lengths=diff([starts;size(pairs,1)+1]);
end
directory=containers.Map('KeyType',keyType,'ValueType','double');
if ~isempty(occupied)
    mapKeys=occupied;
    if ~iscell(mapKeys), mapKeys=num2cell(mapKeys); end
    directory=containers.Map(mapKeys,num2cell((1:numel(occupied)).'));
end
if isempty(ids) || max(ids)<=intmax('uint32'), ids=uint32(ids); else, ids=double(ids); end
grid=struct('occupiedKeys',{occupied},'offsets',[1;cumsum(double(lengths(:)))+1], ...
    'faceIds',ids(:),'directory',directory);
end
