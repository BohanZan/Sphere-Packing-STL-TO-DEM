function ids=spGridLookup(grid,cellKeys)
%SPGRIDLOOKUP Concatenate occupied static spans in requested cell order.
% Legacy Maps remain readable for frozen baselines and saved contexts.
ids=[];
if isempty(cellKeys), return; end
if ischar(cellKeys), cellKeys={cellKeys}; elseif ~iscell(cellKeys), cellKeys=num2cell(cellKeys); end
if isa(grid,'containers.Map')
    present=isKey(grid,cellKeys);
    if any(present), parts=values(grid,cellKeys(present)); ids=[parts{:}]; end
    return
end
if numel(cellKeys)==1
    key=cellKeys{1};
    if ~isKey(grid.directory,key), return; end
    row=grid.directory(key);
    ids=double(grid.faceIds(grid.offsets(row):grid.offsets(row+1)-1)).';
    return
end
present=isKey(grid.directory,cellKeys);
if ~any(present), return; end
parts=values(grid.directory,cellKeys(present)); rows=[parts{:}];
starts=grid.offsets(rows); lengths=grid.offsets(rows+1)-starts;
starts=starts(:).'; lengths=lengths(:).';
positions=repelem(starts-cumsum(lengths)+lengths,lengths)+(0:sum(lengths)-1);
ids=double(grid.faceIds(positions)).';
end
