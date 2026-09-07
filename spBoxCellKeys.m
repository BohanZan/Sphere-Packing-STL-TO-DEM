function keys=spBoxCellKeys(lo,hi,spec)
%SPBOXCELLKEYS Enumerate a clipped box, retaining legacy x/y/z bucket order.
x=lo(1):hi(1); y=lo(2):hi(2); z=lo(3):hi(3);
indices=[repelem(x,numel(y)*numel(z)).', ...
    repmat(repelem(y,numel(z)),1,numel(x)).', ...
    repmat(z,1,numel(x)*numel(y)).'];
keys=spCellKeys(indices,spec);
if ~iscell(keys), keys=num2cell(keys); end
keys=keys(:).';
end
