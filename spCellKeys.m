function [keys, spec] = spCellKeys(indices, spec)
%SPCELLKEYS Exact sparse-grid keys; never allocate the geometric grid.
% A cached SPEC avoids repeated range/overflow checks on clipped hot queries.
if ~isstruct(spec)
    counts=double(spec(:).');
    if any(~isfinite(counts) | counts<1 | counts~=fix(counts) | counts>flintmax) || ...
            size(indices,2)~=numel(counts) || ...
            any(any(~isfinite(indices) | indices<1 | indices>counts | indices~=fix(indices)))
        error('SpherePacking:InvalidCellIndex','Cell coordinates must be in the integer grid.');
    end
    strides=ones(1,numel(counts),'uint64'); total=uint64(1); numeric=true;
    for axis=1:numel(counts)
        strides(axis)=total;
        if total>idivide(intmax('uint64'),uint64(counts(axis)),'floor')
            numeric=false; break;
        end
        total=total*uint64(counts(axis));
    end
    spec=struct('counts',counts,'strides',strides,'numeric',numeric, ...
        'small',numeric && total<=uint64(flintmax),'keyType','uint64');
    if ~numeric, spec.keyType='char'; end
end
if spec.small
    % A single-precision mesh can produce single indices; doing this
    % arithmetic before promotion would merge adjacent keys above 2^24.
    keys=uint64(1+(double(indices)-1)*double(spec.strides(:)));
elseif spec.numeric
    keys=ones(size(indices,1),1,'uint64');
    for axis=1:size(indices,2)
        keys=keys+(uint64(indices(:,axis))-uint64(1))*spec.strides(axis);
    end
else
    keys=cell(size(indices,1),1);
    format=[repmat('%d,',1,size(indices,2)-1) '%d'];
    for row=1:size(indices,1), keys{row}=sprintf(format,indices(row,:)); end
end
end
