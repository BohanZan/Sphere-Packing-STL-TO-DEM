function inside=spExactPointInsideBatch(context,points)
%SPEXACTPOINTINSIDEBATCH Exact parity with immutable data shared per XY bin.
% No random candidates, quantized coordinates, or Q-by-face temporary arrays.
inside=false(size(points,1),1);
if isempty(points), return; end
%Match the checked C++ point-query contract before indexing any batch row.
if size(points,2)~=3 || ~isreal(points) || any(~isfinite(points(:)))
 error('SpherePacking:InvalidPoint','Point coordinates must be finite and real.');
end
points=double(points);
finiteRows=(1:size(points,1)).';
pointsToIndex=points(finiteRows,1:2);
indices=min(max(floor((pointsToIndex-context.lower(1:2))/context.xySize)+1,1),context.xyCount);
if isfield(context,'xyKeySpec')
    encoded=spCellKeys(indices,context.xyKeySpec);
elseif isa(context.xyCells,'containers.Map') && strcmp(context.xyCells.KeyType,'char')
    encoded=cell(size(indices,1),1);
    for k=1:size(indices,1), encoded{k}=sprintf('%d,%d',indices(k,:)); end
else
    encoded=spCellKeys(indices,context.xyCount);
end
[sorted,order]=sort(encoded);
if isempty(sorted), return; end
if iscell(sorted), change=~strcmp(sorted(2:end),sorted(1:end-1)); else, change=sorted(2:end)~=sorted(1:end-1); end
starts=[1;find(change)+1]; ends=[starts(2:end)-1;numel(sorted)];
ray=context.ray; tolerance=context.tolerance;
for bin=1:numel(starts)
    key=sorted(starts(bin)); if iscell(key), key=key{1}; end
    ids=spGridLookup(context.xyCells,key); ids=ids(:);
    if isempty(ids), continue; end
    valid=abs(ray.determinant(ids))>tolerance^2;
    useBounds=isfield(ray,'xyBoundsCells') || isfield(ray,'xyLower');
    if isfield(ray,'xyBoundsCells')
        bounds=ray.xyBoundsCells(key); bounds=bounds(valid,:);
    elseif isfield(ray,'xyLower')
        bounds=[ray.xyLower(ids(valid),:),ray.xyUpper(ids(valid),:)];
    end
    ids=ids(valid);
    if isempty(ids), continue; end
    a=ray.a(ids,:); ab=ray.ab(ids,:); ac=ray.ac(ids,:);
    inverse=ray.inverseDeterminant(ids); zDelta=ray.zDelta(ids,:);
    for position=starts(bin):ends(bin)
        row=finiteRows(order(position)); point=points(row,:);
        if useBounds
            selected=find(point(1)>=bounds(:,1) & point(1)<=bounds(:,3) & ...
                point(2)>=bounds(:,2) & point(2)<=bounds(:,4));
        else
            selected=(1:numel(ids)).';
        end
        if isempty(selected), continue; end
        ap=point(1:2)-a(selected,1:2);
        % Scalar parity stores these expressions in default-double arrays,
        % including when a saved context contains single ray coefficients.
        alpha=double((ap(:,1).*ac(selected,2)-ap(:,2).*ac(selected,1)).*inverse(selected));
        beta=double((ab(selected,1).*ap(:,2)-ab(selected,2).*ap(:,1)).*inverse(selected));
        hit=alpha>=-tolerance & beta>=-tolerance & alpha+beta<=1+tolerance;
        selected=selected(hit);
        % Preserve the scalar function's left-associated height arithmetic.
        z=double(a(selected,3)+alpha(hit).*zDelta(selected,1)+beta(hit).*zDelta(selected,2));
        heights=sort(z(z<point(3)-tolerance));
        inside(row)=~isempty(heights) && mod(1+sum(diff(heights)>tolerance),2)==1;
    end
end
end
