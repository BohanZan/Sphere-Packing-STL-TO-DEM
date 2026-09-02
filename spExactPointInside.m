function inside = spExactPointInside(context, point)
%SPEXACTPOINTINSIDE Test a point by the exact downward-ray parity rule.
%Locate the point's XY ray cell and reject empty projected regions quickly.
index = floor((point(1:2)-context.lower(1:2))/context.xySize) + 1;
index = min(max(index, 1), context.xyCount);
key = sprintf('%d,%d', index(1), index(2));
if ~isKey(context.xyCells, key)
    inside = false;
    return;
end

%Evaluate the selected finite projected triangles in one face-indexed batch.
ids = context.xyCells(key);
ids = ids(:);
if isempty(ids)
    inside = false;
    return;
end
ray = context.ray;
tolerance = context.tolerance;
determinant = ray.determinant(ids);
valid = abs(determinant) > tolerance^2;
alpha = zeros(size(ids));
beta = zeros(size(ids));
zHit = zeros(size(ids));
if any(valid)
    validIds = ids(valid);
    ap = point(1:2) - ray.a(validIds,1:2);
    alpha(valid) = (ap(:,1).*ray.ac(validIds,2) - ap(:,2).*ray.ac(validIds,1)) ...
        .* ray.inverseDeterminant(validIds);
    beta(valid) = (ray.ab(validIds,1).*ap(:,2) - ray.ab(validIds,2).*ap(:,1)) ...
        .* ray.inverseDeterminant(validIds);
    valid = valid & alpha >= -tolerance & beta >= -tolerance & ...
        alpha + beta <= 1 + tolerance;
    zHit(valid) = ray.a(ids(valid),3) + alpha(valid).*ray.zDelta(ids(valid),1) + ...
        beta(valid).*ray.zDelta(ids(valid),2);
end

%Retain the scalar implementation's strict z cutoff and duplicate-height parity.
heights = sort(zHit(valid & zHit < point(3) - tolerance));
inside = ~isempty(heights) && mod(1 + sum(diff(heights) > tolerance), 2) == 1;
end
