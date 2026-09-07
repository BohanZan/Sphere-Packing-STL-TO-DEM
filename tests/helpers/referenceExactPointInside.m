function inside = referenceExactPointInside(context, point)
%SPEXACTPOINTINSIDE Test a point by the exact downward-ray parity rule.
% Independent all-face oracle: never reuse the production grid or bounds.
ids = (1:size(context.faces,1)).';
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
