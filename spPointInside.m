function inside = spPointInside(context, point)
%SPPOINTINSIDE Test a point against the closed STL by downward ray parity.
index = floor((point(1:2)-context.lower(1:2))/context.xySize) + 1;
index = min(max(index, 1), context.xyCount);
key = sprintf('%d,%d', index(1), index(2));
if ~isKey(context.xyCells, key)
    inside = false;
    return;
end
ids = context.xyCells(key);
heights = zeros(1, numel(ids));
count = 0;
for id = ids
    [hit, zHit] = verticalTriangleHit(point, context.vertices(context.faces(id,:),:), context.tolerance);
    if hit && zHit < point(3) - context.tolerance
        count = count + 1;
        heights(count) = zHit;
    end
end
if count == 0
    inside = false;
    return;
end
heights = sort(heights(1:count));
inside = mod(1 + sum(diff(heights) > context.tolerance), 2) == 1;
end

function [hit, zHit] = verticalTriangleHit(point, triangle, tolerance)
a = triangle(1,1:2);
ab = triangle(2,1:2) - a;
ac = triangle(3,1:2) - a;
ap = point(1:2) - a;
determinant = ab(1)*ac(2) - ab(2)*ac(1);
if abs(determinant) <= tolerance
    hit = false;
    zHit = NaN;
    return;
end
alpha = (ap(1)*ac(2) - ap(2)*ac(1)) / determinant;
beta = (ab(1)*ap(2) - ab(2)*ap(1)) / determinant;
hit = alpha >= -tolerance && beta >= -tolerance && alpha+beta <= 1+tolerance;
zHit = triangle(1,3) + alpha*(triangle(2,3)-triangle(1,3)) + beta*(triangle(3,3)-triangle(1,3));
end
