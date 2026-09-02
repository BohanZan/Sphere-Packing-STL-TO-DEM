function inside = spReferencePointInside(context, point)
%SPREFERENCEPOINTINSIDE Scalar pre-vectorization point-in-solid oracle.
%Locate the point's XY ray cell and reject empty projected regions quickly.
index = floor((point(1:2)-context.lower(1:2))/context.xySize) + 1;
index = min(max(index, 1), context.xyCount);
key = sprintf('%d,%d', index(1), index(2));
if ~isKey(context.xyCells, key)
    inside = false;
    return;
end

%Collect intersections strictly below the point along the vertical ray.
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

%Remove duplicate face-edge intersections before applying odd-even parity.
if count == 0
    inside = false;
    return;
end
heights = sort(heights(1:count));
inside = mod(1 + sum(diff(heights) > context.tolerance), 2) == 1;
end

function [hit, zHit] = verticalTriangleHit(point, triangle, tolerance)
%VERTICALTRIANGLEHIT Intersect a vertical ray with a finite projected triangle.
%Form barycentric coordinates in XY; the determinant has area units.
a = triangle(1,1:2);
ab = triangle(2,1:2) - a;
ac = triangle(3,1:2) - a;
ap = point(1:2) - a;
determinant = ab(1)*ac(2) - ab(2)*ac(1);

%Ignore numerically degenerate projections before dividing by their area.
if abs(determinant) <= tolerance^2
    hit = false;
    zHit = NaN;
    return;
end

%Accept projected barycentric coordinates inside the finite triangle.
alpha = (ap(1)*ac(2) - ap(2)*ac(1)) / determinant;
beta = (ab(1)*ap(2) - ab(2)*ap(1)) / determinant;
hit = alpha >= -tolerance && beta >= -tolerance && alpha+beta <= 1+tolerance;

%Interpolate the vertical intersection height from the same barycentric point.
zHit = triangle(1,3) + alpha*(triangle(2,3)-triangle(1,3)) + beta*(triangle(3,3)-triangle(1,3));
end
