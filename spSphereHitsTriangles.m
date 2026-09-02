function intersects = spSphereHitsTriangles(context, centre, radius, triangleIds)
%SPSPHEREHITSTRIANGLES Check finite-triangle penetration, including edges/vertices.
%Evaluate exact finite-triangle distances in face-ID-indexed batches.
intersects = false;
limit = radius - context.tolerance;
if isempty(triangleIds) || limit <= 0
    return;
end

ids = triangleIds(:);
a = context.triangles.a(ids,:);
b = context.triangles.b(ids,:);
c = context.triangles.c(ids,:);
ab = b - a;
ac = c - a;
ap = bsxfun(@minus, centre, a);
bp = bsxfun(@minus, centre, b);
cp = bsxfun(@minus, centre, c);
d1 = sum(ab .* ap, 2); d2 = sum(ac .* ap, 2);
d3 = sum(ab .* bp, 2); d4 = sum(ac .* bp, 2);
d5 = sum(ab .* cp, 2); d6 = sum(ac .* cp, 2);
vc = d1 .* d4 - d3 .* d2;
vb = d5 .* d2 - d1 .* d6;
va = d3 .* d6 - d5 .* d4;

closest = zeros(numel(ids), 3, 'like', a);
remaining = true(numel(ids), 1);
mask = remaining & d1 <= 0 & d2 <= 0;
closest(mask,:) = a(mask,:); remaining(mask) = false;
mask = remaining & d3 >= 0 & d4 <= d3;
closest(mask,:) = b(mask,:); remaining(mask) = false;
mask = remaining & vc <= 0 & d1 >= 0 & d3 <= 0;
assignAB(mask);
mask = remaining & d6 >= 0 & d5 <= d6;
closest(mask,:) = c(mask,:); remaining(mask) = false;
mask = remaining & vb <= 0 & d2 >= 0 & d6 <= 0;
assignAC(mask);
mask = remaining & va <= 0 & (d4-d3) >= 0 & (d5-d6) >= 0;
assignBC(mask);
assignFace(remaining);

distanceSquared = sum(bsxfun(@minus, centre, closest).^2, 2);
intersects = any(distanceSquared < limit^2);

    function assignAB(mask)
        rows = find(mask);
        if isempty(rows), return; end
        weight = d1(rows) ./ (d1(rows) - d3(rows));
        closest(rows,:) = a(rows,:) + bsxfun(@times, weight, ab(rows,:));
        remaining(rows) = false;
    end

    function assignAC(mask)
        rows = find(mask);
        if isempty(rows), return; end
        weight = d2(rows) ./ (d2(rows) - d6(rows));
        closest(rows,:) = a(rows,:) + bsxfun(@times, weight, ac(rows,:));
        remaining(rows) = false;
    end

    function assignBC(mask)
        rows = find(mask);
        if isempty(rows), return; end
        numerator = d4(rows) - d3(rows);
        denominator = numerator + d5(rows) - d6(rows);
        weight = numerator ./ denominator;
        closest(rows,:) = b(rows,:) + bsxfun(@times, weight, c(rows,:) - b(rows,:));
        remaining(rows) = false;
    end

    function assignFace(mask)
        rows = find(mask);
        if isempty(rows), return; end
        denominator = va(rows) + vb(rows) + vc(rows);
        v = vb(rows) ./ denominator;
        w = vc(rows) ./ denominator;
        closest(rows,:) = a(rows,:) + bsxfun(@times, v, ab(rows,:)) + ...
            bsxfun(@times, w, ac(rows,:));
        remaining(rows) = false;
    end
end
