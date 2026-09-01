function intersects = spSphereHitsTriangles(context, centre, radius, triangleIds)
%SPSPHEREHITSTRIANGLES Check finite-triangle penetration, including edges/vertices.
intersects = false;
for id = triangleIds
    triangle = context.vertices(context.faces(id,:), :);
    if spPointTriangleDistance(centre, triangle) < radius - context.tolerance
        intersects = true;
        return;
    end
end
end
