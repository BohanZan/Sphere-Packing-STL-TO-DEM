function intersects = spSphereHitsTriangles(context, centre, radius, triangleIds)
%SPSPHEREHITSTRIANGLES Check finite-triangle penetration, including edges/vertices.
%Evaluate exact point-to-triangle distances for the local triangle candidates.
intersects = false;
if isempty(triangleIds)
    return;
end

%Stop at the first penetrating face because only feasibility is required.
for id = triangleIds
    triangle = context.vertices(context.faces(id,:), :);
    if spPointTriangleDistance(centre, triangle) < radius - context.tolerance
        intersects = true;
        return;
    end
end
end
