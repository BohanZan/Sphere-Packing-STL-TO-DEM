function distance=spPointTriangleDistance(point,triangle)
%SPPOINTTRIANGLEDISTANCE Minimum distance to a finite triangle.
%Use the same vertex, edge and interior classifier as continuous contacts.
closest=spClosestPoint(point,triangle);
distance=norm(double(point(:).')-closest);
end
