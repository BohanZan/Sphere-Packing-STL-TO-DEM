function intersects=spSphereHitsTriangles(context,centre,radius,triangleIds)
%SPSPHEREHITSTRIANGLES Check finite-triangle penetration in stored face order.
%The shared distance solver also handles repeated and collinear vertices.
intersects=false;
limit=radius-context.tolerance;
if limit<=0, return; end
for id=triangleIds(:).'
 triangle=[context.triangles.a(id,:);context.triangles.b(id,:);context.triangles.c(id,:)];
 if spPointTriangleDistance(centre,triangle)<limit
  intersects=true;
  return;
 end
end
end
