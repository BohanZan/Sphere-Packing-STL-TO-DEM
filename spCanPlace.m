function yes = spCanPlace(context, state, centre, radius, ignoreId)
%SPCANPLACE Test one candidate sphere against bounds, spheres and STL faces.
yes=false;
if nargin<5, ignoreId=0; end

%Reject centres whose enclosing sphere leaves the axis-aligned STL bounds.
if any(centre-radius<context.lower) || any(centre+radius>context.upper), return; end

%Use parity to reject centres outside the closed surface.
if ~spPointInside(context,centre), return; end

%Check only neighbouring sparse cells for overlapping accepted spheres.
idx=spCellIndex(context,centre);
ids=spHashNeighbours(context,state.sphereCells,idx);
for id=ids
 if id==ignoreId, continue; end
 if norm(centre-state.centres(id,:)) < radius+state.radii(id)-context.tolerance, return; end
end

%Check nearby finite triangles so the sphere cannot cross the STL surface.
triIds=spHashNeighbours(context,context.triangleCells,idx);
if spSphereHitsTriangles(context,centre,radius,triIds), return; end
yes=true;
end
