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
ids=spSphereNeighbours(context,state.sphereCells,idx);
ids(ids==ignoreId)=[];
if ~isempty(ids)
    delta=state.centres(ids,:)-centre;
    limit=radius+state.radii(ids)-context.tolerance;
    limit=limit(:);
    valid=limit>0;
    if any(sum(delta(valid,:).^2,2)<limit(valid).^2), return; end
end

%Check nearby finite triangles so the sphere cannot cross the STL surface.
triIds=spHashNeighbours(context,context.triangleCells,idx);
if spSphereHitsTriangles(context,centre,radius,triIds), return; end
yes=true;
end
