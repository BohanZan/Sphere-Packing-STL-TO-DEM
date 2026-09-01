function yes = spCanPlace(context, state, centre, radius, ignoreId)
yes=false;
if nargin<5, ignoreId=0; end
if any(centre-radius<context.lower) || any(centre+radius>context.upper), return; end
if ~spPointInside(context,centre), return; end
idx=spCellIndex(context,centre);
ids=spHashNeighbours(context,state.sphereCells,idx);
for id=ids
 if id==ignoreId, continue; end
 if norm(centre-state.centres(id,:)) < radius+state.radii(id)-context.tolerance, return; end
end
triIds=spHashNeighbours(context,context.triangleCells,idx);
if spSphereHitsTriangles(context,centre,radius,triIds), return; end
yes=true;
end
