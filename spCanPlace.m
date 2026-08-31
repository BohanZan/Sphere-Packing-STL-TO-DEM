function yes = spCanPlace(context, state, centre, radius, ignoreId)
yes=false;
if nargin<5, ignoreId=0; end
if any(centre-radius<context.lower) || any(centre+radius>context.upper), return; end
if ~spGeometry('inside',context,centre), return; end
idx=spGeometry('cell',context,centre);
ids=spGeometry('nearby',context,state.sphereCells,idx);
for id=ids
 if id==ignoreId, continue; end
 if norm(centre-state.centres(id,:)) < radius+state.radii(id)-context.tolerance, return; end
end
triIds=spGeometry('nearby',context,context.triangleCells,idx);
if spGeometry('triangleHit',context,centre,radius,triIds), return; end
yes=true;
end
