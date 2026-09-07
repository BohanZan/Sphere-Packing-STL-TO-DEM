function yes = spCanPlace(context, state, centre, radius, ignoreId)
%SPCANPLACE Test one candidate sphere against bounds, spheres and STL faces.
yes=false;
if nargin<5, ignoreId=0; end

%Reject centres whose enclosing sphere leaves the axis-aligned STL bounds.
if any(centre-radius<context.lower) || any(centre+radius>context.upper), return; end

%Check only neighbouring sparse cells for overlapping accepted spheres.
idx=floor((centre-context.lower)./context.cellSize)+1;
idx=min(max(idx,1),context.cellCount);
ids=spSphereNeighbours(context,state.sphereCells,idx);
ids(ids==ignoreId)=[];
if ~isempty(ids)
    delta=state.centres(ids,:)-centre;
    limit=radius+state.radii(ids)-context.tolerance;
    limit=limit(:);
    valid=limit>0;
    if any(sum(delta(valid,:).^2,2)<limit(valid).^2), return; end
end

if isstruct(context.triangleCells) && isfield(context.triangleCells,'centreCoverage') && ...
        context.triangleCells.centreCoverage
    % The complete finite distance band already contains every face that
    % can touch a sphere with its centre in this cell.
    triIds=spGridLookup(context.triangleCells,spCellKeys(idx,context.cellKeySpec));
else
% Legacy AABB indices require all cells intersecting this sphere's AABB.
cellLower=context.lower+(idx-1).*context.cellSize;
lo=idx-(centre-radius<cellLower);
hi=idx+(centre+radius>=cellLower+context.cellSize);
lo=min(max(lo,1),context.cellCount);
hi=min(max(hi,1),context.cellCount);
triIds=spGridLookup(context.triangleCells,spBoxCellKeys(lo,hi,context.cellKeySpec));
triIds=unique(triIds);
end
if spSphereHitsTriangles(context,centre,radius,triIds), return; end

%Use parity after every local geometric rejection has passed.
if ~spPointInside(context,centre), return; end
yes=true;
end
