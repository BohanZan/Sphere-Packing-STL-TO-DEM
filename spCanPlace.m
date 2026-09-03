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

%Collect only cells intersecting this sphere's AABB (at most 2-by-2-by-2).
cellLower=context.lower+(idx-1).*context.cellSize;
lo=idx-(centre-radius<cellLower);
hi=idx+(centre+radius>=cellLower+context.cellSize);
lo=min(max(lo,1),context.cellCount);
hi=min(max(hi,1),context.cellCount);
parts=cell(8,1);
partCount=0;
triangleCells=context.triangleCells;
for ix=lo(1):hi(1)
    for iy=lo(2):hi(2)
        for iz=lo(3):hi(3)
            key=sprintf('%d,%d,%d',ix,iy,iz);
            if isKey(triangleCells,key)
                partCount=partCount+1;
                parts{partCount}=triangleCells(key);
            end
        end
    end
end
if partCount==0, triIds=[]; else, triIds=unique([parts{1:partCount}]); end
if spSphereHitsTriangles(context,centre,radius,triIds), return; end

%Use parity after every local geometric rejection has passed.
if ~spPointInside(context,centre), return; end
yes=true;
end
