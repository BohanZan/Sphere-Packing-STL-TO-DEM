function distance=spFirstContactDistance(context,state,id,direction)
%SPFIRSTCONTACTDISTANCE Continuous, grid-accelerated translation to contact.
% Traverse every centre cell along the path (DDA), not just its start/end.
% With cellSize >= 2*maxRadius, each cell's 27 neighbours cover all sphere
% and triangle contacts within that cell. Triangle primitives are memoized;
% sphere pairs use the cheap vectorized quadratic without set operations.
% BUFFER retains its grid-padding meaning; zero does not disable movement.
p=state.centres(id,:); r=state.radii(id); u=direction;
positive=u>0; negative=u<0;
limits=inf(1,3);
limits(positive)=(context.upper(positive)-r-p(positive))./u(positive);
limits(negative)=(context.lower(negative)+r-p(negative))./u(negative);
distance=max(0,min(limits));
if distance==0 || ~isfinite(distance), return; end
h=context.cellSize;
index=spCellIndex(context,p);
steps=sign(u);
crossing=inf(1,3); increments=inf(1,3);
for axis=1:3
 if u(axis)>0
  boundary=context.lower(axis)+index(axis)*h;
 elseif u(axis)<0
  boundary=context.lower(axis)+(index(axis)-1)*h;
 else
  continue;
 end
 crossing(axis)=max(0,(boundary-p(axis))/u(axis));
 increments(axis)=h/abs(u(axis));
end
seenTriangles=[];
centreCoverage=isstruct(context.triangleCells) && ...
 isfield(context.triangleCells,'centreCoverage') && context.triangleCells.centreCoverage;
while true
 % Build the same neighbour keys once for both hashes. ID membership is
 % enough here: setdiff's additional sort/unique work is unnecessary.
 [sphereIds,cellKeys]=neighbours(context,state.sphereCells,index);
 ids=sphereIds(sphereIds~=id);
 if ~isempty(ids)
  delta=state.centres(ids,:)-p;
  along=delta*u.';
  radii=r+state.radii(ids);
  contact=sum(delta.^2,2)<=(radii+context.tolerance).^2;
  roundoff=64*eps*sum(abs(delta.*u),2);
  if any(contact & along>roundoff), distance=0; return; end
  keep=~contact & along>0;
  if any(keep)
   projection=along(keep)*u;
   separation=delta(keep,:)-projection;
   discriminant=radii(keep).^2-sum(separation.^2,2);
   % Cancellation error grows with travel length, not just sphere radius.
   componentError=16*eps*(abs(delta(keep,:))+abs(projection));
   errorBound=sum(2*abs(separation).*componentError+componentError.^2,2) ...
    +32*eps(radii(keep).^2);
   valid=discriminant>=-errorBound;
   if any(valid)
    lengths=along(keep); lengths=lengths(valid);
    hit=lengths-sqrt(max(0,discriminant(valid)));
    distance=min(distance,max(0,min(hit)));
   end
  end
 end
 % A sphere contact may already have stopped the move; defer the triangle
 % hash lookup until needed, reusing the neighbour keys above.
 if centreCoverage
  % The band covers all contacts before the next DDA crossing. Every
  % encountered face is still solved along the complete remaining path.
  ids=spGridLookup(context.triangleCells,spCellKeys(index,context.cellKeySpec));
 else
  ids=spGridLookup(context.triangleCells,cellKeys);
 end
 if ~isempty(ids)
  if ~centreCoverage
   ids=sort(ids);
   ids=ids([true,diff(ids)~=0]);
  end
  if ~isempty(seenTriangles), ids=ids(~ismember(ids,seenTriangles)); end
 else
  ids=[];
 end
 if ~isempty(ids)
  seenTriangles=[seenTriangles,ids]; %#ok<AGROW>
  % Cheap swept-AABB rejection precedes finite-feature calculations.
  finish=p+distance*u;
  lower=min(p,finish)-r-context.tolerance;
  upper=max(p,finish)+r+context.tolerance;
  a=context.triangles.a(ids,:); b=context.triangles.b(ids,:); c=context.triangles.c(ids,:);
  keep=all(min(min(a,b),c)<=upper & max(max(a,b),c)>=lower,2);
  for k=find(keep).'
   hit=spTriangleContactDistance(p,r,u,[a(k,:);b(k,:);c(k,:)],distance,context.tolerance);
   distance=min(distance,hit);
   if distance==0, return; end
  end
 end
 nextCrossing=min(crossing);
 if distance<=nextCrossing, return; end
 axes=crossing==nextCrossing;
 index(axes)=index(axes)+steps(axes);
 if any(index<1 | index>context.cellCount), return; end
 crossing(axes)=crossing(axes)+increments(axes);
end
end

function [sphereIds,cellKeys]=neighbours(context,sphereCells,index)
lo=max(index-1,1); hi=min(index+1,context.cellCount);
cellKeys=spBoxCellKeys(lo,hi,context.cellKeySpec);
present=isKey(sphereCells,cellKeys);
if any(present)
 parts=values(sphereCells,cellKeys(present)); sphereIds=[parts{:}];
else
 sphereIds=[];
end
end
