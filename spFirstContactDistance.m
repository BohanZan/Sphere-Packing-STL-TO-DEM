function distance=spFirstContactDistance(context,state,id,direction,maxDistance,halo)
%SPFIRSTCONTACTDISTANCE Continuous, grid-accelerated translation to contact.
% Traverse every centre cell along the path (DDA), not just its start/end.
%With cellSize >= 2*maxRadius, a one-cell halo covers stationary spheres.
%Relaxation supplies a wider halo when sweep-start cells require it.
%Triangle primitives are visited once; sphere neighbours retain bucket order.
% Optional nonnegative MAXDISTANCE bounds the geometric query. Compression
%and shaking supply b_u (paper Eq.22); four inputs query full contact.
if nargin<5, maxDistance=inf; end
if nargin<6, halo=1; end
if numel(direction)~=3 || ~isreal(direction) || any(~isfinite(direction(:))) || abs(norm(direction)-1)>1e-10
 error('SpherePacking:InvalidDirection','Movement direction must have unit length.');
end
direction=direction(:).';
p=state.centres(id,:); r=state.radii(id); u=direction;
positive=u>0; negative=u<0;
limits=inf(1,3);
limits(positive)=(context.upper(positive)-r-p(positive))./u(positive);
limits(negative)=(context.lower(negative)+r-p(negative))./u(negative);
distance=min(maxDistance,max(0,min(limits)));
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
 [sphereIds,cellKeys]=neighbours(context,state.sphereCells,index,halo);
 ids=sphereIds(sphereIds~=id);
 for j=ids
  delta=state.centres(j,:)-p;
  along=dot(delta,u); radii=r+state.radii(j);
  if dot(delta,delta)<=(radii+context.tolerance)*(radii+context.tolerance)
   roundoff=64*eps*sum(abs(delta.*u));
   if along>roundoff, distance=0; return; end
   continue;
  end
  if along<=0, continue; end
  projection=along*u;
  separation=delta-projection;
  discriminant=radii*radii-dot(separation,separation);
  % Cancellation error grows with travel length, not just sphere radius.
  errorBound=32*eps*radii*radii;
  for axis=1:3
   componentError=16*eps*(abs(delta(axis))+abs(projection(axis)));
   errorBound=errorBound+(2*abs(separation(axis))*componentError+componentError*componentError);
  end
  if discriminant>=-errorBound
   distance=min(distance,max(0,along-sqrt(max(0,discriminant))));
  end
  if distance==0, return; end
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
  for j=ids
   %Each earlier hit shortens the swept box used for the next finite face.
   finish=p+distance*u;
   lower=min(p,finish)-r-context.tolerance;
   upper=max(p,finish)+r+context.tolerance;
   tri=[context.triangles.a(j,:);context.triangles.b(j,:);context.triangles.c(j,:)];
   if any(min(tri,[],1)>upper | max(tri,[],1)<lower), continue; end
   hit=spTriangleContactDistance(p,r,u,tri,distance,context.tolerance);
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

function [sphereIds,cellKeys]=neighbours(context,sphereCells,index,halo)
lo=max(index-halo,1); hi=min(index+halo,context.cellCount);
cellKeys=spBoxCellKeys(lo,hi,context.cellKeySpec);
present=isKey(sphereCells,cellKeys);
if any(present)
 parts=values(sphereCells,cellKeys(present)); sphereIds=[parts{:}];
else
 sphereIds=[];
end
end
