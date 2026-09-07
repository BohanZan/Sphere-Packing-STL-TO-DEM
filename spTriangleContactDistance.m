function distance=spTriangleContactDistance(p,r,u,triangle,limit,tolerance)
%SPTRIANGLECONTACTDISTANCE Continuous first contact with a finite triangle.
% Face entry is the fast branch (paper Eq.14). Only a missed face enters
% the finite edge/vertex branches of Eq.17--21. U is a unit direction.
distance=inf;
if numel(p)~=3 || numel(u)~=3 || ~isequal(size(triangle),[3 3]) || ...
 ~isreal([p(:);u(:);triangle(:);r;limit;tolerance]) || ...
 any(~isfinite([p(:);u(:);triangle(:);r;tolerance])) || r<0 || isnan(limit) || tolerance<0
 error('SpherePacking:InvalidTriangleInput','Invalid triangle contact input.');
end
if limit<0, return; end
%C++ receives doubles even when the caller's source mesh used float storage.
p=double(p(:).'); u=double(u(:).'); triangle=double(triangle);
r=double(r); limit=double(limit); tolerance=double(tolerance);
% Work relative to a vertex to avoid dependence on the coordinate origin.
p=p-triangle(1,:); triangle=triangle-triangle(1,:);
scale=max([r,abs(p),abs(triangle(2,:)),abs(triangle(3,:))]);
if ~isfinite(scale), error('SpherePacking:TriangleOverflow','Triangle contact arithmetic overflow.'); end
if scale>0 && (scale<1e-50 || scale>1e50)
 distance=scale*spTriangleContactDistance(p/scale,r/scale,u,triangle/scale,limit/scale,tolerance/scale);
 return;
end
ab=triangle(2,:); ac=triangle(3,:);
normal=cross(ab,ac); normalLength=norm(normal);
lo=0; hi=limit;
if normalLength>0
 n=normal/normalLength;
 h=dot(p,n); velocity=dot(u,n);
 if velocity==0
  if abs(h)>r+tolerance, return; end
 else
  times=[(-r-h)/velocity,(r-h)/velocity];
  lo=max(0,min(times)); hi=min(limit,max(times));
  if lo>hi, return; end
 end
 % A settled particle may slide on or separate from this convex obstacle.
 % Blocking every zero root here would incorrectly disable local shaking.
 if abs(h)<=r+tolerance
  q=spClosestPoint(p,triangle);
  delta=p-q;
  if norm(delta)<=r+tolerance
   roundoff=64*eps*sum(abs(delta.*u));
   if dot(delta,u)<-roundoff, distance=0; end
   return;
  end
 end
 if lo>0
  atEntry=p+lo*u;
  alpha=dot(cross(atEntry,ac),n)/normalLength;
  beta=dot(cross(ab,atEntry),n)/normalLength;
  if alpha>=0 && beta>=0 && alpha+beta<=1
   distance=lo;
   return; % No edge/vertex can precede a valid plane-slab entry.
  end
 end
else
 % Degenerate facets reduce to their finite edges/vertices.
 q=spClosestPoint(p,triangle);
 delta=p-q;
 if norm(delta)<=r+tolerance
  if dot(delta,u)<-64*eps*sum(abs(delta.*u)), distance=0; end
  return;
 end
end

% Each edge is a capsule. Its closest-feature parameter s(t)=s0+s1*t
% selects one quadratic on each nonempty time interval, not all features.
for edge=1:3
 a=triangle(edge,:); b=triangle(mod(edge,3)+1,:);
 hit=edgeContact(p-a,u,b-a,r,lo,min(hi,distance));
 distance=min(distance,hit);
end
end

function distance=edgeContact(offset,u,edge,r,lo,hi)
distance=inf;
if lo>hi, return; end
lengthSquared=dot(edge,edge);
if lengthSquared==0
 distance=entryRoot(offset,u,r,lo,hi);
 return;
end
s0=dot(offset,edge)/lengthSquared;
s1=dot(u,edge)/lengthSquared;
% Partition time at s=0 and s=1, in chronological order.
cuts=[lo hi];
if s1~=0
 events=[-s0/s1,(1-s0)/s1];
 cuts=sort([cuts,events(events>lo & events<hi)]);
end
for k=1:numel(cuts)-1
 start=cuts(k); finish=cuts(k+1);
 s=s0+s1*(start+(finish-start)/2);
 if s<0
  hit=entryRoot(offset,u,r,start,finish);
 elseif s>1
  hit=entryRoot(offset-edge,u,r,start,finish);
 else
  offsetProjection=s0*edge; velocityProjection=s1*edge;
  % Projection can cancel long along-edge coordinates before root solving.
  offsetError=16*eps*(abs(offset)+abs(offsetProjection));
  velocityError=16*eps*(abs(u)+abs(velocityProjection));
  hit=entryRoot(offset-offsetProjection,u-velocityProjection,r,start,finish, ...
   offsetError,velocityError);
 end
 if isfinite(hit), distance=hit; return; end
end
end

function distance=entryRoot(offset,velocity,r,lo,hi,offsetError,velocityError)
% Solve ||offset+t*velocity||^2=r^2 in [lo,hi], with a stable entry root.
distance=inf;
if nargin<7, offsetError=zeros(1,3); velocityError=zeros(1,3); end
offsetError=offsetError+abs(lo)*velocityError+16*eps*(abs(offset)+abs(lo*velocity));
offset=offset+lo*velocity;
A=dot(velocity,velocity); B=dot(offset,velocity);
C=dot(offset,offset)-r*r;
if C<=0, distance=lo; return; end
if A==0 || B>=0, return; end
% Avoid B^2-A*C: long nearly tangent paths subtract two huge squares.
% Closest-approach residuals retain the small perpendicular gap instead.
projection=(B/A)*velocity;
perpendicular=offset-projection;
perpendicularSquared=dot(perpendicular,perpendicular);
discriminant=r*r-perpendicularSquared;
extra=4*norm(offsetError)+4*norm(offset)*norm(velocityError)/sqrt(A);
componentError=16*eps*(abs(offset)+abs(projection))+extra;
errorBound=dot(2*abs(perpendicular),componentError)+dot(componentError,componentError) ...
 +32*eps(max(r*r,perpendicularSquared));
if discriminant < -errorBound, return; end
root=C/(-B+sqrt(A)*sqrt(max(0,discriminant)));
if root<=hi-lo, distance=lo+root; end
end
