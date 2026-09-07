function distance=spTriangleContactDistance(p,r,u,triangle,limit,tolerance)
%SPTRIANGLECONTACTDISTANCE Continuous first contact with a finite triangle.
% Face entry is the fast branch (paper Eq.14). Only a missed face enters
% the finite edge/vertex branches of Eq.17--21. U is a unit direction.
distance=inf;
if limit<0, return; end
% Work relative to a vertex to avoid dependence on the coordinate origin.
p=p-triangle(1,:); triangle=triangle-triangle(1,:);
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
  q=closestPoint(p,triangle);
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
 q=closestPoint(p,triangle);
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
componentError=16*eps*(abs(offset)+abs(projection)) ...
 +4*norm(offsetError)+4*norm(offset)*norm(velocityError)/sqrt(A);
errorBound=sum(2*abs(perpendicular).*componentError+componentError.^2) ...
 +32*eps(max(r*r,perpendicularSquared));
if discriminant < -errorBound, return; end
root=C/(-B+sqrt(A)*sqrt(max(0,discriminant)));
if root<=hi-lo, distance=lo+root; end
end

function q=closestPoint(p,tri)
% Initial-contact classification only; later hits use analytic branches.
a=tri(1,:); b=tri(2,:); c=tri(3,:); ab=b-a; ac=c-a; ap=p-a;
if dot(cross(ab,ac),cross(ab,ac))==0
 q=a; best=inf;
 for edge=1:3
  x=tri(edge,:); e=tri(mod(edge,3)+1,:)-x; denominator=dot(e,e);
  if denominator>0, x=x+min(1,max(0,dot(p-x,e)/denominator))*e; end
  d=dot(p-x,p-x); if d<best, q=x; best=d; end
 end
 return;
end
d1=dot(ab,ap); d2=dot(ac,ap);
if d1<=0 && d2<=0, q=a; return; end
bp=p-b; d3=dot(ab,bp); d4=dot(ac,bp);
if d3>=0 && d4<=d3, q=b; return; end
vc=d1*d4-d3*d2;
if vc<=0 && d1>=0 && d3<=0, q=a+d1/(d1-d3)*ab; return; end
cp=p-c; d5=dot(ab,cp); d6=dot(ac,cp);
if d6>=0 && d5<=d6, q=c; return; end
vb=d5*d2-d1*d6;
if vb<=0 && d2>=0 && d6<=0, q=a+d2/(d2-d6)*ac; return; end
va=d3*d6-d5*d4;
if va<=0 && d4-d3>=0 && d5-d6>=0
 q=b+(d4-d3)/((d4-d3)+(d5-d6))*(c-b); return;
end
q=a+vb/(va+vb+vc)*ab+vc/(va+vb+vc)*ac;
end
