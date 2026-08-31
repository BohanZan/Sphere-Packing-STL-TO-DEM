function value = spGeometry(action, context, varargin)
%SPGEOMETRY Local spatial queries used by insertion, compression and refill.
switch action
    case 'cell', value = cellOf(context, varargin{1});
    case 'nearby', value = nearby(context, varargin{1}, varargin{2});
    case 'inside', value = inside(context, varargin{1});
    case 'triangleHit', value = triangleHit(context, varargin{1}, varargin{2}, varargin{3});
    otherwise, error('SpherePacking:UnknownGeometryAction','Unknown geometry action.');
end
end

function idx=cellOf(c,p)
idx=min(max(floor((p-c.lower)/c.cellSize)+1,1),c.cellCount);
end
function ids=nearby(c,cells,idx)
ids=[];
for i=max(1,idx(1)-1):min(c.cellCount(1),idx(1)+1)
 for j=max(1,idx(2)-1):min(c.cellCount(2),idx(2)+1)
  for k=max(1,idx(3)-1):min(c.cellCount(3),idx(3)+1)
   ids=[ids,cells{i,j,k}]; %#ok<AGROW>
  end
 end
end
ids=unique(ids);
end
function answer=inside(c,p)
idx=min(max(floor((p(1:2)-c.lower(1:2))/c.xySize)+1,1),c.xyCount);
ids=c.xyCells{idx(1),idx(2)}; z=[];
for id=ids
 tri=c.vertices(c.faces(id,:),:); [hit,h]=verticalHit(p,tri,c.tolerance);
 if hit && h<p(3)-c.tolerance, z(end+1)=h; end %#ok<AGROW>
end
if isempty(z), answer=false; return; end
z=sort(z); answer=mod(1+sum(diff(z)>c.tolerance),2)==1;
end
function [hit,z]=verticalHit(p,t,tol)
a=t(1,1:2); u=t(2,1:2)-a; v=t(3,1:2)-a; q=p(1:2)-a;
d=u(1)*v(2)-u(2)*v(1);
if abs(d)<=tol, hit=false; z=NaN; return; end
x=(q(1)*v(2)-q(2)*v(1))/d; y=(u(1)*q(2)-u(2)*q(1))/d;
hit=x>=-tol && y>=-tol && x+y<=1+tol;
z=t(1,3)+x*(t(2,3)-t(1,3))+y*(t(3,3)-t(1,3));
end
function answer=triangleHit(c,p,r,ids)
answer=false;
for id=ids
 if pointTriDistance(p,c.vertices(c.faces(id,:),:)) < r-c.tolerance
  answer=true; return;
 end
end
end
function d=pointTriDistance(p,t)
% Ericson closest-point construction: face, edges and vertices all included.
a=t(1,:); b=t(2,:); cc=t(3,:); ab=b-a; ac=cc-a; ap=p-a; d1=dot(ab,ap); d2=dot(ac,ap);
if d1<=0 && d2<=0, d=norm(ap); return; end
bp=p-b; d3=dot(ab,bp); d4=dot(ac,bp);
if d3>=0 && d4<=d3, d=norm(bp); return; end
vc=d1*d4-d3*d2; if vc<=0 && d1>=0 && d3<=0, q=a+d1/(d1-d3)*ab; d=norm(p-q); return; end
cp=p-cc; d5=dot(ab,cp); d6=dot(ac,cp);
if d6>=0 && d5<=d6, d=norm(cp); return; end
vb=d5*d2-d1*d6; if vb<=0 && d2>=0 && d6<=0, q=a+d2/(d2-d6)*ac; d=norm(p-q); return; end
va=d3*d6-d5*d4;
if va<=0 && d4-d3>=0 && d5-d6>=0, q=b+(d4-d3)/((d4-d3)+(d5-d6))*(cc-b); d=norm(p-q); return; end
v=vb/(va+vb+vc); w=vc/(va+vb+vc); d=norm(p-(a+v*ab+w*ac));
end
