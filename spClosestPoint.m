function q=spClosestPoint(p,triangle)
%SPCLOSESTPOINT Closest finite triangle feature, following the C++ branch order.
%Scale extreme relative coordinates before the Voronoi-region calculations.
if ~isequal(size(triangle),[3 3]) || numel(p)~=3 || ...
 ~isreal(p) || ~isreal(triangle) || any(~isfinite([p(:);triangle(:)]))
 error('SpherePacking:InvalidTriangleInput','Triangle distance inputs must be finite real coordinates.');
end
p=double(p(:).'); triangle=double(triangle);
a=triangle(1,:); b=triangle(2,:); c=triangle(3,:);
ab=b-a; ac=c-a; ap=p-a;
scale=max(abs([ab ac ap]));
if ~isfinite(scale)
 error('SpherePacking:TriangleOverflow','Triangle distance arithmetic overflow.');
end
if scale>0 && (scale<1e-50 || scale>1e50)
 q=a+scale*spClosestPoint(ap/scale,[zeros(1,3);ab/scale;ac/scale]);
 return;
end

%A degenerate facet consists only of finite segments and their endpoints.
if norm(cross(ab,ac))==0
 q=a; best=Inf;
 for edge=1:3
  x=triangle(edge,:); e=triangle(mod(edge,3)+1,:)-x;
  denominator=dot(e,e);
  if denominator>0, x=x+min(1,max(0,dot(p-x,e)/denominator))*e; end
  distanceSquared=dot(p-x,p-x);
  if distanceSquared<best, q=x; best=distanceSquared; end
 end
 return;
end

%Test vertex A, vertex B and edge AB before the remaining Voronoi regions.
d1=dot(ab,ap); d2=dot(ac,ap);
if d1<=0 && d2<=0, q=a; return; end
bp=p-b; d3=dot(ab,bp); d4=dot(ac,bp);
if d3>=0 && d4<=d3, q=b; return; end
vc=d1*d4-d3*d2;
if vc<=0 && d1>=0 && d3<=0, q=a+d1/(d1-d3)*ab; return; end

%Test vertex C and edges AC/BC, then project onto the face interior.
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
