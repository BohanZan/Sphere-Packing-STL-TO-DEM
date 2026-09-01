function distance = spPointTriangleDistance(point, triangle)
%SPPOINTTRIANGLEDISTANCE Minimum distance to a finite triangle.
%Classify the closest feature as a vertex, edge or triangle interior.
a=triangle(1,:); b=triangle(2,:); c=triangle(3,:); ab=b-a; ac=c-a; ap=point-a;

%Test the Voronoi region of vertex A and the adjacent edge AB.
d1=dot(ab,ap); d2=dot(ac,ap);
if d1<=0 && d2<=0, distance=norm(ap); return; end
bp=point-b; d3=dot(ab,bp); d4=dot(ac,bp);
if d3>=0 && d4<=d3, distance=norm(bp); return; end
vc=d1*d4-d3*d2;
if vc<=0 && d1>=0 && d3<=0, distance=norm(point-(a+d1/(d1-d3)*ab)); return; end

%Test the Voronoi region of vertex C and the adjacent edge AC.
cp=point-c; d5=dot(ab,cp); d6=dot(ac,cp);
if d6>=0 && d5<=d6, distance=norm(cp); return; end
vb=d5*d2-d1*d6;
if vb<=0 && d2>=0 && d6<=0, distance=norm(point-(a+d2/(d2-d6)*ac)); return; end

%Test the remaining edge BC before falling through to the face interior.
va=d3*d6-d5*d4;
if va<=0 && d4-d3>=0 && d5-d6>=0
    distance=norm(point-(b+(d4-d3)/((d4-d3)+(d5-d6))*(c-b))); return;
end

%Project the point onto the triangle interior using barycentric coordinates.
v=vb/(va+vb+vc); w=vc/(va+vb+vc);
distance=norm(point-(a+v*ab+w*ac));
end
