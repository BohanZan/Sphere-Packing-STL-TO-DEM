function intersects = spReferenceSphereHitsTriangles(context, centre, radius, triangleIds)
%SPREFERENCESPHEREHITSTRIANGLES Independent scalar pre-Task-6 collision oracle.
intersects = false;
if isempty(triangleIds)
    return
end
for id = triangleIds
    triangle = context.vertices(context.faces(id,:), :);
    if spReferencePointTriangleDistance(centre, triangle) < radius - context.tolerance
        intersects = true;
        return
    end
end
end

function distance = spReferencePointTriangleDistance(point, triangle)
%SPREFERENCEPOINTTRIANGLEDISTANCE Literal scalar finite-triangle distance oracle.
a=triangle(1,:); b=triangle(2,:); c=triangle(3,:); ab=b-a; ac=c-a; ap=point-a;
d1=dot(ab,ap); d2=dot(ac,ap);
if d1<=0 && d2<=0, distance=norm(ap); return; end
bp=point-b; d3=dot(ab,bp); d4=dot(ac,bp);
if d3>=0 && d4<=d3, distance=norm(bp); return; end
vc=d1*d4-d3*d2;
if vc<=0 && d1>=0 && d3<=0, distance=norm(point-(a+d1/(d1-d3)*ab)); return; end
cp=point-c; d5=dot(ab,cp); d6=dot(ac,cp);
if d6>=0 && d5<=d6, distance=norm(cp); return; end
vb=d5*d2-d1*d6;
if vb<=0 && d2>=0 && d6<=0, distance=norm(point-(a+d2/(d2-d6)*ac)); return; end
va=d3*d6-d5*d4;
if va<=0 && d4-d3>=0 && d5-d6>=0
    distance=norm(point-(b+(d4-d3)/((d4-d3)+(d5-d6))*(c-b))); return;
end
v=vb/(va+vb+vc); w=vc/(va+vb+vc);
distance=norm(point-(a+v*ab+w*ac));
end
