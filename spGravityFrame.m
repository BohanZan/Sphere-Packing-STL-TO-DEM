function frame=spGravityFrame(context,gravity)
%SPGRAVITYFRAME Relative height along -gravity, without rotating the geometry.
if ~(isnumeric(gravity) && isreal(gravity) && numel(gravity)==3 && ...
        all(isfinite(gravity(:))) && any(gravity(:)~=0))
    error('SpherePacking:InvalidGravity', ...
        'gravity must be a finite, real, nonzero three-element vector.');
end
% Scaling first avoids overflow/underflow for extreme magnitudes.
g=double(gravity(:).'); g=g/max(abs(g)); g=g/norm(g);
up=-g;
%Traverse A/B/C of each face in mesh order, as the C++ triangle array does.
%Unused indexed vertices must not change the domain support or its origin.
vertices=context.vertices;
if isfield(context,'faces')
 vertices=vertices(reshape(context.faces.',[],1),:);
end
if isempty(vertices), error('SpherePacking:EmptyMesh','Mesh must be nonempty.'); end
anchor=vertices(1,:);
bottom=0; top=0; origin=anchor;
for id=1:size(vertices,1)
 height=dot(vertices(id,:)-anchor,up);
 if height<bottom, bottom=height; origin=vertices(id,:); end
 top=max(top,height);
end
frame=struct('direction',g,'up',up,'origin',origin,'height',top-bottom);
end
