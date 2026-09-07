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
vertices=context.vertices;
anchor=vertices(1,:);
heights=(vertices-anchor)*up.';
[bottom,bottomId]=min(heights);
frame=struct('direction',g,'up',up,'origin',vertices(bottomId,:), ...
    'height',max(heights)-bottom);
end
