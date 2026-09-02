function state = spRelax(context, state, gravity, options, globalPass)
% Algorithms 2 and 3.  A global pass uses gravity; a shake uses random
% directions perpendicular to gravity before each gravity sweep.
%Repeat directional sweeps until the weighted height measure converges.
previous=inf;
for sweep=1:options.maxCompressionSweeps
 %Perturb laterally before settling to reduce local geometric locking.
 if ~globalPass
  for shake=1:options.shakeSweeps
   direction=randomPerpendicular(gravity); state=moveSweep(context,state,direction);
  end
 end
 %Settle every sphere along gravity, then check the compression criterion.
 state=moveSweep(context,state,gravity);
 validIds=1:state.count;
 energy=sum((state.centres(validIds,:)*gravity.').*state.radii(validIds).^3);
 if isfinite(previous) && abs(energy/previous-1)<options.compressionTolerance, return; end
 previous=energy;
end
end
function state=moveSweep(context,state,direction)
%MOVESWEEP Translate every sphere to its first admissible contact point.
%Rebuild the hash after movement so later spheres see the updated packing.
for id=1:state.count
 distance=firstContactDistance(context,state,id,direction);
 if distance>0
  state.centres(id,:)=state.centres(id,:)+distance*direction;
  state=spReindex(context,state,id);
 end
end
end
function d=firstContactDistance(context,state,id,direction)
%FIRSTCONTACTDISTANCE Find the largest feasible translation by bracketing.
%Use coarse steps first, then fixed-iteration bisection near contact.
centre=state.centres(id,:); radius=state.radii(id);
limit=norm(context.upper-context.lower); step=max(radius/3,context.tolerance*10);
low=0; high=step;
while high<=limit && spCanPlace(context,state,centre+high*direction,radius,id)
 low=high; high=high+step;
end

%Return the last safe coarse step when no contact exists inside the domain.
if high>limit, d=low; return; end

%Refine the first contact distance without crossing into an invalid state.
for k=1:30
 mid=(low+high)/2;
 if spCanPlace(context,state,centre+mid*direction,radius,id), low=mid; else, high=mid; end
end
d=low;
end
function v=randomPerpendicular(g)
%RANDOMPERPENDICULAR Generate a unit random vector perpendicular to gravity.
v=randn(1,3); v=v-dot(v,g)*g; if norm(v)<eps, v=[g(2),-g(1),0]; end; v=v/norm(v);
end
