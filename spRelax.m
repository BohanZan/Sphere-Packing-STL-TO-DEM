function state = spRelax(context, state, gravity, options, globalPass)
% Algorithms 2 and 3.  A global pass uses gravity; a shake uses random
% directions perpendicular to gravity before each gravity sweep.
previous=inf;
for sweep=1:options.maxCompressionSweeps
 if ~globalPass
  for shake=1:options.shakeSweeps
   direction=randomPerpendicular(gravity); state=moveSweep(context,state,direction);
  end
 end
 state=moveSweep(context,state,gravity);
 energy=sum((state.centres*gravity.').*state.radii.^3);
 if isfinite(previous) && abs(energy/previous-1)<options.compressionTolerance, return; end
 previous=energy;
end
end
function state=moveSweep(context,state,direction)
for id=1:state.count
 distance=firstContactDistance(context,state,id,direction);
 if distance>0
  state.centres(id,:)=state.centres(id,:)+distance*direction;
  state=spReindex(context,state);
 end
end
end
function d=firstContactDistance(context,state,id,direction)
centre=state.centres(id,:); radius=state.radii(id);
limit=norm(context.upper-context.lower); step=max(radius/3,context.tolerance*10);
low=0; high=step;
while high<=limit && spCanPlace(context,state,centre+high*direction,radius,id)
 low=high; high=high+step;
end
if high>limit, d=low; return; end
for k=1:30
 mid=(low+high)/2;
 if spCanPlace(context,state,centre+mid*direction,radius,id), low=mid; else, high=mid; end
end
d=low;
end
function v=randomPerpendicular(g)
v=randn(1,3); v=v-dot(v,g)*g; if norm(v)<eps, v=[g(2),-g(1),0]; end; v=v/norm(v);
end
