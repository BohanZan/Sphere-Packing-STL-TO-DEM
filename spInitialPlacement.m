function [state, nextRadius] = spInitialPlacement(context, state, radii, nextRadius, options)
% Algorithm 1: random centres, then cuboid, parity, sphere and triangle tests.
while nextRadius<=numel(radii)
 r=radii(nextRadius); placed=false;
 lo=context.lower+r; hi=context.upper-r;
 if any(lo>hi), return; end
 for attempt=1:options.maxAttempts
  p=lo+rand(1,3).*(hi-lo);
  if spCanPlace(context,state,p,r), state=spAddSphere(context,state,p,r); placed=true; break; end
 end
 if ~placed, return; end
 nextRadius=nextRadius+1;
end
end
