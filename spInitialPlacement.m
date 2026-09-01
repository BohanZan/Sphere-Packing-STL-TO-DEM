function [state, nextRadius] = spInitialPlacement(context, state, radii, nextRadius, options)
% Algorithm 1: random centres, then cuboid, parity, sphere and triangle tests.
%Attempt the ordered radii one at a time until a radius cannot be placed.
while nextRadius<=numel(radii)
 %Define the admissible centre cuboid so the entire sphere remains in bounds.
 r=radii(nextRadius); placed=false;
 lo=context.lower+r; hi=context.upper-r;
 if any(lo>hi), return; end
 %Sample random candidate centres and accept the first geometrically valid one.
 for attempt=1:options.maxAttempts
  p=lo+rand(1,3).*(hi-lo);
  if spCanPlace(context,state,p,r)
   state=spAddSphere(context,state,p,r);
   state=spReportFillProgress(state,numel(radii));
   placed=true;
   break;
  end
 end
 %Stop this phase when the current prescribed radius exhausts its trials.
 if ~placed, return; end
 nextRadius=nextRadius+1;
end
end
