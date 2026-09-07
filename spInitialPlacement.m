function [state, nextRadius] = spInitialPlacement(context, state, radii, nextRadius, options, frame, layerBase)
%SPINITIALPLACEMENT Algorithm 1: generate the next ordered insertion batch.
%The optional frame restricts centres to one h-high slab opposite gravity.
%Five-input direct calls retain the C++ initialPlacement full-cuboid helper.
if nargin<6, frame=[]; end
if nargin<7, layerBase=0; end
if ~isempty(frame), [~,heightAxis]=max(abs(frame.up)); end
%Attempt the ordered radii one at a time until a radius cannot be placed.
while nextRadius<=numel(radii)
 %Define the admissible centre cuboid so the entire sphere remains in bounds.
 r=radii(nextRadius); placed=false;
 lo=context.lower+r; hi=context.upper-r;
 if any(lo>hi), return; end
 %Sample random candidate centres and accept the first geometrically valid one.
 for attempt=1:options.maxAttempts
  draw=rand(1,3);
  p=lo+draw.*(hi-lo);
  if ~isempty(frame)
   %Intersect the dominant-axis line with both the slab and centre cuboid.
   %Reuse the same three draws even when that line misses the current slab.
   otherHeight=0;
   for axis=1:3
    if axis~=heightAxis
     otherHeight=otherHeight+(p(axis)-frame.origin(axis))*frame.up(axis);
    end
   end
   a=frame.origin(heightAxis)+(layerBase-otherHeight)/frame.up(heightAxis);
   b=frame.origin(heightAxis)+(min(frame.height,layerBase+context.cellSize)-otherHeight)/frame.up(heightAxis);
   lower=max(lo(heightAxis),min(a,b)); upper=min(hi(heightAxis),max(a,b));
   if lower>upper, continue; end
   p(heightAxis)=lower+draw(heightAxis)*(upper-lower);
  end
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
