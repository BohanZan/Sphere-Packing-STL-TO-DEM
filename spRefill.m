function [state, nextRadius] = spRefill(context, state, radii, nextRadius, options)
% Algorithm 4: place remaining prescribed radii tangent to active faces.
%maxRefillPasses is the cumulative empty-sweep limit (paper Mr), not a cap
%on successful insertion/compression cycles. Keep the legacy option name.
if options.maxRefillPasses < 1
 fprintf('[Refilling] Skipped: disabled (rejection limit < 1).\n');
 return;
end
if nextRadius > numel(radii)
 fprintf('[Refilling] Finished: prescribed radii exhausted; total %d/%d.\n',state.count,numel(radii));
 return;
end
% Public callers already cache this; direct callers may supply columns or
% unnormalised magnitudes too. Keep no-op exits above free of side effects.
if ~isfield(context,'gravityFrame') || ~isequal(context.gravityFrame.direction,options.gravity)
 context.gravityFrame=spGravityFrame(context,options.gravity);
end
options.gravity=context.gravityFrame.direction;
% Ignore roundoff-sized projections of faces perpendicular to gravity.
active = find(context.inwardNormals * options.gravity.' > 64*eps).';
if isempty(active)
 fprintf('[Refilling] Finished: no active faces; total %d/%d.\n',state.count,numel(radii));
 return;
end
rejects=0; round=0; faceCount=numel(active);
while rejects<options.maxRefillPasses && nextRadius<=numel(radii)
 %N0/N1 are the actual newly inserted range across this complete traversal.
 batchStart=state.count+1;
 round=round+1; facesDone=0;
 fprintf('[Refilling] Round %d | placement | faces 0/%d; added 0; total %d/%d\n', ...
  round,faceCount,state.count,numel(radii));
 progressClock=tic;
 for faceIndex=1:faceCount
  if nextRadius>numel(radii), break; end
  id=active(faceIndex);
  tri=context.vertices(context.faces(id,:),:); n=context.inwardNormals(id,:);
  while nextRadius<=numel(radii)
   r=radii(nextRadius); placed=false;
   for attempt=1:options.maxAttempts
    %Check time sparsely in long rejection runs; never print per candidate.
    if (attempt==1 || mod(attempt,64)==0) && toc(progressClock)>=1
     fprintf('[Refilling] Round %d | placement | faces %d/%d; added %d; total %d/%d\n', ...
      round,facesDone,faceCount,state.count-batchStart+1,state.count,numel(radii));
     progressClock=tic;
    end
    q=randomTrianglePoint(tri);
    if spCanPlace(context,state,q+r*n,r)
     state=spAddSphere(context,state,q+r*n,r);
     state=spReportFillProgress(state,numel(radii));
     placed=true;
     break;
    end
   end
   %Failure on this face leaves the same prescribed radius for the next face.
   if ~placed, break; end
   nextRadius=nextRadius+1;
  end
  facesDone=faceIndex;
 end
 fprintf('[Refilling] Round %d | placement done | faces %d/%d; added %d; total %d/%d\n', ...
  round,facesDone,faceCount,state.count-batchStart+1,state.count,numel(radii));
 if state.count<batchStart
  %An empty batch has nothing to compress. Retry until Mr empty traversals.
  rejects=rejects+1;
  fprintf('[Refilling] Round %d | no new spheres; rejections %d/%d\n', ...
   round,rejects,options.maxRefillPasses);
 else
  %Algorithm 2 then 3 for this batch only; old spheres remain obstacles.
  progressPrefix=sprintf('[Refilling] Round %d',round);
  state=spRelax(context,state,options.gravity,options,true,batchStart,progressPrefix);
  state=spRelax(context,state,options.gravity,options,false,batchStart,progressPrefix);
 end
end
if nextRadius>numel(radii)
 reason='prescribed radii exhausted';
else
 reason='rejection limit reached';
end
fprintf('[Refilling] Finished: %s; rounds %d; rejections %d/%d; total %d/%d\n', ...
 reason,round,rejects,options.maxRefillPasses,state.count,numel(radii));
end
function q=randomTrianglePoint(tri)
%RANDOMTRIANGLEPOINT Draw a uniform point from a triangular face.
%Reflect samples across u+v=1 to preserve uniform barycentric density.
u=rand; v=rand; if u+v>1, u=1-u; v=1-v; end
q=tri(1,:)+u*(tri(2,:)-tri(1,:))+v*(tri(3,:)-tri(1,:));
end
