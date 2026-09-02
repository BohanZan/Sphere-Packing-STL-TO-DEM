function [state, nextRadius] = spRefill(context, state, radii, nextRadius, options)
% Algorithm 4: place remaining prescribed radii tangent to active faces.
%Use gravity-facing surface triangles to refill the difficult boundary region.
if options.maxRefillPasses < 1 || nextRadius > numel(radii), return; end
active = find(context.inwardNormals * options.gravity.' > 0).';
if isempty(active), active = 1:size(context.faces,1); end
for pass=1:options.maxRefillPasses
 %Record whether this refill sweep found any new feasible sphere centres.
 before=state.count;
 while nextRadius<=numel(radii)
  r=radii(nextRadius); placed=false;
  %Sample points on active triangles and offset them along the inward normal.
  for attempt=1:options.maxAttempts
   id=active(randi(numel(active))); tri=context.vertices(context.faces(id,:),:);
   q=randomTrianglePoint(tri); n=context.inwardNormals(id,:);
   if spCanPlace(context,state,q+r*n,r)
    state=spAddSphere(context,state,q+r*n,r);
    state=spReportFillProgress(state,numel(radii));
    placed=true;
    break;
   end
  end
  %End the current sweep if the unresolved radius cannot be placed.
  if ~placed, break; end
  nextRadius=nextRadius+1;
 end
 %Avoid further relaxation when this pass made no geometric progress.
 if state.count==before, return; end

 %Settle the newly added boundary spheres before the next refill pass.
 state=spRelax(context,state,options.gravity,options,true);
 state=spRelax(context,state,options.gravity,options,false);
 if nextRadius>numel(radii), return; end
end
end
function q=randomTrianglePoint(tri)
%RANDOMTRIANGLEPOINT Draw a uniform point from a triangular face.
%Reflect samples across u+v=1 to preserve uniform barycentric density.
u=rand; v=rand; if u+v>1, u=1-u; v=1-v; end
q=tri(1,:)+u*(tri(2,:)-tri(1,:))+v*(tri(3,:)-tri(1,:));
end
