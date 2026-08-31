function [state, nextRadius] = spRefill(context, state, radii, nextRadius, options)
% Algorithm 4: place remaining prescribed radii tangent to active faces.
active=spActiveTriangles(context,options.gravity);
for pass=1:options.maxRefillPasses
 before=state.count;
 while nextRadius<=numel(radii)
  r=radii(nextRadius); placed=false;
  for attempt=1:options.maxAttempts
   id=active(randi(numel(active))); tri=context.vertices(context.faces(id,:),:);
   q=randomTrianglePoint(tri); n=inwardNormal(context,tri);
   if spCanPlace(context,state,q+r*n,r)
    state=spAddSphere(context,state,q+r*n,r); placed=true; break;
   end
  end
  if ~placed, break; end
  nextRadius=nextRadius+1;
 end
 if state.count==before, return; end
 state=spRelax(context,state,options.gravity,options,true);
 state=spRelax(context,state,options.gravity,options,false);
 if nextRadius>numel(radii), return; end
end
end
function ids=spActiveTriangles(context,gravity)
ids=[];
for id=1:size(context.faces,1)
 tri=context.vertices(context.faces(id,:),:); n=inwardNormal(context,tri);
 if dot(n,gravity)>0, ids(end+1)=id; end %#ok<AGROW>
end
if isempty(ids), ids=1:size(context.faces,1); end
end
function q=randomTrianglePoint(tri)
u=rand; v=rand; if u+v>1, u=1-u; v=1-v; end
q=tri(1,:)+u*(tri(2,:)-tri(1,:))+v*(tri(3,:)-tri(1,:));
end
function n=inwardNormal(context,tri)
n=cross(tri(2,:)-tri(1,:),tri(3,:)-tri(1,:)); n=n/norm(n);
probe=mean(tri,1)+n*max(context.tolerance*100,1e-8*context.cellSize);
if ~spGeometry('inside',context,probe), n=-n; end
end
