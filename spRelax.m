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
contextLower=context.lower; contextUpper=context.upper;
cellSize=context.cellSize; cellCount=context.cellCount; tolerance=context.tolerance;
sphereCells=state.sphereCells; stateCentres=state.centres; stateRadii=state.radii;
triangleCells=context.triangleCells;
centre=stateCentres(id,:); radius=stateRadii(id);
lowerCentre=contextLower+radius; upperCentre=contextUpper-radius;
limit=norm(contextUpper-contextLower); step=max(radius/3,tolerance*10);
low=0; high=step;
cachedSphereIndex=[0 0 0]; cachedSphereIds=[];
cachedTriangleLo=[0 0 0]; cachedTriangleHi=[0 0 0]; cachedTriangleIds=[];
cachedTriangleA=[]; cachedTriangleB=[]; cachedTriangleC=[];
cachedTriangleAB=[]; cachedTriangleAC=[];
triangleLimitSquared=(radius-tolerance)^2;
cachedPointXY=[NaN NaN]; cachedRayHeights=[];
while high<=limit
 candidate=centre+high*direction;
 if any(candidate<lowerCentre | candidate>upperCentre), break; end
 if ~canPlaceCached(candidate), break; end
 low=high; high=high+step;
end

%Return the last safe coarse step when no contact exists inside the domain.
if high>limit, d=low; return; end

%Refine the first contact distance without crossing into an invalid state.
for k=1:30
 mid=(low+high)/2;
 candidate=centre+mid*direction;
 if any(candidate<lowerCentre | candidate>upperCentre)
  high=mid;
 elseif canPlaceCached(candidate)
  low=mid;
 else
  high=mid;
 end
end
d=low;

 function yes=canPlaceCached(candidate)
 %Reuse hash queries while candidate centres remain in the same grid cells.
 yes=false;
 if any(candidate-radius<contextLower) || any(candidate+radius>contextUpper), return; end

 sphereIndex=floor((candidate-contextLower)./cellSize)+1;
 sphereIndex=min(max(sphereIndex,1),cellCount);
 if ~isequal(sphereIndex,cachedSphereIndex)
  cachedSphereIndex=sphereIndex;
  cachedSphereIds=spSphereNeighbours(context,sphereCells,sphereIndex);
  cachedSphereIds(cachedSphereIds==id)=[];
 end
 if ~isempty(cachedSphereIds)
  delta=stateCentres(cachedSphereIds,:)-candidate;
  overlapLimit=radius+stateRadii(cachedSphereIds)-tolerance;
  overlapLimit=overlapLimit(:);
  valid=overlapLimit>0;
  if any(sum(delta(valid,:).^2,2)<overlapLimit(valid).^2), return; end
 end

 cellLower=contextLower+(sphereIndex-1).*cellSize;
 triangleLo=sphereIndex-(candidate-radius<cellLower);
 triangleHi=sphereIndex+(candidate+radius>=cellLower+cellSize);
 triangleLo=min(max(triangleLo,1),cellCount);
 triangleHi=min(max(triangleHi,1),cellCount);
 if ~isequal(triangleLo,cachedTriangleLo) || ~isequal(triangleHi,cachedTriangleHi)
  cachedTriangleLo=triangleLo; cachedTriangleHi=triangleHi;
  parts=cell(8,1); partCount=0;
  for ix=triangleLo(1):triangleHi(1)
   for iy=triangleLo(2):triangleHi(2)
    for iz=triangleLo(3):triangleHi(3)
     key=sprintf('%d,%d,%d',ix,iy,iz);
     if isKey(triangleCells,key)
      partCount=partCount+1;
      parts{partCount}=triangleCells(key);
     end
    end
   end
  end
  if partCount==0
   cachedTriangleIds=[];
   cachedTriangleA=[]; cachedTriangleB=[]; cachedTriangleC=[];
   cachedTriangleAB=[]; cachedTriangleAC=[];
  else
   cachedTriangleIds=unique([parts{1:partCount}]);
   cachedTriangleA=context.triangles.a(cachedTriangleIds,:);
   cachedTriangleB=context.triangles.b(cachedTriangleIds,:);
   cachedTriangleC=context.triangles.c(cachedTriangleIds,:);
   cachedTriangleAB=cachedTriangleB-cachedTriangleA;
   cachedTriangleAC=cachedTriangleC-cachedTriangleA;
  end
 end
 if ~isempty(cachedTriangleIds) && sphereHitsCached(candidate), return; end
 if ~pointInsideCached(candidate), return; end
 yes=true;
 end

 function intersects=sphereHitsCached(candidate)
 %Reuse triangle vertices and edge vectors for every candidate in one cell range.
 a=cachedTriangleA; b=cachedTriangleB; c=cachedTriangleC;
 ab=cachedTriangleAB; ac=cachedTriangleAC;
 ap=candidate-a;
 bp=candidate-b;
 cp=candidate-c;
 d1=sum(ab.*ap,2); d2=sum(ac.*ap,2);
 d3=sum(ab.*bp,2); d4=sum(ac.*bp,2);
 d5=sum(ab.*cp,2); d6=sum(ac.*cp,2);
 vc=d1.*d4-d3.*d2;
 vb=d5.*d2-d1.*d6;
 va=d3.*d6-d5.*d4;

 triangleCount=size(a,1);
 closest=zeros(triangleCount,3,'like',a);
 remaining=true(triangleCount,1);
 mask=remaining & d1<=0 & d2<=0;
 closest(mask,:)=a(mask,:); remaining(mask)=false;
 mask=remaining & d3>=0 & d4<=d3;
 closest(mask,:)=b(mask,:); remaining(mask)=false;
 mask=remaining & vc<=0 & d1>=0 & d3<=0;
 if any(mask)
  weight=d1(mask)./(d1(mask)-d3(mask));
  closest(mask,:)=a(mask,:)+weight.*ab(mask,:);
  remaining(mask)=false;
 end
 mask=remaining & d6>=0 & d5<=d6;
 closest(mask,:)=c(mask,:); remaining(mask)=false;
 mask=remaining & vb<=0 & d2>=0 & d6<=0;
 if any(mask)
  weight=d2(mask)./(d2(mask)-d6(mask));
  closest(mask,:)=a(mask,:)+weight.*ac(mask,:);
  remaining(mask)=false;
 end
 mask=remaining & va<=0 & (d4-d3)>=0 & (d5-d6)>=0;
 if any(mask)
  numerator=d4(mask)-d3(mask);
  denominator=numerator+d5(mask)-d6(mask);
  weight=numerator./denominator;
  closest(mask,:)=b(mask,:)+weight.*(c(mask,:)-b(mask,:));
  remaining(mask)=false;
 end
 if any(remaining)
  denominator=va(remaining)+vb(remaining)+vc(remaining);
  v=vb(remaining)./denominator; w=vc(remaining)./denominator;
  closest(remaining,:)=a(remaining,:)+v.*ab(remaining,:)+w.*ac(remaining,:);
 end
 distanceSquared=sum((candidate-closest).^2,2);
 intersects=any(distanceSquared<triangleLimitSquared);
 end

 function inside=pointInsideCached(candidate)
 %Reuse vertical-ray intersections while a directional path keeps XY fixed.
 if isfield(context,'occupancy') && context.occupancy.enabled
  occupancyIndex=spOccupancyCellIndex(context.occupancy,candidate);
  label=context.occupancy.labels(occupancyIndex(1),occupancyIndex(2),occupancyIndex(3));
  if label==0, inside=false; return; end
  if label==1, inside=true; return; end
 end
 pointXY=candidate(1:2);
 if ~isequal(pointXY,cachedPointXY)
  cachedPointXY=pointXY;
  rayIndex=floor((pointXY-contextLower(1:2))/context.xySize)+1;
  rayIndex=min(max(rayIndex,1),context.xyCount);
  rayKey=sprintf('%d,%d',rayIndex(1),rayIndex(2));
  if ~isKey(context.xyCells,rayKey)
   cachedRayHeights=[];
  else
   rayIds=context.xyCells(rayKey); rayIds=rayIds(:);
   if isempty(rayIds)
    cachedRayHeights=[];
   else
    ray=context.ray;
    determinant=ray.determinant(rayIds);
    valid=abs(determinant)>tolerance^2;
    alpha=zeros(size(rayIds)); beta=zeros(size(rayIds)); zHit=zeros(size(rayIds));
    if any(valid)
     validIds=rayIds(valid);
     ap=pointXY-ray.a(validIds,1:2);
     alpha(valid)=(ap(:,1).*ray.ac(validIds,2)-ap(:,2).*ray.ac(validIds,1)) ...
      .*ray.inverseDeterminant(validIds);
     beta(valid)=(ray.ab(validIds,1).*ap(:,2)-ray.ab(validIds,2).*ap(:,1)) ...
      .*ray.inverseDeterminant(validIds);
     valid=valid & alpha>=-tolerance & beta>=-tolerance & ...
      alpha+beta<=1+tolerance;
     zHit(valid)=ray.a(rayIds(valid),3)+alpha(valid).*ray.zDelta(rayIds(valid),1) ...
      +beta(valid).*ray.zDelta(rayIds(valid),2);
    end
    cachedRayHeights=sort(zHit(valid));
   end
  end
 end
 heights=cachedRayHeights(cachedRayHeights<candidate(3)-tolerance);
 inside=~isempty(heights) && mod(1+sum(diff(heights)>tolerance),2)==1;
 end
end
function v=randomPerpendicular(g)
%RANDOMPERPENDICULAR Generate a unit random vector perpendicular to gravity.
v=randn(1,3); v=v-dot(v,g)*g; if norm(v)<eps, v=[g(2),-g(1),0]; end; v=v/norm(v);
end
