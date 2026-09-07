function tests=testSurfaceRaster
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.oldPath=path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(mfilename('fullpath')),'helpers'));
end
function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end
function testFiniteTriangleCoverageAllCells(testCase)
grid=struct('lower',[-1 -1 -1],'cellSize',.4,'cellCount',[9 8 7]);
[x,y,z]=ndgrid(1:9,1:8,1:7); allIndices=[x(:) y(:) z(:)];
centres=grid.lower+(allIndices-.5)*grid.cellSize;
triangles={[-.4 -.3 .2;1.8 -.3 .2;-.4 1.5 .2], ...
 [-.7 -.6 -.5;1.9 1.6 .8;-.4 1.1 1.3], ...
 [0 0 0;0 0 0;1.7 1.3 1.1], ...
 [0 0 0;1.7 1.3 1.1;.85 .65+.000000001 .55], ...
 [.6 .6 .6;.6 .6 .6;.6 .6 .6]};
for tri=triangles
 for radius=[sqrt(3)/2,(1+sqrt(3))/2]*grid.cellSize
  [actual,~,stats]=spTriangleCellCandidates(tri{1},grid,radius,1e-12);
  distance=oracleDistance(centres,tri{1});
  expected=allIndices(distance<=radius,:);
  verifyTrue(testCase,all(ismember(expected,actual,'rows')));
  verifyEqual(testCase,size(unique(actual,'rows'),1),size(actual,1));
  verifyTrue(testCase,all(all(actual>=1 & actual<=grid.cellCount)));
  verifyGreaterThanOrEqual(testCase,stats.enumerated,size(actual,1));
 end
end
end
function testTranslationsScalesAndPartialLastCell(testCase)
for scale=[1e-4 1 1e4]
 for shift={[0 0 0],[1e6 -2e6 3e6]*scale}
  grid=struct('lower',shift{1},'cellSize',.4*scale,'cellCount',[7 6 5]);
  tri=scale*[.2 .3 .4;2.61 1.12 1.91;.13 2.39 1.24]+shift{1};
  [x,y,z]=ndgrid(1:7,1:6,1:5); ijk=[x(:) y(:) z(:)];
  p=grid.lower+(ijk-.5)*grid.cellSize; R=(1+sqrt(3))*.2*scale;
  expected=ijk(oracleDistance(p,tri)<=R,:);
  actual=spTriangleCellCandidates(tri,grid,R,1e-12*scale);
  verifyTrue(testCase,all(ismember(expected,actual,'rows')));
 end
end
end
function testObliqueEnumerationIsNotVolumeSized(testCase)
grid=struct('lower',[0 0 0],'cellSize',1,'cellCount',[80 80 80]);
[~,~,stats]=spTriangleCellCandidates([1 1 1;78 78 3;3 76 78],grid,(1+sqrt(3))/2,1e-12);
verifyLessThan(testCase,stats.enumerated,prod(grid.cellCount)/5);
end
function testLongThinTriangleEnumeration(testCase)
grid=struct('lower',[0 0 0],'cellSize',1,'cellCount',[100 100 100]);
[~,~,stats]=spTriangleCellCandidates([1 1 1;98 98 98;49 49+1e-10 49],grid,1,1e-12);
verifyLessThan(testCase,stats.enumerated,20000);
end
function testRepeatedVertexDistance(testCase)
verifyEqual(testCase,spPointTriangleDistance([.5 1 0],[0 0 0;0 0 0;1 0 0]),1,'AbsTol',1e-14);
end
function testSkinnySlantedFaceDistanceBand(testCase)
tri=[0 0 0;1 2 3;.5+3*2^-26 1 1.5]; n=[0 3 -2]/sqrt(13);
h=5e-8; R=(1+sqrt(3))*h/2;
for side=[-1 1]
 p=mean(tri,1)+side*n*R*(1-1e-8);
 grid=struct('lower',p-.5*h,'cellSize',h,'cellCount',[1 1 1]);
 actual=spTriangleCellCandidates(tri,grid,R,0);
 verifyEqual(testCase,actual,[1 1 1]);
end
end
function testFiniteLargeScaleDoesNotSquareToInfinity(testCase)
s=1e155; tri=s*[0 0 1.3;4 0 1.3;0 4 1.3];
grid=struct('lower',[0 0 0],'cellSize',s,'cellCount',[4 4 4]);
actual=spTriangleCellCandidates(tri,grid,.9*s,0);
verifyTrue(testCase,ismember([1 1 2],actual,'rows'));
end
function testObliqueOccupancyPreservesConservativeParity(testCase)
model=physicsBoxMesh;
theta=.63; R=[cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
model.vertices=model.vertices*R.';
options=struct('enabled',true,'cellSize',.125,'maxCells',2e6);
evalc('context=spBuildContext(model,.25,0,1e-9,options);');
p=[2 2 2]*R.'; index=spOccupancyCellIndex(context.occupancy,p);
verifyTrue(testCase,spPointInside(context,p));
verifyNotEqual(testCase,context.occupancy.labels(index(1),index(2),index(3)),uint8(0));
rng(991,'twister'); points=context.lower+rand(250,3).*(context.upper-context.lower);
for k=1:size(points,1)
 verifyEqual(testCase,spPointInside(context,points(k,:)),referenceExactPointInside(context,points(k,:)));
end
end
function testOccupancyPreservesEnclosedCavity(testCase)
model=physicsBoxMesh; inner=physicsBoxMesh([1 1 1],[3 3 3]);
model.vertices=[model.vertices;inner.vertices];
model.faces=[model.faces;inner.faces(:,[1 3 2])+8];
options=struct('enabled',true,'cellSize',.125,'maxCells',2e6);
evalc('context=spBuildContext(model,.25,0,1e-9,options);');
index=spOccupancyCellIndex(context.occupancy,[2 2 2]);
verifyEqual(testCase,context.occupancy.labels(index(1),index(2),index(3)),uint8(0));
verifyFalse(testCase,spPointInside(context,[2 2 2]));
verifyTrue(testCase,spPointInside(context,[.5 .5 .5]));
end
function distance=oracleDistance(p,tri)
% All finite segments plus orthogonal face projection: independent of bins.
a=tri(1,:); b=tri(2,:); c=tri(3,:); distance=inf(size(p,1),1);
for edge=[1 2;2 3;3 1].'
 q=tri(edge(1),:); d=tri(edge(2),:)-q; length2=sum(d.^2);
 if length2==0, residual=p-q; else
  t=max(0,min(1,(p-q)*d.'/length2)); residual=p-q-t*d;
 end
 distance=min(distance,sqrt(sum(residual.^2,2)));
end
n=cross(b-a,c-a); n2=sum(n.^2); if n2==0, return; end
projection=p-((p-a)*n.'/n2)*n;
side1=cross(repmat(b-a,size(p,1),1),projection-a,2)*n.';
side2=cross(repmat(c-b,size(p,1),1),projection-b,2)*n.';
side3=cross(repmat(a-c,size(p,1),1),projection-c,2)*n.';
inside=side1>=0 & side2>=0 & side3>=0;
distance(inside)=min(distance(inside),abs((p(inside,:)-a)*n.')/sqrt(n2));
end
