function tests=testPackingPhysics
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

function testFineGridIgnoresLargestTriangle(testCase)
context=spBuildContext(physicsBoxMesh([0 0 0],[6 6 6]),.25,0,1e-9);
verifyEqual(testCase,context.cellSize,.5);
verifyEqual(testCase,context.cellCount,[12 12 12]);
hitX=[];
cellKeys=context.triangleCells.occupiedKeys;
if ~iscell(cellKeys), cellKeys=num2cell(cellKeys); end
for k=1:numel(cellKeys)
 if ismember(1,spGridLookup(context.triangleCells,cellKeys(k)))
  index=spCellIndices(cellKeys(k),context.cellCount); hitX(end+1)=index(1); %#ok<AGROW>
 end
end
verifyGreaterThanOrEqual(testCase,numel(unique(hitX)),12);
verifyEqual(testCase,context.xySize,context.cellSize);
end

function testFineGridCavityAndAllSphereNeighbors(testCase)
model=physicsBoxMesh([0 0 0],[6 6 6]);
inner=physicsBoxMesh([2 2 2],[4 4 4]);
model.vertices=[model.vertices;inner.vertices];
model.faces=[model.faces;inner.faces(:,[1 3 2])+8];
context=spBuildContext(model,.25,0,1e-9);
state=spEmptyState(context,1);
verifyFalse(testCase,spCanPlace(context,state,[1.9 3.2 3.2],.2));
verifyTrue(testCase,spCanPlace(context,state,[1.75 3.2 3.2],.2));
context=spBuildContext(physicsBoxMesh([0 0 0],[6 6 6]),.25,0,1e-9);
p=[3.25 3.25 3.25];
for ix=-1:1
 for iy=-1:1
  for iz=-1:1
   offset=[ix iy iz]; if ~any(offset), continue; end
   state=spEmptyState(context,1);
   state=spAddSphere(context,state,p+.26*offset,.25);
   verifyFalse(testCase,spCanPlace(context,state,p,.25));
  end
 end
end
state.centres(1,:)=[.3 .3 .3]; state=spReindex(context,state,1);
state.centres(1,:)=[5.7 5.7 5.7]; state=spReindex(context,state,[1 1]);
members=values(state.sphereCells);
verifyEqual(testCase,[members{:}],1);
verifyEqual(testCase,state.cellIndices(1,:),spCellIndex(context,state.centres(1,:)));
end

function testInvalidGridBufferRejected(testCase)
verifyError(testCase,@() spBuildContext(physicsBoxMesh,.25,-1,1e-9), ...
 'SpherePacking:InvalidGridBuffer');
end

function testReindexPreservesSharedCellMembersAndRebuildsOnce(testCase)
context=spBuildContext(physicsBoxMesh,.25,0,1e-9);
state=spEmptyState(context,3);
state=spAddSphere(context,state,[1.1 1.1 1.1],.05);
state=spAddSphere(context,state,[1.3 1.3 1.3],.05);
state=spAddSphere(context,state,[2.3 2.3 2.3],.05);
originalCell=state.cellIndices(1,:);
state.centres(1,:)=[3.7 3.7 3.7]; state=spReindex(context,state,[1 1]);
members=state.sphereCells(spCellKeys(originalCell,context.cellKeySpec));
verifyEqual(testCase,members,2);
members=values(state.sphereCells); verifyEqual(testCase,sort([members{:}]),1:3);
state=spReindex(context,state);
members=values(state.sphereCells); verifyEqual(testCase,sort([members{:}]),1:3);
for id=1:3
 verifyEqual(testCase,state.cellIndices(id,:),spCellIndex(context,state.centres(id,:)));
end
end

function testObliqueGravityPackingExportsWorldCoordinates(testCase)
R=axisRotation([1 2 3],.63); model=physicsBoxMesh;
model.vertices=model.vertices*R.'+[11 -7 3];
options=smallOptions; options.gravity=([0 0 -1]*R.').';
options.randomSeed=79; options.maxRefillPasses=1; options.coordinateFrame='world';
options.occupancyAcceleration=true;
root=fileparts(fileparts(mfilename('fullpath')));
options.outputDirectory=tempname(fullfile(root,'tmp')); options.outputPrefix='oblique_physics';
evalc('[assembly,~,~,~,report]=spawnSpheres(model,.2*ones(12,1),12,0,options);');
verifyEqual(testCase,report.acceptedCount,12);
verifyTrue(testCase,all(cellfun(@isfile,report.outputFiles)));
rows=readtable(report.outputFiles{1});
verifyEqual(testCase,[rows.x rows.y rows.z],assembly(1:3,:).','AbsTol',1e-10);
context=spBuildContext(model,.2,0,1e-9);
state=spEmptyState(context,12);
for id=1:12
 state=spAddSphere(context,state,assembly(1:3,id).',assembly(4,id));
end
verifyPacking(testCase,context,state);
end

function testGravityFrameUsesMeshSupportAndRelativeOrigin(testCase)
model=physicsBoxMesh([0 0 0],[2 4 6]);
q=[0 -1 0;1 0 0;0 0 1];
for rotation={eye(3),q,axisRotation([1 2 3],.63)}
 R=rotation{1}; rotated=model; rotated.vertices=model.vertices*R.'+[101 -73 1000];
 context=struct('vertices',rotated.vertices);
 frame=spGravityFrame(context,(-[1 0 0]*R.').');
 verifyEqual(testCase,frame.height,2,'AbsTol',1e-10);
 p=[1.3 2 3]*R.'+[101 -73 1000];
 verifyEqual(testCase,(p-frame.origin)*frame.up.',1.3,'AbsTol',1e-10);
 verifyEqual(testCase,frame.direction,-[1 0 0]*R.','AbsTol',1e-14);
end
end

function testCsvLegacyMapMatchesStaticGrid(testCase)
context=spBuildContext(physicsBoxMesh,.5,0,1e-9);
state=spEmptyState(context,1); state=spAddSphere(context,state,[2 2 2],.5);
root=fileparts(fileparts(mfilename('fullpath')));
options=struct('outputDirectory',tempname(fullfile(root,'tmp')),'outputPrefix','static');
report=struct('requestedCount',1,'acceptedCount',1,'unplacedCount',0, ...
 'stopReason','completed','capacityWarning',false);
a=spWriteCsv([], [2;2;2;.5],1,pi/6,eye(3),report,options,context,state,[0 0 0]);
legacy=context; legacy.triangleCells=containers.Map('KeyType','char','ValueType','any');
keys=context.triangleCells.occupiedKeys;
indices=spCellIndices(keys,context.cellCount);
for k=1:numel(keys)
 legacy.triangleCells(sprintf('%d,%d,%d',indices(k,:)))=spGridLookup(context.triangleCells,keys(k));
end
legacyState=state; legacyState.sphereCells=containers.Map('KeyType','char','ValueType','any');
legacyState.sphereCells(sprintf('%d,%d,%d',state.cellIndices(1,:)))=1;
options.outputPrefix='legacy';
b=spWriteCsv([], [2;2;2;.5],1,pi/6,eye(3),report,options,legacy,legacyState,[0 0 0]);
for k=1:4, verifyEqual(testCase,fileread(a{k}),fileread(b{k})); end
end

function testGravityValidation(testCase)
model=physicsBoxMesh; context=struct('vertices',model.vertices);
for g={[0 0 0],[NaN 0 1],[Inf 0 0],[1 0],[1 0 1i]}
 verifyError(testCase,@() spGravityFrame(context,g{1}),'SpherePacking:InvalidGravity');
end
a=spGravityFrame(context,[-1;0;0]); b=spGravityFrame(context,[-1e300 0 0]);
verifyEqual(testCase,a.direction,b.direction);
end

function testCompressionIsTranslationInvariant(testCase)
options=smallOptions; options.maxCompressionSweeps=8; options.compressionTolerance=.01;
output=cell(1,2);
for k=1:2
 shift=(k-1)*[101 -73 1000];
 context=spBuildContext(physicsBoxMesh(shift,shift+4),.25,0,1e-9);
 state=spEmptyState(context,3);
 for z=[3 2 1], state=spAddSphere(context,state,[2 2 z]+shift,.25); end
 state=spRelax(context,state,[0 0 -1],options,true);
 output{k}=state.centres-shift;
end
verifyEqual(testCase,output{1},output{2},'AbsTol',2e-8);
verifyEqual(testCase,output{1}(:,3),[1.25;.75;.25],'AbsTol',2e-8);
end

function testGravityCompressionAllAxesAndOblique(testCase)
model=physicsBoxMesh;
rotations={eye(3),axisRotation([0 1 0],pi/2),axisRotation([0 1 0],-pi/2), ...
 axisRotation([1 0 0],pi/2),axisRotation([1 0 0],-pi/2), ...
 axisRotation([1 0 0],pi),axisRotation([1 2 3],.71)};
options=smallOptions; options.maxCompressionSweeps=8;
for k=1:numel(rotations)
 R=rotations{k}; rotated=model; rotated.vertices=model.vertices*R.';
 context=spBuildContext(rotated,.25,0,1e-9);
 state=spEmptyState(context,3);
 for z=[3 2 1], state=spAddSphere(context,state,[2 2 z]*R.',.25); end
 state=spRelax(context,state,[0 0 -1]*R.',options,true);
 verifyEqual(testCase,state.centres*R,[2 2 1.25;2 2 .75;2 2 .25],'AbsTol',2e-8);
 verifyPacking(testCase,context,state);
end
end

function testShakingDrawsOneDirectionPerMovingParticle(testCase)
context=spBuildContext(physicsBoxMesh,.25,0,1e-9);
options=smallOptions; options.maxCompressionSweeps=1;
state=spEmptyState(context,3);
for x=[1 2 3], state=spAddSphere(context,state,[x 2 2],.25); end
rng(56,'twister'); for id=1:3, randn(1,3); end; expectedRng=rng;
rng(56,'twister'); actual=spRelax(context,state,[0 0 -1],options,false);
verifyEqual(testCase,rng,expectedRng);
verifyPacking(testCase,context,actual);
end

function testSphereGrazingIntervalCannotBeSkipped(testCase)
context=spBuildContext(physicsBoxMesh([0 0 0],[20 20 20]),1,0,1e-9);
state=spEmptyState(context,2);
state=spAddSphere(context,state,[5+1/6 5+1.999 5],1);
state=spAddSphere(context,state,[5 5 5],1);
options=smallOptions; options.maxCompressionSweeps=1;
actual=spRelax(context,state,[1 0 0],options,true,2);
expected=1/6-sqrt(4-1.999^2);
verifyEqual(testCase,actual.centres(2,1)-5,expected,'AbsTol',1e-8);
verifyEqual(testCase,actual.centres(1,:),state.centres(1,:),'AbsTol',0);
end

function testTriangleContactFeatureBranches(testCase)
tri=[0 0 0;4 0 0;0 4 0];
points=[1 1 2;1 -1 1;-1 -1 1;2 -2 0;-1 1 .5;2 -1 2;1 -1 .1];
directions=[0 0 -1;0 1 -1;1 1 -1;0 1 0;1 0 0;0 0 -1;0 sqrt(.99) .1];
expected=[1.5;sqrt(2)-.5;sqrt(3)-.5;1.5;1;inf;.5265114276384391];
for k=1:size(points,1)
 u=directions(k,:)/norm(directions(k,:));
 for winding={tri,tri([1 3 2],:)}
  distance=spTriangleContactDistance(points(k,:),.5,u,winding{1},10,1e-10);
  verifyEqual(testCase,distance,expected(k),'AbsTol',1e-9);
 end
end
end

function testTriangleStartingContactsDoNotLockTangentialMotion(testCase)
tri=[0 0 0;4 0 0;0 4 0];
verifyEqual(testCase,spTriangleContactDistance([1 1 .5],.5,[0 0 -1],tri,4,1e-9),0);
verifyEqual(testCase,spTriangleContactDistance([1 1 .5],.5,[0 0 1],tri,4,1e-9),inf);
verifyEqual(testCase,spTriangleContactDistance([1 1 .5],.5,[1 0 0],tri,4,1e-9),inf);
verifyEqual(testCase,spTriangleContactDistance([-1 0 0],.5,[1 0 0],tri,4,1e-9),.5,'AbsTol',1e-9);
end

function testLongNearlyGrazingTriangleMissIsNotFalseContact(testCase)
tri=[0 0 0;4 0 0;0 4 0];
distance=spTriangleContactDistance([-1e8 -.5001 0],.5,[1 0 0],tri,2e8,1e-9);
verifyEqual(testCase,distance,inf);
end

function testLongGrazingSphereContactIsNotMissed(testCase)
context=spBuildContext(physicsBoxMesh([-2 -2 -2],[10002 10002 2]),.5,9999,1e-9);
state=spEmptyState(context,2);
state=spAddSphere(context,state,[5999.2 8000.6 0],.5);
state=spAddSphere(context,state,[0 0 0],.5);
distance=spFirstContactDistance(context,state,2,[.6 .8 0]);
verifyEqual(testCase,distance,10000,'AbsTol',1e-5);
end

function testLongTriangleEdgeProjectionRetainsGrazingContact(testCase)
tri=[0 0 0;50000 120000 0;-12 5 0];
p=[25001.846153846152 59999.230769230766 4];
distance=spTriangleContactDistance(p,2,[0 0 -1],tri,8,1e-9);
verifyEqual(testCase,distance,4,'AbsTol',1e-5);
end

function testTriangleContactsAcrossScalesAndTranslations(testCase)
for scale=[1e-4 1 1e4]
 shift=scale*[1e3 -2e3 3e3];
 tri=scale*[0 0 0;4 0 0;0 4 0]+shift;
 distance=spTriangleContactDistance(scale*[1 -1 1]+shift,.5*scale, ...
  [0 1 -1]/sqrt(2),tri,10*scale,1e-10*scale);
 verifyEqual(testCase,distance,(sqrt(2)-.5)*scale,'AbsTol',1e-8*scale);
end
end

function testPublicBatchStopFollowsGravity(testCase)
options=smallOptions; options.shakeSweeps=0; options.maxRefillPasses=0;
options.randomSeed=19; options.coordinateFrame='world'; options.occupancyAcceleration=false;
options.gravity=[-1 0 0];
old=warning('off','SpherePacking:CapacityReached'); cleanup=onCleanup(@() warning(old)); %#ok<NASGU>
evalc('[a,~,~,~,xReport]=spawnSpheres(physicsBoxMesh([0 0 0],[10 1 1]),[.2;.6],1,0,options);');
verifyEqual(testCase,xReport.acceptedCount,1);
verifyEqual(testCase,xReport.initialFailures,6);
verifyEqual(testCase,a(1,1),.2,'AbsTol',1e-8);
options.gravity=[0;0;-1];
evalc('[a,~,~,~,zReport]=spawnSpheres(physicsBoxMesh([0 0 0],[10 1 1]),[.2;.6],1,0,options);');
verifyEqual(testCase,zReport.acceptedCount,1);
verifyEqual(testCase,zReport.initialFailures,3);
verifyEqual(testCase,a(3,1),.2,'AbsTol',1e-8);
end

function testDirectRefillNormalizesGravityAndIgnoresPerpendicularFaces(testCase)
R=axisRotation([1 2 3],.63); model=physicsBoxMesh;
model.vertices=model.vertices*R.';
context=spBuildContext(model,3,0,1e-9);
options=smallOptions; options.maxAttempts=1; options.maxRefillPasses=2;
options.gravity=1e-300*([0 0 -1]*R.').';
rng(66,'twister'); rand(1,8); expectedRng=rng;
rng(66,'twister');
evalc('[state,next]=spRefill(context,spEmptyState(context,1),3,1,options);');
verifyEqual(testCase,state.count,0); verifyEqual(testCase,next,1);
verifyEqual(testCase,rng,expectedRng);
end

function testDistantSphereAndNegativeCornerCrossing(testCase)
context=spBuildContext(physicsBoxMesh([0 0 0],[8 8 8]),.2,0,1e-9);
u=-[1 1 1]/sqrt(3);
state=spEmptyState(context,2);
state=spAddSphere(context,state,[1 1 1],.2);
state=spAddSphere(context,state,[6 6 6],.2); % Exactly on three grid planes.
distance=spFirstContactDistance(context,state,2,u);
verifyEqual(testCase,distance,5*sqrt(3)-.4,'AbsTol',1e-8);
end

function testTriangleSolverAgainstIndependentDistanceMinimization(testCase)
rng(872,'twister');
for k=1:80
 tri=randn(3,3); r=.1+.4*rand; p=3*randn(1,3);
 if spPointTriangleDistance(p,tri)<=r+1e-6, continue; end
 u=mean(tri,1)-p+.8*randn(1,3); u=u/norm(u); limit=10;
 actual=spTriangleContactDistance(p,r,u,tri,limit,1e-10);
 distance=@(t) spPointTriangleDistance(p+t*u,tri);
 [tMin,dMin]=fminbnd(distance,0,limit,optimset('TolX',1e-11));
 if dMin<r-1e-7
  lo=0; hi=tMin;
  for step=1:45
   mid=(lo+hi)/2;
   if distance(mid)>r, lo=mid; else, hi=mid; end
  end
  verifyEqual(testCase,actual,hi,'AbsTol',1e-7);
 elseif dMin>r+1e-7
  verifyEqual(testCase,actual,inf);
 end
end
end

function testSweptCellsSeeDistantCavityWall(testCase)
model=physicsBoxMesh([0 0 0],[8 4 4]);
inner=physicsBoxMesh([3 1 1],[3.02 3 3]);
model.vertices=[model.vertices;inner.vertices];
model.faces=[model.faces;inner.faces(:,[1 3 2])+8];
context=spBuildContext(model,.2,0,1e-9);
state=spEmptyState(context,1); state=spAddSphere(context,state,[.5 2 2],.2);
distance=spFirstContactDistance(context,state,1,[1 0 0]);
verifyEqual(testCase,distance,2.3,'AbsTol',1e-8);
end

function verifyPacking(testCase,context,state)
for id=1:state.count
 p=state.centres(id,:); r=state.radii(id);
 verifyTrue(testCase,referenceExactPointInside(context,p));
 verifyFalse(testCase,referenceSphereHitsTriangles(context,p,r,1:size(context.faces,1)));
 others=id+1:state.count;
 verifyTrue(testCase,all(sqrt(sum((state.centres(others,:)-p).^2,2))>= ...
  state.radii(others)+r-2*context.tolerance));
 verifyEqual(testCase,state.cellIndices(id,:),spCellIndex(context,p));
end
end
function options=smallOptions
options=struct('maxAttempts',6,'maxCompressionSweeps',2,'shakeSweeps',1, ...
 'compressionTolerance',1e-7);
end
function R=axisRotation(axis,angle)
axis=axis/norm(axis); x=axis(1);y=axis(2);z=axis(3);
K=[0 -z y;z 0 -x;-y x 0];
R=eye(3)*cos(angle)+(1-cos(angle))*(axis.'*axis)+sin(angle)*K;
end
