function tests=testCppParity
%TESTCPPPARITY Analytic cases and replay sequences from the C++ implementation.
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

function testFirstCandidateUsesOneLayerAndAutomaticBuffer(testCase)
cleanup=useReplay([.5 .5 .5],[]); %#ok<NASGU>
options=packingOptions; options.maxCompressionSweeps=0;
evalc('[a,~,~,~,report]=spawnSpheres(box,1,1,0,options);');
verifyEqual(testCase,a,[5;5;2;1],'AbsTol',1e-12);
verifyEqual(testCase,report.initialRelaxationScope,'layer');
end

function testEachSuccessfulCallSettlesBeforeTheNextCandidate(testCase)
cleanup=useReplay(repmat(.5,1,9),[]); %#ok<NASGU>
evalc('a=spawnSpheres(box,[1;1],1,8,packingOptions);');
verifyEqual(testCase,a,[5 5;5 5;1 3;1 1],'AbsTol',1e-12);
verifyDraws(testCase,9,0);
end

function testEmptyCallsStayCumulativeAcrossSuccessfulCalls(testCase)
cleanup=useReplay(repelem([.125 .125 .125 .875 .875 .875 .875],3),[]); %#ok<NASGU>
options=packingOptions; options.maxCompressionSweeps=0;
evalc('[a,~,~,~,report]=spawnSpheres(box,[1;1;1],1,8,options);');
verifyEqual(testCase,a,[2 8;2 8;2 8;1 1],'AbsTol',1e-12);
verifyEqual(testCase,report.initialFailures,3);
verifyEqual(testCase,report.nextUnplacedRadiusIndex,3);
verifyDraws(testCase,21,0);
end

function testInitialScopesWithinOneLayer(testCase)
for scope={'batch','layer','all'}
 cleanup=useReplay([.125 .225 0 .125 .225 0 .875 .225 0],repmat([1 0 0],1,3));
 options=packingOptions; options.shakeSweeps=1; options.initialRelaxationScope=scope{1};
 evalc('a=spawnSpheres(box,[1;1],1,.25,options);');
 expected=2.5; normals=9;
 if strcmp(scope{1},'batch'), expected=2.25; normals=6; end
 verifyEqual(testCase,a(1,:),[expected 8.25],'AbsTol',1e-12);
 verifyEqual(testCase,a(3,:),[1 1],'AbsTol',1e-12);
 verifyDraws(testCase,9,normals);
 clear cleanup
end
end

function testAllScopeRevisitsCompletedLayers(testCase)
for scope={'batch','layer','all'}
 cleanup=useReplay([repmat([.125 .225 0],1,5),.875 .225 0],repmat([1 0 0],1,3));
 options=packingOptions; options.shakeSweeps=1; options.initialRelaxationScope=scope{1};
 evalc('[a,~,~,~,report]=spawnSpheres(box,[1;1],1,.25,options);');
 expected=2.25; if strcmp(scope{1},'all'), expected=2.5; end
 verifyEqual(testCase,a(1,:),[expected 8.25],'AbsTol',1e-12);
 verifyEqual(testCase,report.initialFailures,3);
 verifyEqual(testCase,report.initialRelaxationScope,scope{1});
 clear cleanup
end
end

function testReversedAndObliqueGravityUseFirstSlab(testCase)
for gravity=[0 0 1;1 1 1;-1 2 -3].'
 options=packingOptions; options.gravity=gravity.'; options.maxCompressionSweeps=0;
 options.randomSeed=47; options.maxAttempts=10000;
 evalc('[a,~,~,~,report]=spawnSpheres(box,1,10000,2,options);');
 evalc('context=spBuildContext(box,1,2,1e-9);');
 frame=spGravityFrame(context,gravity);
 height=dot(a(1:3).'-frame.origin,frame.up);
 verifyEqual(testCase,report.acceptedCount,1);
 verifyGreaterThanOrEqual(testCase,height,0);
 verifyLessThanOrEqual(testCase,height,4+1e-12);
end
end

function testSweepStartCellsCoverPreviouslyMovedSphere(testCase)
evalc('context=spBuildContext(box,1,0,1e-9);');
state=spEmptyState(context,2);
state=spAddSphere(context,state,[2 2 1],1);
state=spAddSphere(context,state,[9 2 1],1);
options=packingOptions; options.buffer=4; options.shakeSweeps=1; options.maxCompressionSweeps=1;
cleanup=useReplay([],[1 0 0 -1 0 0]); %#ok<NASGU>
[actual,status]=spRelax(context,state,[0 0 -1],options,false);
verifyEqual(testCase,actual.centres,[6 2 1;8 2 1],'AbsTol',1e-12);
verifyTrue(testCase,status.converged);
verifyEqual(testCase,actual.cellIndices,[spCellIndex(context,actual.centres(1,:));spCellIndex(context,actual.centres(2,:))]);
end

function testEnergyConvergenceIsIndependentOfModelUnits(testCase)
for scale=[1e-100 1 1e80]
 [~,spec]=spCellKeys([1 1 1],[1 1 1]);
 mesh=box; mesh.vertices=mesh.vertices*scale;
 context=struct('vertices',mesh.vertices,'faces',mesh.faces,'lower',[0 0 0], ...
  'upper',[10 10 10]*scale,'cellSize',20*scale,'cellCount',[1 1 1], ...
  'cellKeySpec',spec,'tolerance',0,'triangleCells',spBuildStaticGrid([],[],'uint64'));
 context.triangleCells.centreCoverage=true;
 state=spEmptyState(context,1); state=spAddSphere(context,state,[5 5 8]*scale,scale);
 options=packingOptions; options.buffer=scale; options.maxCompressionSweeps=20;
 [actual,status]=spRelax(context,state,[0 0 -1],options,true);
 verifyEqual(testCase,actual.centres/scale,[5 5 1],'AbsTol',1e-12);
 verifyTrue(testCase,status.converged);
end
end

function testUnconvergedPackingDoesNotWriteOutput(testCase)
cleanup=useReplay([.5 .5 .5],[]); %#ok<NASGU>
options=packingOptions; options.maxCompressionSweeps=1; options.buffer=.25;
options.outputDirectory=tempname;
verifyError(testCase,@() spawnSpheres(box,1,1,.25,options),'SpherePacking:CompressionNotConverged');
verifyFalse(testCase,isfolder(options.outputDirectory));
end

function testDefaultSweepsReachEnergyConvergence(testCase)
cleanup=useReplay([.5 .5 .5],[]); %#ok<NASGU>
options=rmfield(packingOptions,'maxCompressionSweeps');
evalc('a=spawnSpheres(box,1,1,.001,options);');
verifyEqual(testCase,a(3),1,'AbsTol',1e-12);
end

function testDegenerateFacetStillRejectsPenetration(testCase)
context=struct('tolerance',0,'triangles',struct('a',[0 0 0],'b',[0 0 0],'c',[1 0 0]));
verifyTrue(testCase,spSphereHitsTriangles(context,[.5 .01 0],.1,1));
end

function testTriangleFeaturesSurviveExtremeScales(testCase)
for scale=[1e-100 1 1e80]
 triangle=[0 0 0;2 0 0;0 2 0]*scale;
 verifyEqual(testCase,spPointTriangleDistance([.5 .5 2]*scale,triangle)/scale,2,'AbsTol',1e-12);
 distance=spTriangleContactDistance([.5 .5 2]*scale,.25*scale,[0 0 -1],triangle,4*scale,0);
 verifyEqual(testCase,distance/scale,1.75,'AbsTol',1e-12);
end
end

function testMeshTraversalIgnoresUnusedVerticesAndKeepsFaceOrder(testCase)
mesh=box; mesh.vertices(end+1,:)=[-100 -100 -100];
evalc('context=spBuildContext(mesh,1,1,1e-9);');
verifyEqual(testCase,context.lower,[0 0 0]);
mesh.faces=mesh.faces([3:12 1:2],:);
frame=spGravityFrame(mesh,[1 1 1]);
verifyEqual(testCase,frame.origin,[10 10 10]);
end

function testZeroAreaNormalsRemainFinite(testCase)
mesh=box; mesh.faces(end+1,:)=[1 1 2];
evalc('context=spBuildContext(mesh,1,1,1e-9);');
verifyEqual(testCase,context.inwardNormals(end,:),[0 0 0]);
end

function testTranslatedMeshRetainsItsVolume(testCase)
mesh=box; mesh.vertices=mesh.vertices+1e8;
cleanup=useReplay([.5 .5 .5],[]); %#ok<NASGU>
options=packingOptions; options.maxCompressionSweeps=0;
evalc('[~,~,~,~,report]=spawnSpheres(mesh,1,1,1,options);');
verifyEqual(testCase,report.stlVolume,1000,'AbsTol',1e-10);
end

function testInvalidSphereIsRejectedBeforeGridLookup(testCase)
evalc('context=spBuildContext(box,1,1,1e-9);');
state=spEmptyState(context,0);
verifyFalse(testCase,spCanPlace(context,state,[5 5 5],NaN));
verifyFalse(testCase,spCanPlace(context,state,[NaN 5 5],1));
verifyFalse(testCase,spCanPlace(context,state,[5 5 5],0));
end

function testEmptySlabsReachDisconnectedUpperComponent(testCase)
model=physicsBoxMesh([0 0 0],[2 2 2]);
upper=physicsBoxMesh([0 0 6],[2 2 8]);
model.vertices=[model.vertices;upper.vertices];
model.faces=[model.faces;upper.faces+8];
options=packingOptions; options.maxAttempts=100; options.maxCompressionSweeps=0;
evalc('assembly=spawnSpheres(model,.3*ones(100,1),100,.3,options);');
verifyTrue(testCase,any(assembly(3,:)>6));
end

function testContiguousLayersDoNotSkipNarrowChamber(testCase)
model=physicsBoxMesh([0 0 1.85],[2.2 2.2 4.05]);
other=physicsBoxMesh([5 0 2.1],[7.2 2.2 4.2]);
thin=physicsBoxMesh([10 0 0],[10.5 .5 10]);
model.vertices=[model.vertices;other.vertices;thin.vertices];
model.faces=[model.faces;other.faces+8;thin.faces+16];
options=packingOptions; options.maxAttempts=10000;
evalc('[assembly,~,~,~,report]=spawnSpheres(model,[1;1],10000,1,options);');
verifyEqual(testCase,report.acceptedCount,2);
verifyTrue(testCase,any(assembly(1,:)>=6 & assembly(1,:)<=6.2));
verifyEqual(testCase,sort(assembly(3,:)),[2.85 3.1],'AbsTol',1e-10);
end

function testContactPromotesSingleCoordinatesBeforeArithmetic(testCase)
triangle=single([0 0 0;2 0 0;0 2 0]*1e20);
p=single([.5 .5 2]*1e20); r=single(.25e20); limit=single(4e20);
expected=spTriangleContactDistance(double(p),double(r),[0 0 -1],double(triangle),double(limit),0);
actual=spTriangleContactDistance(p,r,single([0 0 -1]),triangle,limit,single(0));
verifyEqual(testCase,actual,expected);
end

function testNonfiniteQueriesAreRejectedConsistently(testCase)
evalc('context=spBuildContext(box,1,1,1e-9);');
for point=[NaN 5 5;5 Inf 5;5 5 -Inf].'
 verifyError(testCase,@() spCellIndex(context,point.'),'SpherePacking:InvalidPoint');
 verifyError(testCase,@() spPointInside(context,point.'),'SpherePacking:InvalidPoint');
 verifyError(testCase,@() spExactPointInside(context,point.'),'SpherePacking:InvalidPoint');
 verifyError(testCase,@() spExactPointInsideBatch(context,point.'),'SpherePacking:InvalidPoint');
end
end

function testOutputPrefixIsCheckedBeforeCreatingDirectory(testCase)
options=struct('outputDirectory',tempname,'outputPrefix','../packing');
verifyError(testCase,@() spWriteCsv([],[],[],[],[],[],options,[],[],[]),'SpherePacking:InvalidOutputPrefix');
verifyFalse(testCase,isfolder(options.outputDirectory));
end

function testCsvRetainsDoublePrecisionCoordinates(testCase)
evalc('context=spBuildContext(box,1,1,1e-9);');
p=[1+eps(1),2+eps(2),3-eps(3)];
state=spEmptyState(context,1); state=spAddSphere(context,state,p,.25);
options=struct('outputDirectory',tempname,'outputPrefix','precision');
report=struct('requestedCount',1,'acceptedCount',1,'unplacedCount',0, ...
 'stopReason','completed','capacityWarning',false);
files=spWriteCsv([], [p.';.25],.1,.2,eye(3),report,options,context,state,[0 0 0]);
rows=readtable(files{1});
verifyEqual(testCase,[rows.x rows.y rows.z],p);
end

function options=packingOptions
options=struct('maxAttempts',1,'maxCompressionSweeps',Inf,'shakeSweeps',0, ...
 'maxRefillPasses',0,'compressionTolerance',1e-8,'occupancyAcceleration',false, ...
 'coordinateFrame','world','randomSeed',42);
end
function model=box
model=physicsBoxMesh([0 0 0],[10 10 10]);
end
function cleanup=useReplay(uniforms,normals)
global SP_CPP_REPLAY
oldPath=path; oldRng=rng;
SP_CPP_REPLAY=struct('uniforms',uniforms,'normals',normals,'uniformCount',0,'normalCount',0);
warningState=warning('off','MATLAB:dispatcher:nameConflict');
addpath(fullfile(fileparts(mfilename('fullpath')),'helpers','cppReplay'),'-begin');
warning(warningState);
cleanup=onCleanup(@restore);
 function restore
  path(oldPath); rng(oldRng); SP_CPP_REPLAY=[];
 end
end
function verifyDraws(testCase,uniformCount,normalCount)
global SP_CPP_REPLAY
verifyEqual(testCase,SP_CPP_REPLAY.uniformCount,uniformCount);
verifyEqual(testCase,SP_CPP_REPLAY.normalCount,normalCount);
end
