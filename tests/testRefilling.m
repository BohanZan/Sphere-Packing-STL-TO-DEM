function tests=testRefilling
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

function testEveryFaceIsRetriedUntilRejectionBudget(testCase)
% A first-empty return or random-face sampling misses these exact 24 trials.
context=spBuildContext(boxMesh,3,0,1e-9);
options=smallOptions;
options.maxAttempts=3; options.maxRefillPasses=4;
verifyEqual(testCase,nnz(context.inwardNormals*options.gravity.'>0),2);
state=spEmptyState(context,1);
cleanup=onCleanup(@() profile('off')); %#ok<NASGU>
profile clear; profile on;
[state,next]=spRefill(context,state,3,1,options);
profile off; info=profile('info');
entry=strcmp({info.FunctionTable.FunctionName},'spCanPlace');
verifyEqual(testCase,sum([info.FunctionTable(entry).NumCalls]),24);
verifyEqual(testCase,state.count,0);
verifyEqual(testCase,next,1);
end

function testGravityMovesOnlyNewBatchAndKeepsOldSphereAsObstacle(testCase)
context=spBuildContext(boxMesh,.25,1,1e-9);
options=smallOptions;
state=spEmptyState(context,2);
state=spAddSphere(context,state,[2 2 1],.25);
state=spAddSphere(context,state,[2 2 3],.25);
state=spRelax(context,state,options.gravity,options,true,2);
verifyEqual(testCase,state.centres(1,:),[2 2 1],'AbsTol',0);
verifyEqual(testCase,state.centres(2,:),[2 2 1.5],'AbsTol',2e-8);
verifyGeometry(testCase,context,state);
end

function testShakingMovesOnlyNewBatch(testCase)
context=spBuildContext(boxMesh,.25,1,1e-9);
options=smallOptions;
state=spEmptyState(context,2);
state=spAddSphere(context,state,[1 1 2],.25);
state=spAddSphere(context,state,[3 3 3],.25);
rng(64,'twister');
state=spRelax(context,state,options.gravity,options,false,2);
verifyEqual(testCase,state.centres(1,:),[1 1 2],'AbsTol',0);
verifyLessThan(testCase,state.centres(2,3),3);
verifyGeometry(testCase,context,state);
end

function testEmptyBatchDoesNotMoveOldSpheresOrConsumeRandomness(testCase)
context=spBuildContext(boxMesh,.25,1,1e-9);
options=smallOptions;
state=spEmptyState(context,1);
state=spAddSphere(context,state,[2 2 3],.25);
before=rng;
state=spRelax(context,state,options.gravity,options,false,2);
verifyEqual(testCase,state.centres,[2 2 3],'AbsTol',0);
verifyEqual(testCase,rng,before);
end

function testDisabledAndExhaustedRefillingAreNoOps(testCase)
context=spBuildContext(boxMesh,.25,0,1e-9);
options=smallOptions;
state=spEmptyState(context,1);
options.maxRefillPasses=0;
before=rng;
[actual,next]=spRefill(context,state,.25,1,options);
verifyEqual(testCase,actual.count,0); verifyEqual(testCase,next,1);
options.maxRefillPasses=3;
[actual,next]=spRefill(context,state,.25,2,options);
verifyEqual(testCase,actual.count,0); verifyEqual(testCase,next,2);
verifyEqual(testCase,rng,before);
end

function testSuccessfulCyclesDoNotConsumeRejectionBudget(testCase)
context=spBuildContext(boxMesh,.2,1,1e-9);
options=smallOptions;
options.maxAttempts=1; options.maxRefillPasses=1;
rng(1,'twister');
state=spEmptyState(context,12);
cleanup=onCleanup(@() profile('off')); %#ok<NASGU>
profile clear; profile on;
[state,next]=spRefill(context,state,.2*ones(12,1),1,options);
profile off; info=profile('info');
entry=strcmp({info.FunctionTable.FunctionName},'spRelax');
% With Mr=1, multiple nonempty traversals must still be allowed.
verifyGreaterThan(testCase,sum([info.FunctionTable(entry).NumCalls]),2);
verifyEqual(testCase,state.count,12); verifyEqual(testCase,next,13);
verifyGeometry(testCase,context,state);
end

function testEmptySweepCanRecoverAndRejectionsAreNotResetBySuccess(testCase)
context=spBuildContext(boxMesh,3,1,1e-9);
options=smallOptions;
options.maxAttempts=1; options.maxRefillPasses=1;
radii=[.2*ones(12,1);3]; % Last radius guarantees eventual rejection.
% Seed 2 samples near top-face edges: both first-sweep candidates fail.
rng(2,'twister');
[empty,next]=spRefill(context,spEmptyState(context,13),radii,1,options);
verifyEqual(testCase,empty.count,0); verifyEqual(testCase,next,1);
% Resume with one remaining rejection. This must equal Mr=2 from the start,
% even though successful insertion/compression cycles intervene.
[expected,expectedNext]=spRefill(context,empty,radii,next,options);
expectedRng=rng;
verifyGreaterThan(testCase,expected.count,0);
options.maxRefillPasses=2;
rng(2,'twister');
[actual,actualNext]=spRefill(context,spEmptyState(context,13),radii,1,options);
verifyEqual(testCase,actual.centres,expected.centres,'AbsTol',0);
verifyEqual(testCase,actual.radii,expected.radii,'AbsTol',0);
verifyEqual(testCase,actualNext,expectedNext);
verifyEqual(testCase,rng,expectedRng);
verifyGeometry(testCase,context,actual);
end

function testRefillPreservesOldSpheresInConvexConcaveAndCavityDomains(testCase)
for kind=1:3
    context=spBuildContext(refillMesh(kind),.24,1,1e-9);
    options=smallOptions;
    options.maxAttempts=6;
    radii=linspace(.18,.24,9).';
    state=spEmptyState(context,9);
    state=spAddSphere(context,state,[.5 .5 .5],radii(1));
    rng(71,'twister');
    [state,next]=spRefill(context,state,radii,2,options);
    verifyGreaterThan(testCase,state.count,1);
    verifyEqual(testCase,state.centres(1,:),[.5 .5 .5],'AbsTol',0);
    verifyEqual(testCase,state.radii(1:state.count),radii(1:state.count));
    verifyEqual(testCase,next,state.count+1);
    verifyGeometry(testCase,context,state);
end
end

function testSpawnSpheresExportsPartialAssemblyAfterRefilling(testCase)
options=smallOptions;
options.maxAttempts=1; options.randomSeed=7;
options.coordinateFrame='world';
root=fileparts(fileparts(mfilename('fullpath')));
options.outputDirectory=tempname(fullfile(root,'tmp'));
options.outputPrefix='refill_smoke';
previousWarning=warning('off','SpherePacking:CapacityReached');
cleanup=onCleanup(@() warning(previousWarning)); %#ok<NASGU>
[assembly,masses,volume,~,report]=spawnSpheres(boxMesh,[.2*ones(8,1);3],1,1,options);
verifyEqual(testCase,report.acceptedCount,8);
verifyEqual(testCase,report.nextUnplacedRadiusIndex,9);
verifyEqual(testCase,report.stopReason,'capacity_reached');
verifySize(testCase,assembly,[4 8]);
verifyEqual(testCase,numel(report.outputFiles),4);
verifyTrue(testCase,all(cellfun(@isfile,report.outputFiles)));
rows=readtable(report.outputFiles{1});
verifyEqual(testCase,height(rows),8);
verifyEqual(testCase,[rows.x rows.y rows.z],assembly(1:3,:).','AbsTol',1e-12);
verifyEqual(testCase,rows.radius,assembly(4,:).','AbsTol',1e-12);
verifyEqual(testCase,sum(masses),volume,'AbsTol',1e-12);
end

function testProgressReportsEmptyRoundsAndRejectionExit(testCase)
context=spBuildContext(boxMesh,3,0,1e-9);
options=smallOptions; options.maxRefillPasses=2; options.maxAttempts=1;
state=spEmptyState(context,1);
output=evalc('[state,next]=spRefill(context,state,3,1,options);');
verifyTrue(testCase,contains(output,'Round 1'));
verifyTrue(testCase,contains(output,'Round 2'));
verifyTrue(testCase,contains(output,'faces 2/2'));
verifyTrue(testCase,contains(output,'added 0; total 0/1'));
verifyTrue(testCase,contains(output,'rejections 1/2'));
verifyTrue(testCase,contains(output,'rejections 2/2'));
verifyTrue(testCase,contains(output,'rejection limit reached'));
verifyFalse(testCase,contains(output,'Gravity compression'));
verifyFalse(testCase,contains(output,'Shaking/compression'));
verifyEqual(testCase,next,1);
end

function testProgressReportsStagesBeforeFinalCompletion(testCase)
output=evalc('result=refillProgressFixture;');
gravity=strfind(output,'Gravity compression');
shaking=strfind(output,'Shaking/compression');
finished=strfind(output,'prescribed radii exhausted');
verifyNotEmpty(testCase,gravity); verifyNotEmpty(testCase,shaking);
verifyNotEmpty(testCase,finished);
if ~isempty(gravity) && ~isempty(shaking) && ~isempty(finished)
    verifyLessThan(testCase,gravity(1),shaking(1));
    verifyLessThan(testCase,shaking(end),finished(end));
end
verifyTrue(testCase,contains(output,'sweep 1/Inf'));
verifyTrue(testCase,contains(output,'total 12/12'));
verifyEqual(testCase,result.count,12);
end

function testProgressReportsSkippedRefillingReasons(testCase)
context=spBuildContext(boxMesh,.2,0,1e-9);
options=smallOptions;
state=spEmptyState(context,1);
options.maxRefillPasses=0;
output=evalc('spRefill(context,state,.2,1,options);');
verifyTrue(testCase,contains(output,'disabled'));
options.maxRefillPasses=1;
output=evalc('spRefill(context,state,.2,2,options);');
verifyTrue(testCase,contains(output,'prescribed radii exhausted'));
end

function verifyGeometry(testCase,context,state)
for id=1:state.count
    p=state.centres(id,:); r=state.radii(id);
    verifyTrue(testCase,referenceExactPointInside(context,p));
    verifyFalse(testCase,referenceSphereHitsTriangles(context,p,r,1:size(context.faces,1)));
    others=id+1:state.count;
    distances=sqrt(sum((state.centres(others,:)-p).^2,2));
    verifyTrue(testCase,all(distances>=r+state.radii(others)-2*context.tolerance));
    verifyEqual(testCase,state.cellIndices(id,:),spCellIndex(context,p));
    neighbours=spSphereNeighbours(context,state.sphereCells,state.cellIndices(id,:));
    verifyTrue(testCase,ismember(id,neighbours));
end
end

function options=smallOptions
options=struct('maxAttempts',4,'maxRefillPasses',3,'gravity',[0 0 -1], ...
    'maxCompressionSweeps',Inf,'shakeSweeps',1,'compressionTolerance',1e-7,'buffer',1);
end

function model=refillMesh(kind)
model=boxMesh;
if kind==2
    polygon=[0 0;4 0;4 1;1 1;1 4;0 4];
    model.vertices=[polygon zeros(6,1);polygon 4*ones(6,1)];
    cap=[1 2 3;1 3 4;1 4 6;4 5 6];
    model.faces=[cap(:,[1 3 2]);cap+6];
    for id=1:6
        next=mod(id,6)+1;
        model.faces=[model.faces;id next next+6;id next+6 id+6]; %#ok<AGROW>
    end
elseif kind==3
    model.vertices=[model.vertices;1+model.vertices/2];
    model.faces=[model.faces;model.faces(:,[1 3 2])+8];
end
end

function model=boxMesh
model.vertices=[0 0 0;4 0 0;4 4 0;0 4 0;0 0 4;4 0 4;4 4 4;0 4 4];
model.faces=[1 3 2;1 4 3;5 6 7;5 7 8;1 2 6;1 6 5; ...
    2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8];
end
