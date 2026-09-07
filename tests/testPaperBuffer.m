function tests=testPaperBuffer
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

function testOneGravitySweepUsesTheSameBufferAsGrid(testCase)
for buffer=[.1 .3]
 context=makeContext(buffer);
 state=oneSphere(context,[2 2 2]);
 actual=spRelax(context,state,[0 0 -1],oneSweep(buffer),true);
 verifyEqual(testCase,context.cellSize,.5+buffer,'AbsTol',1e-14);
 verifyEqual(testCase,actual.centres,[2 2 2-buffer],'AbsTol',1e-12);
end
end

function testEachLocalShakeIsCappedAndGravityStillFollows(testCase)
context=makeContext(.1); state=spEmptyState(context,3);
points=[1 1 2;2 2 2;3 3 2];
for id=1:3, state=spAddSphere(context,state,points(id,:),.25); end
options=oneSweep; options.shakeSweeps=2;
rng(56,'twister'); expected=points;
for shake=1:2
 for id=1:3
  draw=randn(1,3); xy=draw(1:2)/norm(draw(1:2));
  expected(id,1:2)=expected(id,1:2)+.1*xy;
 end
end
expected(:,3)=1.9; expectedRng=rng;
rng(56,'twister'); actual=spRelax(context,state,[0 0 -1],options,false);
verifyEqual(testCase,actual.centres,expected,'AbsTol',1e-12);
verifyEqual(testCase,rng,expectedRng);
end

function testCloserSphereStopsBeforeBufferLimit(testCase)
context=makeContext(.3); state=spEmptyState(context,2);
state=spAddSphere(context,state,[2 2 1.3],.25);
state=spAddSphere(context,state,[2 2 2],.25);
actual=spRelax(context,state,[0 0 -1],oneSweep(.3),true,2);
verifyEqual(testCase,actual.centres,[2 2 1.3;2 2 1.8],'AbsTol',1e-12);
end

function testCloserTriangleStopsBeforeBufferLimit(testCase)
context=makeContext(.3); state=oneSphere(context,[2 2 .3]);
actual=spRelax(context,state,[0 0 -1],oneSweep(.3),true);
verifyEqual(testCase,actual.centres,[2 2 .25],'AbsTol',1e-12);
end

function testRepeatedBoundedStepsStillCrossCellsAndSettle(testCase)
context=makeContext(.1); state=oneSphere(context,[2 2 3]);
options=oneSweep; options.maxCompressionSweeps=40;
actual=spRelax(context,state,[0 0 -1],options,true);
verifyEqual(testCase,actual.centres,[2 2 .25],'AbsTol',1e-12);
verifyEqual(testCase,actual.cellIndices,spCellIndex(context,actual.centres));
verifyGreaterThan(testCase,state.cellIndices(3)-actual.cellIndices(3),3);
end

function testObliqueStepUsesLengthNotAnAxisComponent(testCase)
context=makeContext(.1); state=oneSphere(context,[2 2 2]);
u=[1 2 -2]/3;
actual=spRelax(context,state,u,oneSweep,true);
verifyEqual(testCase,actual.centres,[2 2 2]+.1*u,'AbsTol',1e-12);
end

function testZeroBufferSelectsLargestAcceptedRadius(testCase)
context=makeContext(0);
for globalPass=[true false]
 state=oneSphere(context,[2 2 2]);
 options=oneSweep(0); options.shakeSweeps=0;
 actual=spRelax(context,state,[0 0 -1],options,globalPass);
 verifyEqual(testCase,actual.centres,[2 2 1.75],'AbsTol',1e-12);
end
end

function testPublicBufferOverrideReachesBothCompressionPhases(testCase)
options=oneSweep; options.shakeSweeps=0; options.maxRefillPasses=0;
options.randomSeed=42; options.coordinateFrame='world'; options.buffer=.2;
options.occupancyAcceleration=false;
%A loose energy threshold makes each phase finish after its first sweep.
options.compressionTolerance=1;
rng(42,'twister'); draw=rand(1,3); p=.25+3.5*draw;
p(3)=.25+(.7-.25)*draw(3);
verifyGreaterThan(testCase,p(3),.45);
evalc('assembly=spawnSpheres(physicsBoxMesh,.25,1,.1,options);');
expected=p; expected(3)=max(.25,p(3)-.4);
verifyEqual(testCase,assembly(1:3).',expected,'AbsTol',1e-12);
end

function testPublicZeroBufferEnablesCompression(testCase)
options=oneSweep(0); options.randomSeed=42; options.maxRefillPasses=0;
options.maxCompressionSweeps=Inf; options.shakeSweeps=0;
options.occupancyAcceleration=false;
options.coordinateFrame='world';
evalc('assembly=spawnSpheres(physicsBoxMesh,.25,1,0,options);');
verifyEqual(testCase,assembly(3),.25,'AbsTol',1e-12);
end

function context=makeContext(buffer)
evalc('context=spBuildContext(physicsBoxMesh,.25,buffer,1e-9);');
end
function state=oneSphere(context,p)
state=spEmptyState(context,1); state=spAddSphere(context,state,p,.25);
end
function options=oneSweep(buffer)
if nargin<1, buffer=.1; end
options=struct('maxCompressionSweeps',1,'shakeSweeps',1,'compressionTolerance',1e-7,'buffer',buffer);
end
