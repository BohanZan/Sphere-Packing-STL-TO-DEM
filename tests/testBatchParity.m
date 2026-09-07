function tests=testBatchParity
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
function testSameColumnStrictCutoffAndChainedHeights(testCase)
tau=2^-10; context=horizontalContext(tau*[0 .75 1.5 3],tau);
points=[repmat([.2 .2],7,1),tau*[1;1.25;1.75;2.5;2.75;4;4.25]];
verifyEqual(testCase,spExactPointInsideBatch(context,points),logical([0;1;1;1;1;1;0]));
end
function testRetainsFloatingArithmeticOrder(testCase)
tau=2^-10; context=horizontalContext(1,tau);
context.ray.zDelta=[-4,eps(1)];
verifyFalse(testCase,spExactPointInsideBatch(context,[.5 .5 tau+2^-55]));
end
function testMixedBinsOrderEmptyAndRandomState(testCase)
evalc('context=spBuildContext(physicsBoxMesh,.25,0,1e-9);');
rng(438,'twister'); points=[rand(200,3)*4;context.vertices;context.faceCentres];
expected=false(size(points,1),1);
for k=1:size(points,1), expected(k)=referenceExactPointInside(context,points(k,:)); end
state=rng;
verifyEqual(testCase,spExactPointInsideBatch(context,points),expected);
verifyEqual(testCase,spExactPointInsideBatch(context,zeros(0,3)),false(0,1));
verifyEqual(testCase,rng,state);
context.ray=rmfield(context.ray,{'xyLower','xyUpper'});
verifyEqual(testCase,spExactPointInsideBatch(context,points),expected);
end
function testDeterminantThresholdBothSigns(testCase)
tau=2^-10;
for multiplier=[.5 1 2]
 for direction=[-1 1]
  context=horizontalContext(0,tau);
  context.ray.ab=[tau 0]; context.ray.ac=[0 direction*tau*multiplier];
  context.ray.determinant=direction*tau^2*multiplier;
  context.ray.inverseDeterminant=1/context.ray.determinant;
  point=[tau/4 direction*tau*multiplier/4 1];
  verifyEqual(testCase,spExactPointInsideBatch(context,point),multiplier>1);
 end
end
end
function testScalarLegacyCharacterContext(testCase)
context=horizontalContext(0,1e-9);
context=rmfield(context,'xyKeySpec');
context.xyCells=containers.Map({'1,1'},{1});
verifyTrue(testCase,spExactPointInside(context,[.2 .2 1]));
verifyTrue(testCase,spExactPointInsideBatch(context,[.2 .2 1]));
end
function testOrientationManyProbesMatchesExactScalar(testCase)
evalc('context=spBuildContext(physicsBoxMesh,.25,0,1e-9);');
centres=repmat(context.faceCentres,88,1);
normals=repmat(context.inwardNormals,88,1);
distance=max(context.tolerance*100,1e-8*context.cellSize);
expected=normals;
for k=1:size(centres,1)
 if ~spExactPointInside(context,centres(k,:)+distance*normals(k,:)), expected(k,:)=-normals(k,:); end
end
[actual,parallel]=spOrientInwardNormals(centres,normals,distance,context,false);
verifyEqual(testCase,actual,expected);
verifyFalse(testCase,parallel);
end
function testOccupancyContextKeepsShortcutSemantics(testCase)
context=horizontalContext(0,1e-9);
context.occupancy=struct('enabled',true,'lower',[-4 -4 -4], ...
 'cellCount',[1 1 1],'cellSize',100,'labels',uint8(0));
normal=[0 0 1]; centre=[.2 .2 1];
actual=spOrientInwardNormals(centre,normal,.001,context,false);
verifyEqual(testCase,actual,-normal);
end
function testNonfiniteHeightsRetainScalarCutoff(testCase)
context=horizontalContext(0,1e-9);
points=[.2 .2 inf;.2 .2 -inf;.2 .2 NaN;NaN .2 1];
expected=false(4,1);
for k=1:4, expected(k)=spExactPointInside(context,points(k,:)); end
verifyEqual(testCase,spExactPointInsideBatch(context,points),expected);
end
function testSingleCoefficientsKeepDoubleTemporaries(testCase)
context=horizontalContext(0,1e-9);
context.ray.ab=[1 0]; context.ray.ac=[0 1];
context.ray.determinant=1; context.ray.inverseDeterminant=1;
names=fieldnames(context.ray);
for k=1:numel(names), context.ray.(names{k})=single(context.ray.(names{k})); end
point=[.5 .5+2^-24 1];
verifyFalse(testCase,spExactPointInside(context,point));
verifyEqual(testCase,spExactPointInsideBatch(context,point),spExactPointInside(context,point));
end
function testSameBinDifferentXYHaveDifferentSlopingHeights(testCase)
context=horizontalContext(1,1e-9);
context.ray.zDelta=[2 0];
points=[.25 .2 1.5;.75 .2 1.5];
verifyEqual(testCase,spExactPointInsideBatch(context,points),logical([1;0]));
end
function context=horizontalContext(heights,tolerance)
n=numel(heights); [~,spec]=spCellKeys([1 1],[1 1]);
ray=struct('a',[zeros(n,2),heights(:)],'ab',repmat([2 0],n,1), ...
 'ac',repmat([0 2],n,1),'determinant',4*ones(n,1), ...
 'inverseDeterminant',ones(n,1)/4,'zDelta',zeros(n,2));
context=struct('lower',[-4 -4 -4],'xySize',100,'xyCount',[1 1], ...
 'xyKeySpec',spec,'xyCells',spBuildStaticGrid(ones(n,1,'uint64'),(1:n).','uint64'), ...
 'ray',ray,'tolerance',tolerance,'faces',ones(n,3));
end
