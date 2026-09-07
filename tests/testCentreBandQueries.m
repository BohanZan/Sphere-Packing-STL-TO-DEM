function tests=testCentreBandQueries
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.oldPath=path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(mfilename('fullpath')),'helpers'));
model=physicsBoxMesh; inner=physicsBoxMesh([1 1 1],[3 3 3]);
model.vertices=[model.vertices;inner.vertices]; model.faces=[model.faces;inner.faces(:,[1 3 2])+8];
evalc('testCase.TestData.context=spBuildContext(model,.25,0,1e-9);');
end
function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end
function testCentreBandMetadataAndWallBeforeCellCrossing(testCase)
context=testCase.TestData.context;
verifyTrue(testCase,context.triangleCells.centreCoverage);
state=spEmptyState(context,1); state=spAddSphere(context,state,[.7 2 2],.25);
verifyEqual(testCase,spFirstContactDistance(context,state,1,[1 0 0]),.05,'AbsTol',1e-10);
verifyFalse(testCase,spCanPlace(context,spEmptyState(context,1),[.9 2 2],.25));
end
function testPlacementMatchesAllFaceOracle(testCase)
context=testCase.TestData.context; state=spEmptyState(context,1);
rng(478,'twister');
for k=1:300
 p=4*rand(1,3); r=.05+.2*rand;
 expected=all(p-r>=context.lower & p+r<=context.upper) && ...
  referenceExactPointInside(context,p) && ...
  ~referenceSphereHitsTriangles(context,p,r,1:size(context.faces,1));
 verifyEqual(testCase,spCanPlace(context,state,p,r),expected);
end
end
function testContinuousPathMatchesAllTriangleOracle(testCase)
context=testCase.TestData.context;
rng(872,'twister');
for k=1:80
 p=.25+3.5*rand(1,3); r=.1;
 if ~referenceExactPointInside(context,p) || ...
  referenceSphereHitsTriangles(context,p,r,1:size(context.faces,1)), continue; end
 u=randn(1,3); u=u/norm(u);
 limits=inf(1,3); positive=u>0; negative=u<0;
 limits(positive)=(context.upper(positive)-r-p(positive))./u(positive);
 limits(negative)=(context.lower(negative)+r-p(negative))./u(negative);
 expected=max(0,min(limits));
 for face=1:size(context.faces,1)
  tri=context.vertices(context.faces(face,:),:);
  expected=min(expected,spTriangleContactDistance(p,r,u,tri,expected,context.tolerance));
 end
 state=spEmptyState(context,1); state=spAddSphere(context,state,p,r);
 verifyEqual(testCase,spFirstContactDistance(context,state,1,u),expected,'AbsTol',1e-10);
 for cap=[0 .05 .2 2]
  verifyEqual(testCase,spFirstContactDistance(context,state,1,u,cap),min(expected,cap),'AbsTol',1e-10);
 end
end
end
