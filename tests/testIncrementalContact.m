function tests=testIncrementalContact
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.oldPath=path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(mfilename('fullpath')),'helpers'));
evalc('testCase.TestData.context=spBuildContext(physicsBoxMesh([0 0 0],[20 20 20]),.45,.1,1e-10);');
end
function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end
function testNewFacesEdgesAndCorners(testCase)
context=testCase.TestData.context;
for axes={[1 0 0],[1 1 0],[1 1 1],[-1 -1 -1]}
 v=axes{1}; u=v/norm(v);
 if v(1)>0, p=[2.5 2.5 2.5]; else, p=[6.5 6.5 6.5]; end
 q=p+1.7*v; state=spEmptyState(context,2);
 state=spAddSphere(context,state,q,.45); state=spAddSphere(context,state,p,.45);
 verifyEqual(testCase,spFirstContactDistance(context,state,2,u),1.7*norm(v)-.9,'AbsTol',1e-10);
end
end
function testNearAxisRetainsRearSideInitialContact(testCase)
context=testCase.TestData.context; context.tolerance=1e-8;
p=[10 10.1 10.5]; q=[10-1e-10 11.0000000025 10.5];
state=spEmptyState(context,2);
state=spAddSphere(context,state,q,.45); state=spAddSphere(context,state,p,.45);
u=[1 1e-8 0]; u=u/norm(u);
verifyEqual(testCase,spFirstContactDistance(context,state,2,u),0);
end
function testExactAxisSeparatesFromRearTouch(testCase)
context=testCase.TestData.context; p=[10 10.1 10.5];
state=spEmptyState(context,3);
state=spAddSphere(context,state,p-[.9 0 0],.45);
state=spAddSphere(context,state,p+[3 0 0],.45);
state=spAddSphere(context,state,p,.45);
verifyEqual(testCase,spFirstContactDistance(context,state,3,[1 0 0]),2.1,'AbsTol',1e-10);
end
