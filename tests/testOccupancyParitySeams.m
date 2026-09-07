function tests=testOccupancyParitySeams
tests=functiontests(localfunctions);
end
function testInteriorProjectionSeamMatchesExactParity(testCase)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root,fullfile(root,'tests','helpers'));
model=physicsBoxMesh; theta=.63;
R=[cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
model.vertices=model.vertices*R.';
options=struct('enabled',true,'cellSize',.125,'maxCells',2e6);
evalc('context=spBuildContext(model,.25,0,1e-9,options);');
tau=context.tolerance;
points=[4*cos(theta)-2*tau*sin(theta),2,-.4;2,2,-.4];
expected=false(2,1); actual=expected;
for k=1:2
 expected(k)=referenceExactPointInside(context,points(k,:));
 actual(k)=spPointInside(context,points(k,:));
end
verifyEqual(testCase,expected,[false;true]);
verifyEqual(testCase,actual,expected);
end
