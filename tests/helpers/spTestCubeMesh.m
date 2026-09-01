function mesh = spTestCubeMesh(sideLength)
%SPTESTCUBEMESH Construct a consistently oriented cube for shared tests.
v = [0 0 0; sideLength 0 0; sideLength sideLength 0; 0 sideLength 0; ...
     0 0 sideLength; sideLength 0 sideLength; sideLength sideLength sideLength; 0 sideLength sideLength];
f = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
     2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8];
mesh = struct('vertices', v, 'faces', f);
end
