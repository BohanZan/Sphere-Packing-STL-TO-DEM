function mesh = spTestConcavePrism()
%SPTESTCONCAVEPRISM Construct a watertight L-shaped prism.
v = [0 0 0; 3 0 0; 3 1 0; 1 1 0; 1 3 0; 0 3 0; ...
     0 0 3; 3 0 3; 3 1 3; 1 1 3; 1 3 3; 0 3 3];
bottom = [1 4 2; 2 4 3; 1 6 4; 4 6 5];
top = [7 8 10; 8 9 10; 7 10 12; 10 11 12];
sides = [1 2 8; 1 8 7; 2 3 9; 2 9 8; 3 4 10; 3 10 9; ...
         4 5 11; 4 11 10; 5 6 12; 5 12 11; 6 1 7; 6 7 12];
mesh = struct('vertices', v, 'faces', [bottom; top; sides]);
end
