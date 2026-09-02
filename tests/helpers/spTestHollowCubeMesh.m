function mesh = spTestHollowCubeMesh()
%SPTESTHOLLOWCUBEMESH Construct an outer cube with a reversed inner cavity.
outer = spTestCubeMesh(20);
innerVertices = [5 5 5; 15 5 5; 15 15 5; 5 15 5; ...
                 5 5 15; 15 5 15; 15 15 15; 5 15 15];
innerFaces = outer.faces(:, [1 3 2]) + 8;
mesh = struct('vertices', [outer.vertices; innerVertices], ...
    'faces', [outer.faces; innerFaces]);
end
