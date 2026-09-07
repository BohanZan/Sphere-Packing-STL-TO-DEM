function model=physicsBoxMesh(lower,upper)
% Small closed box shared by analytic packing-physics regressions.
if nargin==0, lower=[0 0 0]; upper=[4 4 4]; end
x=[lower(1),upper(1)]; y=[lower(2),upper(2)]; z=[lower(3),upper(3)];
model.vertices=[x(1) y(1) z(1);x(2) y(1) z(1);x(2) y(2) z(1);x(1) y(2) z(1); ...
 x(1) y(1) z(2);x(2) y(1) z(2);x(2) y(2) z(2);x(1) y(2) z(2)];
model.faces=[1 3 2;1 4 3;5 6 7;5 7 8;1 2 6;1 6 5; ...
 2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8];
end
