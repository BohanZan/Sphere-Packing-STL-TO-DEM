function result=refillProgressFixture
% Deterministic small full-refill output for before/after instrumentation QA.
model.vertices=[0 0 0;4 0 0;4 4 0;0 4 0;0 0 4;4 0 4;4 4 4;0 4 4];
model.faces=[1 3 2;1 4 3;5 6 7;5 7 8;1 2 6;1 6 5; ...
    2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8];
context=spBuildContext(model,.2,1,1e-9);
options=struct('maxAttempts',1,'maxRefillPasses',1,'gravity',[0 0 -1], ...
    'maxCompressionSweeps',Inf,'shakeSweeps',1,'compressionTolerance',1e-7,'buffer',1);
rng(1,'twister');
[state,next]=spRefill(context,spEmptyState(context,12),.2*ones(12,1),1,options);
result=struct('centres',state.centres,'radii',state.radii,'count',state.count, ...
    'next',next,'randomState',rng,'cellIndices',state.cellIndices);
end
