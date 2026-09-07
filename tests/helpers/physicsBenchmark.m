function result=physicsBenchmark
% Fixed small workload: compare elapsed time separately from profiler timing.
clock=tic;
context=spBuildContext(physicsBoxMesh,.25,0,1e-9);
result.preprocessingSeconds=toc(clock);
options=struct('maxAttempts',40,'maxCompressionSweeps',3, ...
 'shakeSweeps',2,'compressionTolerance',1e-7);
rng(182,'twister');
initial=spEmptyState(context,24);
[initial,~]=spInitialPlacement(context,initial,.25*ones(24,1),1,options);
result.initialCentres=initial.centres;
times=zeros(1,4);
for repetition=1:4
 state=spReindex(context,initial);
 rng(317,'twister'); clock=tic;
 state=spRelax(context,state,[0 0 -1],options,true);
 state=spRelax(context,state,[0 0 -1],options,false);
 times(repetition)=toc(clock);
end
result.relaxSeconds=times(2:end);
result.medianRelaxSeconds=median(result.relaxSeconds);
result.centres=state.centres; result.radii=state.radii;
result.count=state.count;
result.minPairGap=inf; result.minWallGap=inf;
for id=1:state.count
 others=id+1:state.count;
 gaps=sqrt(sum((state.centres(others,:)-state.centres(id,:)).^2,2)) ...
  -state.radii(others)-state.radii(id);
 result.minPairGap=min([result.minPairGap;gaps]);
 p=state.centres(id,:); r=state.radii(id);
 result.minWallGap=min([result.minWallGap,p-context.lower-r,context.upper-p-r]);
end
state=spReindex(context,initial); rng(317,'twister');
cleanup=onCleanup(@() profile('off')); %#ok<NASGU>
profile clear; profile on;
state=spRelax(context,state,[0 0 -1],options,true);
state=spRelax(context,state,[0 0 -1],options,false); %#ok<NASGU>
profile off; result.profile=profile('info');
end
