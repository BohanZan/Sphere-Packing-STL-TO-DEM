function state = spRelax(context, state, gravity, options, globalPass, batchStart, progressPrefix)
% Algorithms 2 and 3: continuous gravity contacts and independent local shakes.
% Only the explicit new batch moves; older spheres remain fixed obstacles.
if nargin<6, batchStart=1; end
validIds=batchStart:state.count;
if isempty(validIds), return; end
if isfield(context,'gravityFrame') && isequal(context.gravityFrame.direction,gravity)
 frame=context.gravityFrame;
else
 frame=spGravityFrame(context,gravity);
end
gravity=frame.direction;
showProgress=nargin>=7 && ~isempty(progressPrefix);
if showProgress
 if globalPass, stage='Gravity compression'; else, stage='Shaking/compression'; end
 progressClock=tic;
 completedSweeps=0;
end
% Relative heights remove dependence on the model's world-coordinate origin.
weights=state.radii(validIds).^3;
previous=sum(((state.centres(validIds,:)-frame.origin)*frame.up.').*weights);
for sweep=1:options.maxCompressionSweeps
 if showProgress
  completedSweeps=sweep;
  if sweep==1 || toc(progressClock)>=1
   fprintf('%s | %s | sweep %d/%d; batch %d spheres\n', ...
    progressPrefix,stage,sweep,options.maxCompressionSweeps,numel(validIds));
   progressClock=tic;
  end
 end
 if ~globalPass
  for shake=1:options.shakeSweeps
   state=moveSweep(context,state,gravity,validIds,true);
  end
 end
 state=moveSweep(context,state,gravity,validIds,false);
 energy=sum(((state.centres(validIds,:)-frame.origin)*frame.up.').*weights);
 % Compare against the initial energy on the first sweep too (Algorithm 2).
 % The zero case must converge, not divide by zero and exhaust the sweep cap.
 if abs(energy-previous)<=options.compressionTolerance*max(abs(previous),realmin)
  if showProgress
   fprintf('%s | %s done: converged after %d sweeps\n',progressPrefix,stage,sweep);
  end
  return;
 end
 previous=energy;
end
if showProgress
 fprintf('%s | %s done: sweep limit reached (%d/%d)\n', ...
  progressPrefix,stage,completedSweeps,options.maxCompressionSweeps);
end
end

function state=moveSweep(context,state,gravity,validIds,shake)
for id=validIds
 if shake
  direction=randomPerpendicular(gravity);
 else
  direction=gravity;
 end
 distance=spFirstContactDistance(context,state,id,direction);
 if distance>0
  state.centres(id,:)=state.centres(id,:)+distance*direction;
  % Later particles must see this particle's updated cell immediately.
  state=spReindex(context,state,id);
 end
end
end

function v=randomPerpendicular(g)
% A fresh isotropic direction in the plane normal to gravity, per particle.
v=randn(1,3); v=v-dot(v,g)*g;
if norm(v)<eps
 [~,axis]=min(abs(g)); v=zeros(1,3); v(axis)=1;
 v=v-dot(v,g)*g;
end
v=v/norm(v);
end
