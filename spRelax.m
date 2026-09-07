function [state,status] = spRelax(context, state, gravity, options, globalPass, batchStart, progressPrefix)
%SPRELAX Algorithms 2 and 3: gravity contacts and independent lateral shakes.
%Only the selected range moves; all other spheres remain fixed obstacles.
% Each move obeys Eq.22: min(sphere contact, triangle contact, buffer=b_u).
if nargin<6, batchStart=1; end
status=struct('sweeps',0,'converged',true);
validIds=batchStart:state.count;
if isempty(validIds), return; end
spValidatePackingOptions(options);
maxRadius=max(state.radii(1:state.count));
buffer=0;
if isfield(options,'buffer'), buffer=options.buffer; end
if buffer==0, buffer=maxRadius; end
%Coordinates change within a sweep while cell membership stays at its start.
%Include the maximum earlier displacement when querying those recorded cells.
halo=min(ceil((2*maxRadius+buffer)/context.cellSize),max(context.cellCount));
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
end
%Normalise common height/radius scales before the relative-energy comparison.
%Accumulating in insertion order avoids underflow and matches the C++ loop.
previous=energy(context,state,validIds,frame,maxRadius);
while status.sweeps<options.maxCompressionSweeps
 status.sweeps=status.sweeps+1;
 sweep=status.sweeps;
 if showProgress
  if sweep==1 || toc(progressClock)>=1
   fprintf('%s | %s | sweep %d/%d; batch %d spheres\n', ...
    progressPrefix,stage,sweep,options.maxCompressionSweeps,numel(validIds));
   progressClock=tic;
  end
 end
 if ~globalPass
  for shake=1:options.shakeSweeps
   state=moveSweep(context,state,gravity,validIds,true,buffer,halo);
  end
 end
 state=moveSweep(context,state,gravity,validIds,false,buffer,halo);
 current=energy(context,state,validIds,frame,maxRadius);
 % Compare against the initial energy on the first sweep too (Algorithm 2).
 % The zero case must converge, not divide by zero and exhaust the sweep cap.
 if abs(current-previous)<=options.compressionTolerance*max(abs(previous),realmin)
  if showProgress
   fprintf('%s | %s done: converged after %d sweeps\n',progressPrefix,stage,sweep);
  end
  return;
 end
 previous=current;
end
status.converged=options.maxCompressionSweeps==0;
if showProgress
 if status.converged
  fprintf('%s | %s done: disabled (0 sweeps)\n',progressPrefix,stage);
 else
  fprintf('%s | %s done: sweep limit reached (%d/%d)\n', ...
   progressPrefix,stage,status.sweeps,options.maxCompressionSweeps);
 end
end
end

function value=energy(~,state,validIds,frame,maxRadius)
%ENERGY Calculate dimensionless potential energy for exactly the moving range.
value=0;
for id=validIds
 relativeRadius=state.radii(id)/maxRadius;
 height=dot(state.centres(id,:)-frame.origin,frame.up)/frame.height;
 value=value+height*relativeRadius*relativeRadius*relativeRadius;
end
if ~isfinite(value), error('SpherePacking:EnergyOverflow','Nonfinite relaxation energy.'); end
end

function state=moveSweep(context,state,gravity,validIds,shake,buffer,halo)
%MOVESWEEP Update coordinates in order, then update cell lists in the same order.
for id=validIds
 if shake
  direction=randomPerpendicular(gravity);
 else
  direction=gravity;
 end
 distance=spFirstContactDistance(context,state,id,direction,buffer,halo);
 if distance>0
  state.centres(id,:)=state.centres(id,:)+distance*direction;
 end
end
%Algorithms 2/3 finish the entire gravity or shake sweep before reindexing.
state=spReindex(context,state,validIds);
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
