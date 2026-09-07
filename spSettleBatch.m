function state=spSettleBatch(context,state,options,activeStart,progressPrefix)
%SPSETTLEBATCH Apply gravity compression, then shaking/compression to one range.
%A finite guard is diagnostic: unconverged states must not reach CSV output.
if nargin<5, progressPrefix=''; end
for globalPass=[true false]
 [state,status]=spRelax(context,state,options.gravity,options,globalPass,activeStart,progressPrefix);
 if ~status.converged
  if globalPass, stage='Gravity compression'; else, stage='Shaking/compression'; end
  error('SpherePacking:CompressionNotConverged', ...
   '%s did not converge before maxCompressionSweeps; no packing output written. Use Inf to run until energy convergence.',stage);
 end
end
end
