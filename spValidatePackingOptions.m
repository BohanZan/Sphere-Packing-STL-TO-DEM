function spValidatePackingOptions(options)
%SPVALIDATEPACKINGOPTIONS Reject invalid controls before changing the packing.
%Direct helpers may receive a partial struct; validate each supplied control.
positive={'maxAttempts','density'};
for k=1:numel(positive)
 name=positive{k};
 if isfield(options,name)
  validateattributes(options.(name),{'numeric'},{'real','finite','scalar','positive'},'',name);
 end
end
nonnegative={'buffer','tolerance','compressionTolerance'};
for k=1:numel(nonnegative)
 name=nonnegative{k};
 if isfield(options,name)
  validateattributes(options.(name),{'numeric'},{'real','finite','scalar','nonnegative'},'',name);
 end
end
integer={'maxAttempts','maxCompressionSweeps','shakeSweeps','maxRefillPasses','randomSeed'};
for k=1:numel(integer)
 name=integer{k};
 if ~isfield(options,name), continue; end
 value=options.(name);
 %An infinite compression guard means iterate until energy convergence.
 if strcmp(name,'maxCompressionSweeps') && isequal(value,Inf), continue; end
 validateattributes(value,{'numeric'},{'real','finite','scalar','integer','nonnegative'},'',name);
end
if isfield(options,'randomSeed') && options.randomSeed>double(intmax('uint32'))
 error('SpherePacking:InvalidRandomSeed','randomSeed must fit an unsigned 32-bit integer.');
end
end
