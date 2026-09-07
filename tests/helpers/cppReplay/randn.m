function value=randn(varargin)
%RANDN Test-only normal replay; production uses MATLAB's native generator.
global SP_CPP_REPLAY
if isempty(SP_CPP_REPLAY)
 value=builtin('randn',varargin{:}); return;
end
if isfield(SP_CPP_REPLAY,'record') && SP_CPP_REPLAY.record
 value=builtin('randn',varargin{:});
 SP_CPP_REPLAY.normals=[SP_CPP_REPLAY.normals,value(:).'];
 SP_CPP_REPLAY.normalCount=SP_CPP_REPLAY.normalCount+numel(value);
 return;
end
shape=zeros(varargin{:}); count=numel(shape);
ids=SP_CPP_REPLAY.normalCount+(1:count);
if ~isempty(ids) && ids(end)>numel(SP_CPP_REPLAY.normals)
 error('SpherePacking:ReplayExhausted','Normal replay exhausted.');
end
value=reshape(SP_CPP_REPLAY.normals(ids),size(shape));
SP_CPP_REPLAY.normalCount=SP_CPP_REPLAY.normalCount+count;
end
