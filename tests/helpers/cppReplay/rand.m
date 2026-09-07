function value=rand(varargin)
%RAND Test-only uniform replay; this folder is added only by parity tests.
global SP_CPP_REPLAY
if isempty(SP_CPP_REPLAY)
 value=builtin('rand',varargin{:}); return;
end
if isfield(SP_CPP_REPLAY,'record') && SP_CPP_REPLAY.record
 value=builtin('rand',varargin{:});
 SP_CPP_REPLAY.uniforms=[SP_CPP_REPLAY.uniforms,value(:).'];
 SP_CPP_REPLAY.uniformCount=SP_CPP_REPLAY.uniformCount+numel(value);
 return;
end
shape=zeros(varargin{:}); count=numel(shape);
ids=SP_CPP_REPLAY.uniformCount+(1:count);
if ~isempty(ids) && ids(end)>numel(SP_CPP_REPLAY.uniforms)
 error('SpherePacking:ReplayExhausted','Uniform replay exhausted.');
end
value=reshape(SP_CPP_REPLAY.uniforms(ids),size(shape));
SP_CPP_REPLAY.uniformCount=SP_CPP_REPLAY.uniformCount+count;
end
