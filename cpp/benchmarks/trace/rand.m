function varargout=rand(varargin)
% Test-only recorder; do not put this directory on the normal MATLAB path.
[varargout{1:nargout}]=builtin('rand',varargin{:});
global SP_CPP_UNIFORMS
if nargout==1 && isnumeric(varargout{1})
    SP_CPP_UNIFORMS{end+1}=varargout{1}(:);
end
end
