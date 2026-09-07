function varargout=randn(varargin)
% Test-only recorder; records MATLAB's own normal variates unchanged.
[varargout{1:nargout}]=builtin('randn',varargin{:});
global SP_CPP_NORMALS
if nargout==1 && isnumeric(varargout{1})
    SP_CPP_NORMALS{end+1}=varargout{1}(:);
end
end
