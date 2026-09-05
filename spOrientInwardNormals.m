function [inwardNormals, usedParallel] = spOrientInwardNormals(faceCentres, rawNormals, probeDistance, probeContext, allowParallel)
%SPORIENTINWARDNORMALS Orient face normals with an optional parallel probe pass.
% Keep PARFOR in a separate helper so installations without Parallel
% Computing Toolbox can still parse and run the serial implementation.

if nargin < 5 || isempty(allowParallel)
    allowParallel = true;
end

usedParallel = false;
if allowParallel && canUseParallelOrientation(size(rawNormals, 1))
    try
        inwardNormals = spOrientInwardNormalsParallel( ...
            faceCentres, rawNormals, probeDistance, probeContext);
        usedParallel = true;
        return;
    catch parallelError
        % A worker/pool failure must not prevent preprocessing on a MATLAB
        % installation where serial execution remains fully supported.
        warning('SpherePacking:ParallelNormalOrientationFallback', ...
            'Parallel normal orientation was unavailable (%s). Using serial execution.', ...
            parallelError.message);
    end
end

inwardNormals = orientSerial(faceCentres, rawNormals, probeDistance, probeContext);
end

function canUse = canUseParallelOrientation(faceCount)
% Avoid the pool startup and data-transfer cost on small meshes. This code
% intentionally contains no PARFOR syntax, preserving no-toolbox support.
minimumParallelFaces = 1024;
canUse = faceCount >= minimumParallelFaces;
if ~canUse
    return;
end

try
    canUse = license('test', 'Distrib_Computing_Toolbox') && ...
        exist('parpool', 'file') == 2 && exist('gcp', 'file') == 2;
catch
    canUse = false;
end
end

function inwardNormals = orientSerial(faceCentres, rawNormals, probeDistance, probeContext)
inwardNormals = rawNormals;
for id = 1:size(rawNormals, 1)
    probe = faceCentres(id,:) + probeDistance * rawNormals(id,:);
    if ~spPointInside(probeContext, probe)
        inwardNormals(id,:) = -rawNormals(id,:);
    end
end
end
