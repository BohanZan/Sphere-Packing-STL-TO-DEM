function inwardNormals = spOrientInwardNormalsParallel(faceCentres, rawNormals, probeDistance, probeContext)
%SPORIENTINWARDNORMALSPARALLEL Parallel worker implementation.
% This file is called only after the wrapper has confirmed that the
% Parallel Computing Toolbox runtime is available.

inwardNormals = rawNormals;
parfor id = 1:size(rawNormals, 1)
    probe = faceCentres(id,:) + probeDistance * rawNormals(id,:);
    if ~spPointInside(probeContext, probe)
        inwardNormals(id,:) = -rawNormals(id,:);
    end
end
end
