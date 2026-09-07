function inwardNormals = spOrientInwardNormalsParallel(faceCentres, rawNormals, probeDistance, probeContext)
%SPORIENTINWARDNORMALSPARALLEL Parallel worker implementation.
% This file is called only after the wrapper has confirmed that the
% Parallel Computing Toolbox runtime is available.

inwardNormals = rawNormals;
if isfield(probeContext,'occupancy') && probeContext.occupancy.enabled
    parfor id = 1:size(rawNormals, 1)
        probe = faceCentres(id,:) + probeDistance * rawNormals(id,:);
        if ~spPointInside(probeContext, probe)
            inwardNormals(id,:) = -rawNormals(id,:);
        end
    end
    return
end
probes=faceCentres+probeDistance*rawNormals;
indices=min(max(floor((probes(:,1:2)-probeContext.lower(1:2))/probeContext.xySize)+1,1),probeContext.xyCount);
[~,~,groups]=unique(indices,'rows');
[sorted,order]=sort(groups);
if isempty(order), return; end
starts=[1;find(diff(sorted)~=0)+1;numel(order)+1];
groupCount=numel(starts)-1; blockCount=min(32,groupCount);
edges=round(linspace(1,groupCount+1,blockCount+1));
rows=cell(blockCount,1); answers=rows;
for block=1:blockCount
    rows{block}=order(starts(edges(block)):starts(edges(block+1))-1);
end
% A bin belongs wholly to one block. Workers read the same immutable context;
% no giant probe-by-triangle matrix or random-state mutation is introduced.
parfor block=1:blockCount
    answers{block}=spExactPointInsideBatch(probeContext,probes(rows{block},:));
end
for block=1:blockCount
    flip=rows{block}(~answers{block});
    inwardNormals(flip,:)=-rawNormals(flip,:);
end
end
