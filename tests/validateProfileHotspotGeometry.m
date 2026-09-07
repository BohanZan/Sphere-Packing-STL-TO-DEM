function validation=validateProfileHotspotGeometry
%VALIDATEPROFILEHOTSPOTGEOMETRY Full STL differential check, independent of packing.
project=fileparts(fileparts(mfilename('fullpath')));
addpath(project);
addpath(fullfile(project,'tests','helpers'));
[vertices,faces]=stlRead(fullfile(project,'inputs','greatBudda','greatBudda.stl'));
model=struct('vertices',vertices,'faces',faces);
options=struct('enabled',true,'cellSize',[],'maxCells',2e6);
started=tic;
context=spBuildContext(model,.625,0,1e-9,options);
f=model.faces; v=model.vertices;
raw=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
raw=raw./vecnorm(raw,2,2);
probeDistance=max(context.tolerance*100,1e-8*context.cellSize);
referenceNormals=raw;
candidateCounts=zeros(size(f,1),2);
for id=1:size(f,1)
    point=context.faceCentres(id,:)+probeDistance*raw(id,:);
    if ~referenceExactPointInside(context,point)
        referenceNormals(id,:)=-raw(id,:);
    end
    idx=min(max(floor((point(1:2)-context.lower(1:2))/context.xySize)+1,1),context.xyCount);
    key=spCellKeys(idx,context.xyKeySpec);
    if iscell(key), key=key{1}; end
    ids=spGridLookup(context.xyCells,key);
    if ~isempty(ids)
        ray=context.ray;
        keep=point(1)>=ray.xyLower(ids,1) & point(1)<=ray.xyUpper(ids,1) & ...
            point(2)>=ray.xyLower(ids,2) & point(2)<=ray.xyUpper(ids,2);
        candidateCounts(id,:)=[numel(ids),nnz(keep)];
    end
    if mod(id,20000)==0, fprintf('NORMAL_REFERENCE %d / %d\n',id,size(f,1)); end
end
assert(isequaln(context.inwardNormals,referenceNormals),'Inward normals changed');
legacy=context;
legacy.ray=rmfield(legacy.ray,{'xyLower','xyUpper'});
if isfield(legacy.ray,'xyBoundsCells')
    legacy.ray=rmfield(legacy.ray,'xyBoundsCells');
end
legacyOccupancy=spBuildOccupancyGrid(legacy,.625,options);
assert(isequal(context.occupancy,legacyOccupancy),'Occupancy changed');
rng(328,'twister');
points=context.lower+rand(2000,3).*(context.upper-context.lower);
for id=1:size(points,1)
    assert(spExactPointInside(context,points(id,:))==referenceExactPointInside(context,points(id,:)));
end
packed=load(fullfile(project,'tmp','perf_20260905','v0_ascii_plain1.mat'));
assembly=packed.result.assembly;
for id=1:size(assembly,2)
    p=assembly(1:3,id).'; r=assembly(4,id);
    assert(referenceExactPointInside(context,p));
    assert(~referenceSphereHitsTriangles(context,p,r,1:size(f,1)));
    others=id+1:size(assembly,2);
    distances=sqrt(sum((assembly(1:3,others).'-p).^2,2));
    assert(all(distances>=r+assembly(4,others).'-2*context.tolerance));
end
validation=struct('normalCount',size(f,1),'pointCount',size(points,1), ...
    'normalsExactlyEqual',true,'occupancyExactlyEqual',true, ...
    'sphereCount',size(assembly,2),'sphereGeometryValid',true, ...
    'meanRayCandidates',mean(candidateCounts,1),'seconds',toc(started));
save(fullfile(project,'tmp','perf_20260905','geometry_validation_ascii.mat'), ...
    'validation','referenceNormals','legacyOccupancy');
disp(validation);
end
