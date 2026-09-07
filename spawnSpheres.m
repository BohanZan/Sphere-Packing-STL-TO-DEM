function [assembly, masses, totalVolume, inertia, report] = spawnSpheres(model, radii, maxAttempts, buffer, options)
%SPAWNSPHERES Static non-overlapping packing in a closed STL domain.
%   RADII is an ordered radius sequence. Initial placement visits contiguous
%   layers opposite gravity and settles every successful insertion batch.
%   MODEL may be an STL filename or a mesh struct.
%   BUFFER is b_u in the paper: h=2*max(RADII)+BUFFER and every compression
%   or shake move is capped by BUFFER. Zero selects max(RADII), as in C++.
%   Set options.maxCompressionSweeps=0 to disable movement explicitly.

%Normalise optional inputs and validate the prescribed sphere-size sequence.
if nargin < 5, options = struct; end
if nargin < 4 || isempty(buffer), buffer = 0; end
if nargin < 3 || isempty(maxAttempts), maxAttempts = 1000; end
validateattributes(radii, {'numeric'}, {'vector','real','finite','positive','nonempty'});
validateattributes(maxAttempts, {'numeric'}, {'scalar','integer','positive'});
radii = radii(:);
options = spDefaultOptions(options, maxAttempts, buffer, model);
if options.buffer == 0, options.buffer = max(radii); end

%Read the STL geometry and initialise the sparse spatial-hash state.
occupancyOptions = struct('enabled', options.occupancyAcceleration, ...
    'cellSize', options.occupancyCellSize, 'maxCells', options.occupancyMaxCells);
context = spBuildContext(model, max(radii), options.buffer, options.tolerance, occupancyOptions);
context.gravityFrame = spGravityFrame(context, options.gravity);
state = spEmptyState(context, numel(radii));

%Generate and settle initial-packing batches following Algorithm 4.
nextRadius = 1;
initialFailures = 0;
frame = context.gravityFrame;
layerBase = 0;
while nextRadius <= numel(radii) && layerBase < frame.height
    layerStart = state.count + 1;
    lastBatchStart = layerStart;
    layerRejects = 0;

    %Each call inserts until the current radius fails. Settle a successful
    %call immediately; only empty calls consume the cumulative layer budget.
    while nextRadius <= numel(radii) && layerRejects <= 2
        before = state.count;
        [state, nextRadius] = spInitialPlacement(context, state, radii, nextRadius, options, frame, layerBase);
        if state.count == before
            layerRejects = layerRejects + 1;
            continue;
        end
        lastBatchStart = before + 1;
        %Layer ownership follows insertion order even after downward motion.
        %The selected range supplies motion, energy and the cell-list update.
        switch options.initialRelaxationScope
            case 'all', activeStart = 1;
            case 'layer', activeStart = layerStart;
            case 'batch', activeStart = lastBatchStart;
        end
        state = spSettleBatch(context, state, options, activeStart, '[Initial packing]');
    end

    initialFailures = initialFailures + layerRejects;
    if state.count < layerStart
        %An empty slab may separate disconnected components of the STL.
        layerBase = layerBase + context.cellSize;
        continue;
    end

    %Only the latest generated batch supplies the top-boundary criterion.
    %Advance by h rather than the sphere top so no centre-height band is skipped.
    batchTop = -Inf;
    for id = lastBatchStart:state.count
        batchTop = max(batchTop, dot(state.centres(id,:)-frame.origin, frame.up) + state.radii(id));
    end
    if batchTop > frame.height - context.cellSize
        break;
    end
    layerBase = layerBase + context.cellSize;
end

%Try unresolved radii again by generating candidates on active surface faces.
if nextRadius <= numel(radii)
    [state, nextRadius] = spRefill(context, state, radii, nextRadius, options);
end

%Calculate volume, mass and centre of mass in the original world frame.
stlVolume = spSignedMeshVolume(context.vertices, context.faces);
validIds = 1:state.count;
worldCentres = state.centres(validIds,:);
acceptedRadii = state.radii(validIds);
%Accumulate in insertion order, matching the C++ volume and mass reductions.
masses = zeros(1, state.count);
totalVolume = 0; totalMass = 0; weighted = zeros(3,1);
for id = validIds
    r = acceptedRadii(id);
    volume = (4/3)*pi*r*r*r;
    masses(id) = options.density*volume;
    if ~isfinite(volume) || ~isfinite(masses(id)) || volume<=0 || masses(id)<=0
        error('SpherePacking:MassOverflow','Sphere mass or volume is outside the floating-point range.');
    end
    totalVolume = totalVolume + volume;
    totalMass = totalMass + masses(id);
    weighted = weighted + worldCentres(id,:).'*masses(id);
end
centreOfMass = zeros(3, 1);
if state.count > 0
    centreOfMass = weighted / totalMass;
end
if ~isfinite(totalVolume) || ~isfinite(totalMass) || any(~isfinite(centreOfMass))
    error('SpherePacking:MassOverflow','Assembly mass properties are outside the floating-point range.');
end
centredCentres = worldCentres - centreOfMass.';
centredWeighted = zeros(3,1);
for id = validIds
    centredWeighted = centredWeighted + centredCentres(id,:).'*masses(id);
end

%Choose one coordinate frame consistently for both sphere and grid outputs.
if strcmp(options.coordinateFrame, 'world')
    outputCentres = worldCentres;
    coordinateShift = zeros(1, 3);
else
    outputCentres = centredCentres;
    coordinateShift = centreOfMass.';
end

%Assemble the DEM array and calculate inertia about the physical centre of mass.
assembly = [outputCentres.'; acceptedRadii.'];
inertia = spInertia(centredCentres, acceptedRadii, masses.', centredWeighted, totalMass);
report = struct('requestedCount', numel(radii), 'acceptedCount', state.count, ...
    'unplacedCount', numel(radii) - state.count, 'stopReason', 'completed', ...
    'capacityWarning', false, 'nextUnplacedRadiusIndex', nextRadius, ...
    'initialFailures', initialFailures, 'refillPasses', options.maxRefillPasses, ...
    'initialRelaxationScope', options.initialRelaxationScope, ...
    'boundingBoxDimensions', context.upper-context.lower, 'stlVolume', stlVolume, ...
    'sphereAssemblyVolume', totalVolume, 'totalMass', totalMass, ...
    'centreOfMass', centreOfMass, 'coordinateFrame', options.coordinateFrame, ...
    'coordinateShift', coordinateShift, 'centreOfMassAfterShift', ...
    centredWeighted / max(totalMass, eps), ...
    'outputFiles', {{}});

%Record incomplete filling without discarding the valid partial assembly.
if state.count < numel(radii)
    report.stopReason = 'capacity_reached';
    report.capacityWarning = true;
    warning('SpherePacking:CapacityReached', ...
        'Requested %d spheres; placed %d. Remaining radii do not fit this geometry.', ...
        numel(radii), state.count);
end

%Write persistent results and retain the established final academic summary.
report.outputFiles = spWriteCsv(model, assembly, masses, totalVolume, inertia, report, options, context, state, coordinateShift);
spPrintSummary(report, inertia);
end

function options = spDefaultOptions(options, maxAttempts, buffer, model)
%SPDEFAULTOPTIONS Fill optional packing controls and normalise derived fields.
%Define physical, numerical and output defaults for a reproducible run.
defaults = struct('buffer', buffer, 'maxAttempts', maxAttempts, ...
    'gravity', [0 0 -1], 'tolerance', 1e-9, 'compressionTolerance', 1e-5, ...
    'maxCompressionSweeps', Inf, 'shakeSweeps', 2, 'initialRelaxationScope', 'layer', ...
    'maxRefillPasses', 3, 'randomSeed', 42, 'outputDirectory', '', 'outputPrefix', '', ...
    'density', 1.0, 'coordinateFrame', 'center_of_mass', ...
    'occupancyAcceleration', true, 'occupancyCellSize', [], 'occupancyMaxCells', 2e6);
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(options, names{k}) || isempty(options.(names{k}))
        options.(names{k}) = defaults.(names{k});
    end
end
spValidatePackingOptions(options);
options.initialRelaxationScope = char(string(options.initialRelaxationScope));
if ~ismember(options.initialRelaxationScope, {'batch','layer','all'})
    error('SpherePacking:InvalidInitialRelaxationScope', ...
        'options.initialRelaxationScope must be ''batch'', ''layer'' or ''all''.');
end

%Validate gravity before preprocessing, retaining the supplied vector.
%Each gravity frame normalises that same input once, as in the C++ version.
spGravityFrame(struct('vertices',zeros(1,3)),options.gravity);
options.coordinateFrame = char(lower(string(options.coordinateFrame)));
if ~ismember(options.coordinateFrame, {'world', 'center_of_mass'})
    error('SpherePacking:InvalidCoordinateFrame', ...
        'options.coordinateFrame must be ''world'' or ''center_of_mass''.');
end
if ~(isscalar(options.occupancyAcceleration) && ...
        (islogical(options.occupancyAcceleration) || ...
        (isnumeric(options.occupancyAcceleration) && isfinite(options.occupancyAcceleration) && ...
        any(options.occupancyAcceleration == [0 1]))))
    error('SpherePacking:InvalidOccupancyAcceleration', ...
        'options.occupancyAcceleration must be a scalar logical value.');
end
options.occupancyAcceleration = logical(options.occupancyAcceleration);
if ~isempty(options.occupancyCellSize)
    validateattributes(options.occupancyCellSize, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
end
validateattributes(options.occupancyMaxCells, {'numeric'}, ...
    {'real','finite','scalar','integer','positive','<=',double(intmax('uint32'))});

%Apply an optional random seed and infer output names from an STL filename.
if ~isempty(options.randomSeed), rng(options.randomSeed, 'twister'); end
if isempty(options.outputDirectory) && (ischar(model) || (isstring(model) && isscalar(model)))
    [options.outputDirectory, inferredPrefix] = fileparts(char(model));
    if isempty(options.outputPrefix), options.outputPrefix = inferredPrefix; end
end

%Use a generic prefix when an output directory was supplied explicitly.
if ~isempty(options.outputDirectory) && isempty(options.outputPrefix)
    options.outputPrefix = 'sphere_packing';
end
end

function volume = spSignedMeshVolume(vertices, faces)
%SPSIGNEDMESHVOLUME Calculate oriented STL volume from triangular tetrahedra.
%Use the first face vertex as the tetrahedron origin, as in the C++ version.
%This translation avoids cancellation when an STL is far from world zero.
volume = 0;
origin = vertices(faces(1,1),:);
for id = 1:size(faces, 1)
    triangle = vertices(faces(id,:), :) - origin;
    volume = volume + dot(triangle(1,:), cross(triangle(2,:), triangle(3,:)));
end
volume = volume / 6;
if ~isfinite(volume), error('SpherePacking:MeshVolumeOverflow','Mesh volume overflow.'); end
end

function spPrintSummary(report, inertia)
%SPPRINTSUMMARY Print the retained end-of-run packing report.
%Report global geometry, mass properties and final packing status.
d = report.boundingBoxDimensions;
fprintf('Bounding Box Dimensions Lx=%.8g; Ly=%.8g; Lz=%.8g\n', d(1), d(2), d(3));
fprintf('STL Volume: %.8e , Sphere Assembly Volume: %.8e\n', report.stlVolume, report.sphereAssemblyVolume);
fprintf('Sum of sphere masses: %.8e , total mass: %.8e\n', report.totalMass, report.totalMass);
disp('MI of the cluster:'); disp(inertia);
fprintf('CoM of the cluster: %.8g %.8g %.8g\n', report.centreOfMass);
fprintf('CoM of the cluster after shifting: %.8g %.8g %.8g\n', report.centreOfMassAfterShift);
fprintf('\n========================================\nFinished\n');
fprintf('Number of final spheres = %d\n', report.acceptedCount);
fprintf('Assembly volume = %.8e\n', report.sphereAssemblyVolume);
disp('Moment of inertia tensor:'); disp(inertia);
fprintf('========================================\n');
end

function inertia = spInertia(centres, radii, masses, centredWeighted, totalMass)
%SPINERTIA Calculate the inertia tensor of the sphere assembly about its CoM.
%Return a zero tensor directly for an empty packing.
inertia = zeros(3);
if isempty(radii), return; end

%Apply each sphere's intrinsic inertia and parallel-axis contribution.
com = centredWeighted / totalMass;
for k = 1:numel(radii)
    d = centres(k,:).' - com;
    inertia = inertia + (2/5)*masses(k)*radii(k)*radii(k)*eye(3) + ...
        masses(k)*(dot(d,d)*eye(3) - d*d.');
end
if any(~isfinite(inertia(:)))
    error('SpherePacking:InertiaOverflow','Assembly inertia is outside the floating-point range.');
end
end
