function [assembly, masses, totalVolume, inertia, report] = spawnSpheres(model, radii, maxAttempts, buffer, options)
%SPAWNSPHERES Static non-overlapping packing in a closed STL domain.
%   RADII is an ordered radius sequence.  Each value is attempted exactly
%   once in the initial-placement phase; remaining values are offered to the
%   boundary-refilling phase.  MODEL may be an STL filename or a mesh struct.

%Normalise optional inputs and validate the prescribed sphere-size sequence.
if nargin < 5, options = struct; end
if nargin < 4 || isempty(buffer), buffer = 0; end
if nargin < 3 || isempty(maxAttempts), maxAttempts = 1000; end
validateattributes(radii, {'numeric'}, {'vector','real','finite','positive'});
validateattributes(maxAttempts, {'numeric'}, {'scalar','integer','positive'});
radii = radii(:);
options = spDefaultOptions(options, maxAttempts, buffer, model);

%Read the STL geometry and initialise the sparse spatial-hash state.
occupancyOptions = struct('enabled', options.occupancyAcceleration, ...
    'cellSize', options.occupancyCellSize, 'maxCells', options.occupancyMaxCells);
context = spBuildContext(model, max(radii), options.buffer, options.tolerance, occupancyOptions);
state = spEmptyState(context, numel(radii));

%Populate the interior from random trial centres, followed by relaxation.
nextRadius = 1;
failedInitialBatches = 0;
while nextRadius <= numel(radii) && failedInitialBatches <= options.maxInitialFailures
    before = state.count;
    [state, nextRadius] = spInitialPlacement(context, state, radii, nextRadius, options);
    if state.count == before
        failedInitialBatches = failedInitialBatches + 1;
        break;
    end
    state = spRelax(context, state, options.gravity, options, true);
    state = spRelax(context, state, options.gravity, options, false);
end

%Try unresolved radii again by generating candidates on active surface faces.
if nextRadius <= numel(radii)
    [state, nextRadius] = spRefill(context, state, radii, nextRadius, options);
end

%Calculate volume, mass and centre of mass in the original world frame.
validIds = 1:state.count;
worldCentres = state.centres(validIds,:);
acceptedRadii = state.radii(validIds);
unitVolumes = (4/3) * pi * acceptedRadii.^3;
masses = (options.density * unitVolumes).';
totalVolume = sum(unitVolumes);
centreOfMass = zeros(3, 1);
if state.count > 0
    centreOfMass = worldCentres.' * masses.' / sum(masses);
end
centredCentres = worldCentres - centreOfMass.';

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
inertia = spInertia(centredCentres, acceptedRadii, masses.');
stlVolume = spSignedMeshVolume(context.vertices, context.faces);
report = struct('requestedCount', numel(radii), 'acceptedCount', state.count, ...
    'unplacedCount', numel(radii) - state.count, 'stopReason', 'completed', ...
    'capacityWarning', false, 'nextUnplacedRadiusIndex', nextRadius, ...
    'initialFailures', failedInitialBatches, 'refillPasses', options.maxRefillPasses, ...
    'boundingBoxDimensions', context.upper-context.lower, 'stlVolume', stlVolume, ...
    'sphereAssemblyVolume', totalVolume, 'totalMass', sum(masses), ...
    'centreOfMass', centreOfMass, 'coordinateFrame', options.coordinateFrame, ...
    'coordinateShift', coordinateShift, 'centreOfMassAfterShift', ...
    centredCentres.' * masses.' / max(sum(masses), eps), ...
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
    'maxCompressionSweeps', 100, 'shakeSweeps', 2, 'maxInitialFailures', 2, ...
    'maxRefillPasses', 3, 'randomSeed', [], 'outputDirectory', '', 'outputPrefix', '', ...
    'density', 1.0, 'coordinateFrame', 'center_of_mass', ...
    'occupancyAcceleration', true, 'occupancyCellSize', [], 'occupancyMaxCells', 2e6);
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(options, names{k}) || isempty(options.(names{k}))
        options.(names{k}) = defaults.(names{k});
    end
end

%Normalise the gravity direction and validate the requested output frame.
options.gravity = options.gravity(:).' / norm(options.gravity);
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
        {'real','finite','scalar','positive'});
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
%Sum the signed volumes formed by every triangle and the coordinate origin.
volume = 0;
for id = 1:size(faces, 1)
    triangle = vertices(faces(id,:), :);
    volume = volume + dot(triangle(1,:), cross(triangle(2,:), triangle(3,:))) / 6;
end
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

function inertia = spInertia(centres, radii, masses)
%SPINERTIA Calculate the inertia tensor of the sphere assembly about its CoM.
%Return a zero tensor directly for an empty packing.
inertia = zeros(3);
if isempty(radii), return; end

%Apply each sphere's intrinsic inertia and parallel-axis contribution.
com = centres.' * masses / sum(masses);
for k = 1:numel(radii)
    d = centres(k,:).' - com;
    inertia = inertia + (2/5)*masses(k)*radii(k)^2*eye(3) + ...
        masses(k)*(dot(d,d)*eye(3) - d*d.');
end
end
