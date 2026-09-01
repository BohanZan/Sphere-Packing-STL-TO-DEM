function [assembly, masses, totalVolume, inertia, report] = spawnSpheres(model, radii, maxAttempts, buffer, options)
%SPAWNSPHERES Static non-overlapping packing in a closed STL domain.
%   RADII is an ordered radius sequence.  Each value is attempted exactly
%   once in the initial-placement phase; remaining values are offered to the
%   boundary-refilling phase.  MODEL may be an STL filename or a mesh struct.

if nargin < 5, options = struct; end
if nargin < 4 || isempty(buffer), buffer = 0; end
if nargin < 3 || isempty(maxAttempts), maxAttempts = 1000; end
validateattributes(radii, {'numeric'}, {'vector','real','finite','positive'});
validateattributes(maxAttempts, {'numeric'}, {'scalar','integer','positive'});
radii = radii(:);
options = spDefaultOptions(options, maxAttempts, buffer, model);
context = spBuildContext(model, max(radii), options.buffer, options.tolerance);
state = spEmptyState(context);

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

if nextRadius <= numel(radii)
    [state, nextRadius] = spRefill(context, state, radii, nextRadius, options);
end

assembly = [state.centres.'; state.radii.'];
masses = (4/3) * pi * state.radii.'.^3;
totalVolume = sum(masses);
inertia = spInertia(state.centres, state.radii, masses.');
report = struct('requestedCount', numel(radii), 'acceptedCount', state.count, ...
    'unplacedCount', numel(radii) - state.count, 'stopReason', 'completed', ...
    'capacityWarning', false, 'nextUnplacedRadiusIndex', nextRadius, ...
    'initialFailures', failedInitialBatches, 'refillPasses', options.maxRefillPasses, ...
    'outputFiles', {{}});
if state.count < numel(radii)
    report.stopReason = 'capacity_reached';
    report.capacityWarning = true;
    warning('SpherePacking:CapacityReached', ...
        'Requested %d spheres; placed %d. Remaining radii do not fit this geometry.', ...
        numel(radii), state.count);
end
report.outputFiles = spWriteCsv(model, assembly, masses, totalVolume, inertia, report, options);
end

function options = spDefaultOptions(options, maxAttempts, buffer, model)
defaults = struct('buffer', buffer, 'maxAttempts', maxAttempts, ...
    'gravity', [0 0 -1], 'tolerance', 1e-9, 'compressionTolerance', 1e-5, ...
    'maxCompressionSweeps', 100, 'shakeSweeps', 2, 'maxInitialFailures', 2, ...
    'maxRefillPasses', 3, 'randomSeed', [], 'outputDirectory', '', 'outputPrefix', '');
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(options, names{k}) || isempty(options.(names{k}))
        options.(names{k}) = defaults.(names{k});
    end
end
options.gravity = options.gravity(:).' / norm(options.gravity);
if ~isempty(options.randomSeed), rng(options.randomSeed, 'twister'); end
if isempty(options.outputDirectory) && (ischar(model) || (isstring(model) && isscalar(model)))
    [options.outputDirectory, inferredPrefix] = fileparts(char(model));
    if isempty(options.outputPrefix), options.outputPrefix = inferredPrefix; end
end
if ~isempty(options.outputDirectory) && isempty(options.outputPrefix)
    options.outputPrefix = 'sphere_packing';
end
end

function inertia = spInertia(centres, radii, masses)
inertia = zeros(3);
if isempty(radii), return; end
com = centres.' * masses / sum(masses);
for k = 1:numel(radii)
    d = centres(k,:).' - com;
    inertia = inertia + (2/5)*masses(k)*radii(k)^2*eye(3) + ...
        masses(k)*(dot(d,d)*eye(3) - d*d.');
end
end
