function result = profileOccupancyAcceleration(runCount, config)
%PROFILEOCCUPANCYACCELERATION Compare exact-parity bypass with a fixed seed.
% This manual helper intentionally profiles outside the regular test suite.
if nargin < 1 || isempty(runCount), runCount = 3; end
if nargin < 2 || isempty(config), config = defaultConfig(); end

disabled = config;
disabled.options.occupancyAcceleration = false;
enabled = config;
enabled.options.occupancyAcceleration = true;
result = struct('disabled', profileOne(runCount, disabled), ...
    'enabled', profileOne(runCount, enabled));
result.outcomesMatch = isequal(result.disabled.acceptedCount, result.enabled.acceptedCount) && ...
    all(abs(result.disabled.sphereVolumeFraction - result.enabled.sphereVolumeFraction) <= 1e-12);
end

function result = profileOne(runCount, config)
profile clear
profile on
benchmark = benchmarkSpherePacking(runCount, config);
profile off
info = profile('info');
profile clear
names = string({info.FunctionTable.FunctionName});
matches = endsWith(names, 'spExactPointInside');
calls = sum([info.FunctionTable(matches).NumCalls]);
result = struct('exactPointInsideCalls', calls, ...
    'medianSeconds', benchmark.medianSeconds, ...
    'acceptedCount', benchmark.acceptedCount, ...
    'sphereVolumeFraction', benchmark.sphereVolumeFraction);
end

function config = defaultConfig()
root = fileparts(fileparts(mfilename('fullpath')));
config = struct('model', fullfile(root, 'inputs', 'ironParticle', 'ironParticle.stl'), ...
    'radii', repmat(3.0e-6, 100, 1), 'maxAttempts', 100, 'buffer', 0, ...
    'options', struct('randomSeed', 53, 'maxCompressionSweeps', 100, ...
    'shakeSweeps', 2, 'maxRefillPasses', 3, 'coordinateFrame', 'world'));
end
