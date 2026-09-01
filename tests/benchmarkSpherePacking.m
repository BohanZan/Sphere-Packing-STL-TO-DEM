function result = benchmarkSpherePacking(runCount, config)
%BENCHMARKSPHEREPACKING Run reproducible wall-clock packing measurements.
%   RESULT = BENCHMARKSPHEREPACKING(RUNCOUNT) measures the iron-particle
%   case RUNCOUNT times.  An optional CONFIG is accepted for small smoke
%   tests and must contain model, radii, maxAttempts, buffer, and options.

if nargin < 1 || isempty(runCount), runCount = 3; end
validateattributes(runCount, {'numeric'}, {'scalar', 'integer', 'positive'});
if nargin < 2 || isempty(config), config = defaultConfig(); end

required = {'model', 'radii', 'maxAttempts', 'buffer', 'options'};
for id = 1:numel(required)
    if ~isfield(config, required{id})
        error('SpherePacking:InvalidBenchmarkConfig', ...
            'config.%s is required.', required{id});
    end
end

result = struct('seconds', zeros(1, runCount), 'medianSeconds', NaN, ...
    'acceptedCount', zeros(1, runCount), 'sphereVolumeFraction', zeros(1, runCount), ...
    'report', struct());
for runId = 1:runCount
    outputDirectory = tempname;
    cleanup = onCleanup(@() removeOutputDirectory(outputDirectory));
    options = config.options;
    options.outputDirectory = outputDirectory;
    options.outputPrefix = sprintf('benchmark_run_%d', runId);

    started = tic;
    [~, ~, ~, ~, report] = spawnSpheres(config.model, config.radii, ...
        config.maxAttempts, config.buffer, options);
    result.seconds(runId) = toc(started);
    result.acceptedCount(runId) = report.acceptedCount;
    result.sphereVolumeFraction(runId) = report.sphereAssemblyVolume / report.stlVolume;
    result.report = report;
    clear cleanup
end
result.medianSeconds = median(result.seconds);
end

function config = defaultConfig()
%DEFAULTCONFIG Define the representative case used for before/after timing.
root = fileparts(fileparts(mfilename('fullpath')));
config = struct('model', fullfile(root, 'inputs', 'ironParticle', 'ironParticle.stl'), ...
    'radii', repmat(3.0e-6, 100, 1), 'maxAttempts', 100, 'buffer', 0, ...
    'options', struct('randomSeed', 53, 'maxCompressionSweeps', 100, ...
    'shakeSweeps', 2, 'maxRefillPasses', 3, 'coordinateFrame', 'world'));
end

function removeOutputDirectory(pathName)
%REMOVEOUTPUTDIRECTORY Delete only the explicitly created benchmark directory.
if isfolder(pathName)
    rmdir(pathName, 's');
end
end
