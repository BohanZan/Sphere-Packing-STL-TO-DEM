clear;
clc;
close all;

%% ============================================================
%  Input STL
% =============================================================
fileName = '.\inputs\greatBudda\greatBudda.stl';

%% ============================================================
%  Sphere-packing parameters
% =============================================================
% One entry represents one requested DEM sphere. Replace this array later
% with radii sampled from a particle-size frequency distribution.
nTarget = 50;
uniformRadius = 0.625;
radii = repmat(uniformRadius, nTarget, 1);

maxAttempts = 60;
buffer = 0.0;
density = 1.0;

options = struct;
options.gravity = [0, 0, -1];
options.density = density;
options.maxCompressionSweeps = 200;
options.compressionTolerance = 1e-7;
options.shakeSweeps = 5;
options.maxRefillPasses = 6;
options.outputDirectory = './results/Great_Budda_profiling';
options.outputPrefix = 'Great_Budda_packing';
options.coordinateFrame = 'world'; % Choose 'world' or 'center_of_mass'.
% options.randomSeed = 42; % Uncomment to make a run reproducible.

%% ============================================================
%  Run static non-overlapping sphere packing
% =============================================================
[assembly, masses, totalVolume, inertia, report] = ...
    spawnSpheres(fileName, radii, maxAttempts, buffer, options);

%% ============================================================
%  Result files
% =============================================================
fprintf('\nCSV result files:\n');
fprintf('  %s\n', report.outputFiles{:});
