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
nTarget = 44000; % Match the C++ --preset buddha request.
uniformRadius = 0.3;
radii = repmat(uniformRadius, nTarget, 1);

maxAttempts = 60;
buffer = 0; % Paper b_u: zero selects max(radii) for grid padding and each move.
density = 1.0;

options = struct;
options.gravity = [0, 0, -1];
options.density = density;
options.maxCompressionSweeps = Inf; % Iterate until energy convergence; 0 disables movement.
options.compressionTolerance = 1e-8;
options.shakeSweeps = 4;
options.maxRefillPasses = 5; % Cumulative empty surface traversals, not successful rounds.
options.initialRelaxationScope = 'layer'; % Move all spheres generated in the current layer.
options.outputDirectory = './results/Great_Budda';
options.outputPrefix = 'Great_Budda_packing';
options.coordinateFrame = 'world'; % Choose 'world' or 'center_of_mass'.
options.randomSeed = 42; % Repeatable within MATLAB; C++ uses different native random draws.

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
