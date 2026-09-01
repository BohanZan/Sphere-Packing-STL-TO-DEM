clear;
clc;
close all;

%% ============================================================
%  Input STL
% =============================================================
fileName = '.\inputs\ironParticle\ironParticle.stl';

%% ============================================================
%  Sphere-packing parameters
% =============================================================
% One entry represents one requested DEM sphere. Replace this array later
% with radii sampled from a particle-size frequency distribution.
nTarget = 1000;
uniformRadius = 1.0e-6;
radii = repmat(uniformRadius, nTarget, 1);

maxAttempts = 100;
buffer = 0.0;
density = 1.0;

options = struct;
options.gravity = [0, 0, -1];
options.density = density;
options.maxCompressionSweeps = 100;
options.shakeSweeps = 2;
options.maxRefillPasses = 3;
options.outputDirectory = './results';
options.outputPrefix = 'Lily_Maria_packing';
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
