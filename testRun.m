clear;
clc;
close all;

%% ============================================================
%  Input STL
% =============================================================
fileName = './Objects/seated-buddha-scanned-by-delgrande-and-leak/source/Lily-Maria/Lily-Maria.stl';

%% ============================================================
%  Sphere-packing parameters
% =============================================================
% One entry represents one requested DEM sphere. Replace this array later
% with radii sampled from a particle-size frequency distribution.
nTarget = 500;
uniformRadius = 0.002;
radii = repmat(uniformRadius, nTarget, 1);

maxAttempts = 2000;
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
