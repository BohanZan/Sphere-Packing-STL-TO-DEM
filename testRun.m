%TESTRUN Execute the complete static sphere-packing regression suite.
% Run this file from MATLAB, or use: run('DEMTools/SpherePacking/testRun.m')

thisDirectory = fileparts(mfilename('fullpath'));
addpath(thisDirectory);
addpath(fullfile(thisDirectory, 'tests'));

results = runtests(fullfile(thisDirectory, 'tests'));
table(results)
assert(all([results.Passed]), 'SpherePacking:TestsFailed', ...
    'At least one SpherePacking regression test failed.');

fprintf('SpherePacking tests passed: %d tests.\n', numel(results));
