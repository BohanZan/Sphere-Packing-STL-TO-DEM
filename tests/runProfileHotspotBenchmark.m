function result = runProfileHotspotBenchmark(codeDirectory, label, doProfile, meshReader)
%RUNPROFILEHOTSPOTBENCHMARK Reproduce the current 50-sphere workload, seeded.
project = fileparts(fileparts(mfilename('fullpath')));
output = fullfile(project, 'tmp', 'perf_20260905');
oldDirectory = pwd;
oldPath = path;
cleanup = onCleanup(@() restoreEnvironment(oldDirectory, oldPath)); %#ok<NASGU>
addpath(project);
cd(codeDirectory);
clear spBuildContext spExactPointInside spRelax spPointInside spOrientInwardNormals
assert(strcmpi(strrep(fileparts(which('spBuildContext')),'\','/'), ...
    strrep(codeDirectory,'\','/')));
% Use the project's importer: native stlread rounds this ASCII STL differently.
% Import is outside timing so both variants still measure the core workflow.
if nargin<4, meshReader='project'; end
filename=fullfile(project,'inputs','greatBudda','greatBudda.stl');
if strcmp(meshReader,'project')
    [vertices,faces]=stlRead(filename);
    model=struct('vertices',vertices,'faces',faces);
else
    mesh=stlread(filename);
    model=struct('vertices',mesh.Points,'faces',mesh.ConnectivityList);
end
options = struct('gravity',[0 0 -1], 'density',1, ...
    'maxCompressionSweeps',200, 'compressionTolerance',1e-7, ...
    'shakeSweeps',5, 'maxRefillPasses',6, 'coordinateFrame','world', ...
    'randomSeed',42, 'outputDirectory','');
profile off
profile clear
if doProfile
    profile on -detail builtin -nohistory -timer performance
end
started = tic;
[assembly,masses,volume,inertia,report] = spawnSpheres( ...
    model, 0.625*ones(50,1), 60, 0, options);
elapsed = toc(started);
profile off
randomState = rng;
result = struct('label',label,'seconds',elapsed,'profiled',doProfile, ...
    'assembly',assembly,'masses',masses,'volume',volume, ...
    'inertia',inertia,'report',report,'randomState',randomState, ...
    'matlabVersion',version,'source',codeDirectory,'importer',meshReader);
if doProfile
    result.profile = profile('info');
    entries = result.profile.FunctionTable;
    [~,order] = sort([entries.TotalTime], 'descend');
    for k = 1:min(16,numel(order))
        item = entries(order(k));
        fprintf('PROFILE %s calls=%d total=%.6f\n', ...
            item.FunctionName,item.NumCalls,item.TotalTime);
    end
end
save(fullfile(output,[label '.mat']),'result','-v7.3');
fprintf('BENCHMARK %s seconds=%.6f accepted=%d\n',label,elapsed,report.acceptedCount);
end

function restoreEnvironment(directory, originalPath)
profile off
cd(directory);
path(originalPath);
end
