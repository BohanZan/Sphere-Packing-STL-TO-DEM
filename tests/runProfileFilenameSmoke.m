function runProfileFilenameSmoke
%RUNPROFILEFILENAMESMOKE Exercise the real STL importer and CSV output path.
project=fileparts(fileparts(mfilename('fullpath')));
addpath(project);
cd(project);
folder=fullfile(project,'tmp','perf_20260905');
options=struct('gravity',[0 0 -1],'density',1, ...
    'maxCompressionSweeps',200,'compressionTolerance',1e-7, ...
    'shakeSweeps',5,'maxRefillPasses',6,'coordinateFrame','world', ...
    'randomSeed',42,'outputDirectory',fullfile(folder,'csv_smoke'), ...
    'outputPrefix','verified');
file=fullfile(project,'inputs','greatBudda','greatBudda.stl');
[assembly,masses,volume,inertia,report]=spawnSpheres(file,.625*ones(50,1),60,0,options);
baseline=load(fullfile(folder,'v0_ascii_plain1.mat'));
assert(isequaln(assembly,baseline.result.assembly));
assert(isequaln(masses,baseline.result.masses));
assert(isequaln(volume,baseline.result.volume));
assert(isequaln(inertia,baseline.result.inertia));
assert(all(cellfun(@isfile,report.outputFiles)));
rows=readtable(report.outputFiles{1});
assert(height(rows)==50);
assert(max(abs(rows.x-assembly(1,:).'))<1e-10);
assert(max(abs(rows.y-assembly(2,:).'))<1e-10);
assert(max(abs(rows.z-assembly(3,:).'))<1e-10);
save(fullfile(folder,'csv_smoke_validation.mat'),'assembly','report');
saved=load(fullfile(folder,'v3_tests.mat'));
fprintf('FINAL_TESTS passed=%d total=%d duration=%.6f\n', ...
    sum([saved.results.Passed]),numel(saved.results),sum([saved.results.Duration]));
fprintf('CSV_SMOKE_AND_FILENAME_PATH_PASS\n');
end
