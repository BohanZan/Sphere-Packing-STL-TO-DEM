function runMatlabBaseline(scenario,repetitions)
% Unprofiled compute timings exclude application startup, STL import and CSV IO.
% One separate traced run generates exact random replay and all four baseline CSVs.
if nargin<1, scenario='box'; end
if nargin<2, repetitions=3; end
here=fileparts(mfilename('fullpath')); project=fileparts(fileparts(here));
oldpath=path;oldpwd=pwd;cleanup=onCleanup(@() restore(oldpath,oldpwd)); %#ok<NASGU>
addpath(project);addpath(fullfile(project,'tests','helpers'));cd(project);
out=fullfile(here,'output',scenario);if ~isfolder(out),mkdir(out);end
options=struct('randomSeed',42,'maxCompressionSweeps',3,'shakeSweeps',2, ...
 'maxRefillPasses',3,'compressionTolerance',1e-7,'gravity',[0 0 -1], ...
 'density',1,'coordinateFrame','world','outputDirectory','', ...
 'occupancyAcceleration',true);
if strcmp(scenario,'box')
 model=physicsBoxMesh;radius=.25;count=24;attempts=40;
 meshPath=fullfile(out,'mesh.stl');writeAsciiMesh(meshPath,model);
 importSeconds=0;
elseif strcmp(scenario,'buddha') || strcmp(scenario,'buddha_gravity')
 meshPath=fullfile(project,'inputs','greatBudda','greatBudda.stl');
 clock=tic;[vertices,faces]=stlRead(meshPath);importSeconds=toc(clock);
 model=struct('vertices',vertices,'faces',faces);
 radius=.625;count=50;attempts=60;options.maxCompressionSweeps=200;
 options.shakeSweeps=5;options.maxRefillPasses=6;
 if strcmp(scenario,'buddha_gravity'),options.shakeSweeps=0;options.maxCompressionSweeps=3;end
elseif strcmp(scenario,'refill')
 model=physicsBoxMesh([0 0 0],[2 2 2]);radius=.45;count=20;attempts=12;
 options.maxCompressionSweeps=4;options.gravity=[-1 0 0];
 meshPath=fullfile(out,'mesh.stl');writeAsciiMesh(meshPath,model);importSeconds=0;
else,error('Unknown scenario');end
radii=radius*ones(count,1);
fprintf('BASELINE %s mesh=%d triangles count=%d radius=%.8g\n',scenario,size(model.faces,1),count,radius);
% Warm up MATLAB JIT and optional parallel pool before the measured runs.
evalc('[assembly,masses,volume,inertia,report]=spawnSpheres(model,radii,attempts,0,options);');
times=zeros(repetitions,1);
for iteration=1:repetitions
 clock=tic;
 evalc('[assembly,masses,volume,inertia,report]=spawnSpheres(model,radii,attempts,0,options);');
 times(iteration)=toc(clock);
 fprintf('MATLAB_TIME %s repetition=%d seconds=%.9g accepted=%d\n',scenario,iteration,times(iteration),report.acceptedCount);
end
baseline=struct('scenario',scenario,'meshPath',meshPath,'triangleCount',size(model.faces,1), ...
 'radius',radius,'count',count,'attempts',attempts,'options',options,'computeSeconds',times, ...
 'medianComputeSeconds',median(times),'importSeconds',importSeconds,'matlabVersion',version, ...
 'timingDefinition','one warmup; compute includes context/preprocessing and packing; excludes STL import, CSV IO, process startup; no trace/profiler', ...
 'assembly',assembly,'masses',masses,'volume',volume,'inertia',inertia,'report',report);
writeJson(fullfile(out,'matlab_timing.json'),baseline);
% Serialize source triangles at double precision to verify importer parity independently.
triangles=[model.vertices(model.faces(:,1),:) model.vertices(model.faces(:,2),:) model.vertices(model.faces(:,3),:)];
writematrix(triangles,fullfile(out,'triangles.csv'));
global SP_CPP_UNIFORMS SP_CPP_NORMALS
SP_CPP_UNIFORMS={};SP_CPP_NORMALS={};
addpath(fullfile(here,'trace'),'-begin');
options.outputDirectory=fullfile(out,'matlab');options.outputPrefix='packing';
evalc('[assembly,masses,volume,inertia,report]=spawnSpheres(model,radii,attempts,0,options);');
rmpath(fullfile(here,'trace'));
uniforms=vertcat(SP_CPP_UNIFORMS{:});normals=vertcat(SP_CPP_NORMALS{:});
writeValues(fullfile(out,'uniforms.txt'),uniforms);writeValues(fullfile(out,'normals.txt'),normals);
golden=struct('assembly',assembly,'masses',masses,'volume',volume,'inertia',inertia, ...
 'report',report,'uniformDraws',numel(uniforms),'normalDraws',numel(normals));
writeJson(fullfile(out,'matlab_golden.json'),golden);
clear global SP_CPP_UNIFORMS SP_CPP_NORMALS
fprintf('BASELINE_DONE %s median=%.9g draws=%d/%d\n',scenario,median(times),numel(uniforms),numel(normals));
end
function writeValues(filename,values)
fid=fopen(filename,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%.17g\n',values); %#ok<NASGU>
end
function writeJson(filename,value)
fid=fopen(filename,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(value,PrettyPrint=true)); %#ok<NASGU>
end
function writeAsciiMesh(filename,model)
fid=fopen(filename,'w');assert(fid>=0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'solid box\n');
for id=1:size(model.faces,1)
 fprintf(fid,'facet normal 0 0 0\nouter loop\n');
 fprintf(fid,'vertex %.17g %.17g %.17g\n',model.vertices(model.faces(id,:),:).');
 fprintf(fid,'endloop\nendfacet\n');
end
fprintf(fid,'endsolid box\n');
end
function restore(oldpath,oldpwd)
path(oldpath);cd(oldpwd);
end
