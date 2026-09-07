function results=runCppParity(executable,outputDirectory)
%RUNCPPPARITY Replay MATLAB draws through the current C++ command-line program.
%Compare termination, draw counts, mass properties and all four CSV products.
%Use short shake fixtures; repeated interacting shakes amplify native norm
%roundoff (see CPP_ALIGNMENT.md) and are not a bitwise trajectory contract.
%Generated meshes, traces and reports remain in a temporary output directory.
here=fileparts(mfilename('fullpath')); root=fileparts(here);
cpp=fullfile(fileparts(root),'SpherePacking-cpp');
if nargin<1, executable=fullfile(cpp,'build','release','sphere_packing.exe'); end
if nargin<2, outputDirectory=tempname; end
assert(isfile(executable),'Build the C++ sphere_packing executable first.');
oldPath=path; oldRng=rng;
cleanup=onCleanup(@() restoreState(oldPath,oldRng)); %#ok<NASGU>
addpath(root,fullfile(here,'helpers'));
mkdir(outputDirectory);
fprintf('CPP_PARITY_OUTPUT %s\n',outputDirectory);
names={'box_batch_shake','box_layer_gravity','box_all_gravity','oblique','refill','buddha_gravity','oblique_shake'};
results=struct([]);
for scenario=1:numel(names)
 name=names{scenario}; directory=fullfile(outputDirectory,name); mkdir(directory);
 model=physicsBoxMesh;
 radii=.35*ones(16,1); attempts=3; buffer=0;
 options=struct('randomSeed',42,'maxCompressionSweeps',Inf,'shakeSweeps',1, ...
  'maxRefillPasses',2,'compressionTolerance',1e-7,'gravity',[0 0 -1], ...
  'coordinateFrame','world','density',2.5,'occupancyAcceleration',true, ...
  'outputDirectory',fullfile(directory,'matlab'),'outputPrefix','packing', ...
  'initialRelaxationScope','layer');
 if scenario<=3
  scopes={'batch','layer','all'}; options.initialRelaxationScope=scopes{scenario};
  if scenario>1, options.shakeSweeps=0; end
 elseif scenario==4
  options.gravity=[1 2 -3]; options.shakeSweeps=0;
  options.coordinateFrame='center_of_mass';
  radii=linspace(.2,.35,12).'; attempts=20;
 elseif scenario==5
  model=physicsBoxMesh([0 0 0],[2 2 2]); radii=.45*ones(20,1);
  attempts=8; options.shakeSweeps=0; options.gravity=[-1 0 0];
 elseif scenario==6
  meshPath=fullfile(root,'inputs','greatBudda','greatBudda.stl');
  [vertices,faces]=stlRead(meshPath); model=struct('vertices',vertices,'faces',faces);
  radii=.625*ones(12,1); attempts=60; options.shakeSweeps=0;
 elseif scenario==7
  options.gravity=[1 2 -3]; radii=.25*ones(2,1); attempts=20;
 end
 %Serialize the same triangle order and double coordinates for both readers.
 meshPath=fullfile(directory,'mesh.stl'); writeMesh(meshPath,model);
 radiiPath=fullfile(directory,'radii.txt'); writeNumbers(radiiPath,radii);
 global SP_CPP_REPLAY
 SP_CPP_REPLAY=struct('record',true,'uniforms',[],'normals',[],'uniformCount',0,'normalCount',0);
 warningState=warning('off','MATLAB:dispatcher:nameConflict');
 addpath(fullfile(here,'helpers','cppReplay'),'-begin'); warning(warningState);
 assembly=[]; masses=[]; volume=0; inertia=[]; report=struct;
 evalc('[assembly,masses,volume,inertia,report]=spawnSpheres(model,radii,attempts,buffer,options);');
 rmpath(fullfile(here,'helpers','cppReplay'));
 uniformPath=fullfile(directory,'uniforms.txt'); normalPath=fullfile(directory,'normals.txt');
 writeNumbers(uniformPath,SP_CPP_REPLAY.uniforms); writeNumbers(normalPath,SP_CPP_REPLAY.normals);
 cppDirectory=fullfile(directory,'cpp'); reportPath=fullfile(directory,'cpp.json');
 command=sprintf(['"%s" --model "%s" --radii "%s" --attempts %d --buffer %.17g ' ...
  '--shake-sweeps %d --refill-passes %d --compression-tolerance %.17g ' ...
  '--gravity %.17g %.17g %.17g --density %.17g --frame %s --initial-relaxation %s ' ...
  '--output "%s" --prefix packing --report "%s" --uniform-tape "%s" --normal-tape "%s" --quiet'], ...
  executable,meshPath,radiiPath,attempts,buffer,options.shakeSweeps,options.maxRefillPasses, ...
  options.compressionTolerance,options.gravity,options.density,options.coordinateFrame, ...
  options.initialRelaxationScope,cppDirectory,reportPath,uniformPath,normalPath);
 [status,output]=system(command);
 assert(status==0,'C++ replay failed for %s: %s',name,output);
 reference=jsondecode(fileread(reportPath));
 fields={'requestedCount','acceptedCount','unplacedCount','nextUnplacedRadiusIndex', ...
  'initialFailures','refillPasses','stopReason','capacityWarning','initialRelaxationScope','coordinateFrame'};
 for k=1:numel(fields)
  field=fields{k}; assert(isequal(reference.(field),report.(field)),'Report mismatch: %s/%s',name,field);
 end
 assert(reference.uniformDraws==SP_CPP_REPLAY.uniformCount && reference.normalDraws==SP_CPP_REPLAY.normalCount, ...
  'Random consumption differs in %s.',name);
 numericFields={'totalMass','stlVolume','boundingBoxDimensions','sphereAssemblyVolume','centreOfMass', ...
  'coordinateShift','centreOfMassAfterShift'};
 for k=1:numel(numericFields)
  field=numericFields{k}; assertNear(reference.(field),report.(field),name);
 end
 assertNear(reference.totalVolume,volume,name);
 assertNear(reshape(reference.inertia,3,3).',inertia,name);
 suffixes={'spheres','summary','grid_points','grid_hexahedra'};
 for k=1:numel(suffixes)
  actual=readtable(report.outputFiles{k});
  expected=readtable(fullfile(cppDirectory,['packing_' suffixes{k} '.csv']));
  assert(isequal(actual.Properties.VariableNames,expected.Properties.VariableNames),'CSV headers differ.');
  assert(isequal(size(actual),size(expected)),'CSV dimensions differ in %s/%s.',name,suffixes{k});
  for column=1:width(actual)
   a=actual{:,column}; b=expected{:,column};
   if isnumeric(a) || islogical(a), assertNear(a,b,name);
   else, assert(isequal(a,b),'CSV text differs in %s.',name); end
  end
 end
 spheres=readtable(fullfile(cppDirectory,'packing_spheres.csv'));
 assertNear([spheres.x spheres.y spheres.z spheres.radius].',assembly,name);
 assertNear(spheres.mass,masses,name);
 delta=abs([spheres.x spheres.y spheres.z].'-assembly(1:3,:));
 result=struct('scenario',name,'acceptedCount',report.acceptedCount,'initialFailures',report.initialFailures, ...
  'uniformDraws',reference.uniformDraws,'normalDraws',reference.normalDraws, ...
  'maximumCoordinateDifference',max([0;delta(:)]),'passed',true);
 results=[results;result]; %#ok<AGROW>
 fprintf('CPP_PARITY %s accepted=%d draws=%d/%d max_coordinate_error=%.4g PASS\n', ...
  name,result.acceptedCount,result.uniformDraws,result.normalDraws,result.maximumCoordinateDifference);
 SP_CPP_REPLAY=[];
end
fid=fopen(fullfile(outputDirectory,'parity_results.json'),'w'); assert(fid>=0);
fileCleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(results,PrettyPrint=true));
fprintf('CPP_PARITY_OUTPUT %s\n',outputDirectory);
end

function restoreState(oldPath,oldRng)
global SP_CPP_REPLAY
path(oldPath); rng(oldRng); SP_CPP_REPLAY=[];
end

function assertNear(actual,expected,scenario)
actual=double(actual(:)); expected=double(expected(:));
assert(numel(actual)==numel(expected),'Numeric dimensions differ in %s.',scenario);
assert(all(isfinite(actual)) && all(isfinite(expected)),'Nonfinite result in %s.',scenario);
assert(all(abs(actual-expected)<=1e-9*max(1,abs(expected))),'Numeric mismatch in %s (maximum error %.17g).', ...
 scenario,max([0;abs(actual-expected)]));
end
function writeNumbers(filename,values)
fid=fopen(filename,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%.17g\n',values);
end
function writeMesh(filename,model)
fid=fopen(filename,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'solid parity\n');
for id=1:size(model.faces,1)
 fprintf(fid,'facet normal 0 0 0\nouter loop\n');
 fprintf(fid,'vertex %.17g %.17g %.17g\n',model.vertices(model.faces(id,:),:).');
 fprintf(fid,'endloop\nendfacet\n');
end
fprintf(fid,'endsolid parity\n');
end
