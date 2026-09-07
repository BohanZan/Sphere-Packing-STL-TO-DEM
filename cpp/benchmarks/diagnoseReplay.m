function diagnoseReplay
% Diagnostic only: retain production initial placement and contact solver.
project=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(project);addpath(project);
out=fullfile(project,'tmp');if ~isfolder(out),mkdir(out);end
cache=fullfile(out,'diagnose_context.mat');
if isfile(cache),load(cache,'context');else
 context=spBuildContext(fullfile(project,'inputs','greatBudda','greatBudda.stl'),.625,0,1e-9);
 context.gravityFrame=spGravityFrame(context,[0 0 -1]);save(cache,'context');
end
options=struct('maxAttempts',60,'maxCompressionSweeps',200,'shakeSweeps',5,'compressionTolerance',1e-7);
rng(42,'twister');state=spEmptyState(context,50);[state,next]=spInitialPlacement(context,state,.625*ones(50,1),1,options);
assert(next==51);fid=fopen(fullfile(out,'matlab_moves.csv'),'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
state=runRelax(context,state,options,true,fid);state=runRelax(context,state,options,false,fid);
writematrix(state.centres(1:50,:),fullfile(out,'matlab_final.csv'));
end
function state=runRelax(context,state,o,globalPass,fid)
f=context.gravityFrame;g=f.direction;weights=state.radii(1:state.count).^3;
previous=sum(((state.centres(1:state.count,:)-f.origin)*f.up.').*weights);
for sweep=1:o.maxCompressionSweeps
 if globalPass,shakes=0;else,shakes=o.shakeSweeps;end
 for shake=1:shakes+1
  for id=1:state.count
   if shake<=shakes
    u=randn(1,3);u=u-dot(u,g)*g;
    if norm(u)<eps,[~,axis]=min(abs(g));u=zeros(1,3);u(axis)=1;u=u-dot(u,g)*g;end
    u=u/norm(u);
   else,u=g;end
   p=state.centres(id,:);d=spFirstContactDistance(context,state,id,u);
   if globalPass && sweep==1 && id==7
    detail=[];
    for face=1:size(context.faces,1)
     tri=context.vertices(context.faces(face,:),:);hit=spTriangleContactDistance(p,state.radii(id),u,tri,d+1e-8,context.tolerance);
     if isfinite(hit)&&abs(hit-d)<1e-8
      normal=cross(tri(2,:)-tri(1,:),tri(3,:)-tri(1,:));len=norm(normal);n=normal/len;
      detail(end+1,:)=[face hit normal len n dot(p-tri(1,:),n) dot(u,n)]; %#ok<AGROW>
     end
    end
    df=fopen('tmp/matlab_first_contact.csv','w');fprintf(df,[repmat('%.17g,',1,10),'%.17g\n'],detail.');fclose(df);
   end
   fprintf(fid,'%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n',[globalPass,sweep,shake,id,p,u,d]);
   if d>0,state.centres(id,:)=p+d*u;state=spReindex(context,state,id);end
  end
 end
 energy=sum(((state.centres(1:state.count,:)-f.origin)*f.up.').*weights);
 if abs(energy-previous)<=o.compressionTolerance*max(abs(previous),realmin),return;end
 previous=energy;
end
end
