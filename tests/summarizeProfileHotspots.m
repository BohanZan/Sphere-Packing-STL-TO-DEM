function summary=summarizeProfileHotspots
%SUMMARIZEPROFILEHOTSPOTS Verify output equivalence and print comparable timings.
project=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(project,'tmp','perf_20260905');
baseline=load(fullfile(folder,'v0_plain1.mat'));
projectBaseline=load(fullfile(folder,'v0_ascii_plain1.mat'));
files=dir(fullfile(folder,'v*_*.mat'));
rows=cell(0,3);
fields={'assembly','masses','volume','inertia','report','randomState'};
for id=1:numel(files)
    data=load(fullfile(folder,files(id).name));
    if ~isfield(data,'result'), continue; end
    expected=baseline.result;
    if isfield(data.result,'importer') && strcmp(data.result.importer,'project')
        expected=projectBaseline.result;
    end
    for k=1:numel(fields)
        field=fields{k};
        assert(isequaln(data.result.(field),expected.(field)), ...
            'Baseline mismatch: %s / %s',files(id).name,field);
    end
    rows(end+1,:)={data.result.label,data.result.seconds,data.result.profiled}; %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'Version','Seconds','Profiled'});
disp(summary);
writetable(summary,fullfile(folder,'timings.csv'));
for label={'v0_profile','v2_profile','v3_profile','v0_ascii_profile','v3_ascii_profile'}
    file=fullfile(folder,[label{1} '.mat']);
    if ~isfile(file), continue; end
    data=load(file);
    entries=data.result.profile.FunctionTable;
    fprintf('\n%s\n',label{1});
    targets={'spExactPointInside','spBuildContext','spRelax', ...
        'spRelax>firstContactDistance/sphereHitsCached', ...
        'spRelax>firstContactDistance/pointInsideCached'};
    for k=1:numel(targets)
        ix=find(strcmp({entries.FunctionName},targets{k}),1);
        if isempty(ix), continue; end
        entry=entries(ix);
        childTime=0;
        if ~isempty(entry.Children), childTime=sum([entry.Children.TotalTime]); end
        fprintf('%s calls=%d total=%.6f self=%.6f\n', ...
            entry.FunctionName,entry.NumCalls,entry.TotalTime,entry.TotalTime-childTime);
        lines=entry.ExecutedLines;
        [~,order]=sort(lines(:,3),'descend');
        disp(lines(order(1:min(4,numel(order))),:));
    end
end
fprintf('ALL_RUNS_EXACTLY_MATCH_BASELINE\n');
end
