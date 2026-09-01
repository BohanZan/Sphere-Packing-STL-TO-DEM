function state = spReportFillProgress(state, requestedCount)
%SPREPORTFILLPROGRESS Print each successful 25-percent filling milestone once.
while state.nextProgressPercent <= 100 && ...
        100 * state.count / requestedCount >= state.nextProgressPercent
    fprintf('Filling Progress: %d%% (%d / %d spheres placed)\n', ...
        state.nextProgressPercent, state.count, requestedCount);
    state.nextProgressPercent = state.nextProgressPercent + 25;
end
end
