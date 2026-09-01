function outputFiles = spWriteCsv(~, assembly, masses, totalVolume, inertia, report, options)
%SPWRITECSV Write all persisted results as comma-separated files with headers.
outputFiles = {};
if isempty(options.outputDirectory)
    return;
end
if ~isfolder(options.outputDirectory)
    mkdir(options.outputDirectory);
end
prefix = string(options.outputPrefix);
sphereFile = fullfile(options.outputDirectory, prefix + "_spheres.csv");
summaryFile = fullfile(options.outputDirectory, prefix + "_summary.csv");
sphereTable = table((1:size(assembly,2)).', assembly(1,:).', assembly(2,:).', ...
    assembly(3,:).', assembly(4,:).', 2*assembly(4,:).', masses.', ...
    'VariableNames', {'id','x','y','z','radius','diameter','mass'});
summaryTable = table(report.requestedCount, report.acceptedCount, report.unplacedCount, ...
    string(report.stopReason), report.capacityWarning, totalVolume, ...
    inertia(1,1), inertia(2,2), inertia(3,3), ...
    'VariableNames', {'requested_count','accepted_count','unplaced_count', ...
    'stop_reason','capacity_warning','total_volume','inertia_xx','inertia_yy','inertia_zz'});
writetable(sphereTable, sphereFile, 'Delimiter', ',');
writetable(summaryTable, summaryFile, 'Delimiter', ',');
outputFiles = {char(sphereFile), char(summaryFile)};
end
