function outputFiles = spWriteCsv(~, assembly, masses, totalVolume, inertia, report, options, context, state)
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
cellsFile = fullfile(options.outputDirectory, prefix + "_occupied_cells.csv");
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
spWriteOccupiedCells(cellsFile, context, state);
outputFiles = {char(sphereFile), char(summaryFile), char(cellsFile)};
end

function spWriteOccupiedCells(fileName, context, state)
cellKeys = unique([keys(context.triangleCells), keys(state.sphereCells)]);
rows = zeros(numel(cellKeys), 11);
for row = 1:numel(cellKeys)
    index = sscanf(cellKeys{row}, '%d,%d,%d').';
    lower = context.lower + (index - 1) .* context.cellSize;
    upper = min(lower + context.cellSize, context.upper);
    if isKey(state.sphereCells, cellKeys{row}), sphereCount = numel(state.sphereCells(cellKeys{row})); else, sphereCount = 0; end
    if isKey(context.triangleCells, cellKeys{row}), triangleCount = numel(context.triangleCells(cellKeys{row})); else, triangleCount = 0; end
    rows(row,:) = [index lower(1) upper(1) lower(2) upper(2) lower(3) upper(3) sphereCount triangleCount];
end
cellTable = array2table(rows, 'VariableNames', ...
    {'ix','iy','iz','xmin','xmax','ymin','ymax','zmin','zmax','sphere_count','triangle_count'});
writetable(cellTable, fileName, 'Delimiter', ',');
end
