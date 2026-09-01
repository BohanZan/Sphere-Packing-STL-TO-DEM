function outputFiles = spWriteCsv(~, assembly, masses, totalVolume, inertia, report, options, context, state, coordinateShift)
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
gridPointsFile = fullfile(options.outputDirectory, prefix + "_grid_points.csv");
gridCellsFile = fullfile(options.outputDirectory, prefix + "_grid_hexahedra.csv");
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
spWriteOccupiedGrid(gridPointsFile, gridCellsFile, context, state, coordinateShift);
outputFiles = {char(sphereFile), char(summaryFile), char(gridPointsFile), char(gridCellsFile)};
end

function spWriteOccupiedGrid(pointsFile, cellsFile, context, state, coordinateShift)
cellKeys = unique([keys(context.triangleCells), keys(state.sphereCells)]);
pointRows = zeros(8*numel(cellKeys), 4);
cellRows = zeros(numel(cellKeys), 11);
for cellId = 1:numel(cellKeys)
    index = sscanf(cellKeys{cellId}, '%d,%d,%d').';
    lower = context.lower + (index - 1) .* context.cellSize - coordinateShift;
    upper = min(context.lower + (index - 1) .* context.cellSize + context.cellSize, context.upper) - coordinateShift;
    if isKey(state.sphereCells, cellKeys{cellId}), sphereCount = numel(state.sphereCells(cellKeys{cellId})); else, sphereCount = 0; end
    if isKey(context.triangleCells, cellKeys{cellId}), triangleCount = numel(context.triangleCells(cellKeys{cellId})); else, triangleCount = 0; end
    pointIds = (cellId-1)*8 + (1:8);
    pointRows(pointIds,:) = [pointIds.', hexahedronCorners(lower, upper)];
    cellRows(cellId,:) = [cellId, pointIds, sphereCount, triangleCount];
end
pointTable = array2table(pointRows, 'VariableNames', {'point_id','x','y','z'});
cellTable = array2table(cellRows, 'VariableNames', ...
    {'cell_id','p1','p2','p3','p4','p5','p6','p7','p8','sphere_count','triangle_count'});
writetable(pointTable, pointsFile, 'Delimiter', ',');
writetable(cellTable, cellsFile, 'Delimiter', ',');
end

function corners = hexahedronCorners(lower, upper)
corners = [lower(1) lower(2) lower(3); upper(1) lower(2) lower(3); ...
    lower(1) upper(2) lower(3); upper(1) upper(2) lower(3); ...
    lower(1) lower(2) upper(3); upper(1) lower(2) upper(3); ...
    lower(1) upper(2) upper(3); upper(1) upper(2) upper(3)];
end
