function outputFiles = spWriteCsv(~, assembly, masses, totalVolume, inertia, report, options, context, state, coordinateShift)
%SPWRITECSV Write all persisted results as comma-separated files with headers.
%Return without writing when the caller intentionally omitted an output directory.
outputFiles = {};
if isempty(options.outputDirectory)
    return;
end

%Create the requested result directory and construct consistent file names.
prefix = string(options.outputPrefix);
if ~isscalar(prefix) || ismissing(prefix) || strlength(prefix)==0 || ...
 any(prefix==["." ".."]) || any(ismember(char(prefix),['/' char(92) ':' char(13) char(10)]))
    error('SpherePacking:InvalidOutputPrefix','Output prefix must be a filename without a directory.');
end
if ~isfolder(options.outputDirectory)
    mkdir(options.outputDirectory);
end
sphereFile = fullfile(options.outputDirectory, prefix + "_spheres.csv");
summaryFile = fullfile(options.outputDirectory, prefix + "_summary.csv");
gridPointsFile = fullfile(options.outputDirectory, prefix + "_grid_points.csv");
gridCellsFile = fullfile(options.outputDirectory, prefix + "_grid_hexahedra.csv");

%Preserve double precision with the same 17 significant digits as C++.
sphereRows=[(1:size(assembly,2)).',assembly.',2*assembly(4,:).',masses.'];
writeCsv(sphereFile,'id,x,y,z,radius,diameter,mass',sphereRows);
summaryRow=sprintf('%d,%d,%d,%s,%d,%.17g,%.17g,%.17g,%.17g\n', ...
    report.requestedCount,report.acceptedCount,report.unplacedCount, ...
    report.stopReason,report.capacityWarning,totalVolume,inertia(1,1),inertia(2,2),inertia(3,3));
writeCsv(summaryFile, ...
    'requested_count,accepted_count,unplaced_count,stop_reason,capacity_warning,total_volume,inertia_xx,inertia_yy,inertia_zz',summaryRow);

%Export occupied spatial cells as hexahedra in the selected coordinate frame.
spWriteOccupiedGrid(gridPointsFile, gridCellsFile, context, state, coordinateShift);
outputFiles = {char(sphereFile), char(summaryFile), char(gridPointsFile), char(gridCellsFile)};
end

function spWriteOccupiedGrid(pointsFile, cellsFile, context, state, coordinateShift)
%SPWRITEOCCUPIEDGRID Write occupied sparse hash cells as ParaView-ready hexahedra.
%Use cells containing either STL triangles or accepted sphere centres.
if isstruct(context.triangleCells)
    triangleKeys=context.triangleCells.occupiedKeys;
else
    triangleKeys=keys(context.triangleCells);
end
if ~iscell(triangleKeys), triangleKeys=num2cell(triangleKeys); end
cellKeys = [triangleKeys(:).', keys(state.sphereCells)];
% Preserve the previous lexical export order (not numeric key order).
indices = spCellIndices(cellKeys,context.cellCount);
labels=cell(size(indices,1),1);
for k=1:size(indices,1), labels{k}=sprintf('%d,%d,%d',indices(k,:)); end
[~,order]=unique(labels);
cellKeys=cellKeys(order); indices=indices(order,:);
pointRows = zeros(8*numel(cellKeys), 4);
cellRows = zeros(numel(cellKeys), 11);
for cellId = 1:numel(cellKeys)
    %Recover this cell's physical bounds and apply the common output shift.
    index = indices(cellId,:);
    lower = context.lower + (index - 1) .* context.cellSize - coordinateShift;
    upper = min(context.lower + (index - 1) .* context.cellSize + context.cellSize, context.upper) - coordinateShift;

    %Count local geometry and assign eight unique point identifiers.
    if isKey(state.sphereCells, cellKeys{cellId}), sphereCount = numel(state.sphereCells(cellKeys{cellId})); else, sphereCount = 0; end
    triangleCount=numel(spGridLookup(context.triangleCells,cellKeys(cellId)));
    pointIds = (cellId-1)*8 + (1:8);
    pointRows(pointIds,:) = [pointIds.', hexahedronCorners(lower, upper)];
    cellRows(cellId,:) = [cellId, pointIds, sphereCount, triangleCount];
end

%Persist the point list and hexahedral connectivity with descriptive headers.
writeCsv(pointsFile,'point_id,x,y,z',pointRows);
writeCsv(cellsFile,'cell_id,p1,p2,p3,p4,p5,p6,p7,p8,sphere_count,triangle_count',cellRows);
end

function writeCsv(filename,header,rows)
%WRITECSV Keep established headers and diagnose failed output streams.
[fid,message]=fopen(filename,'w');
if fid<0, error('SpherePacking:OutputFailed','Cannot open CSV: %s',message); end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',header);
if isnumeric(rows)
    format=[repmat('%.17g,',1,size(rows,2)-1),'%.17g\n'];
    fprintf(fid,format,rows.');
else
    fprintf(fid,'%s',rows);
end
[message,number]=ferror(fid);
if number~=0, error('SpherePacking:OutputFailed','Cannot write CSV: %s',message); end
end

function corners = hexahedronCorners(lower, upper)
%HEXAHEDRONCORNERS Return the eight axis-aligned vertices of one grid cell.
corners = [lower(1) lower(2) lower(3); upper(1) lower(2) lower(3); ...
    lower(1) upper(2) lower(3); upper(1) upper(2) lower(3); ...
    lower(1) lower(2) upper(3); upper(1) lower(2) upper(3); ...
    lower(1) upper(2) upper(3); upper(1) upper(2) upper(3)];
end
