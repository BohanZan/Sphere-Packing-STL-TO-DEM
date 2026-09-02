function occupancy = spBuildOccupancyGrid(context, maxRadius, options)
%SPBUILDOCCUPANCYGRID Conservatively label exact-inside shortcut cells.
% Labels are 0 outside, 1 inside, and 2 margin (which always falls back
% to the exact ray-parity test).  The temporary value 3 is never returned.

width = options.cellSize;
if isempty(width)
    width = max(maxRadius / 2, 32 * context.tolerance);
end
count = max(1, ceil((context.upper - context.lower) ./ width));
cellTotal = prod(double(count));
while cellTotal > options.maxCells
    width = width * ceil((cellTotal / options.maxCells)^(1/3));
    count = max(1, ceil((context.upper - context.lower) ./ width));
    cellTotal = prod(double(count));
end

occupancy = struct('enabled', true, 'lower', context.lower, 'cellSize', width, ...
    'cellCount', count, 'labels', []);
labels = 3 * ones(count, 'uint8');
margin = false(count);
halo = max(context.tolerance, 8 * eps(max(abs(context.vertices(:)))));
for triangleId = 1:size(context.faces, 1)
    triangle = context.vertices(context.faces(triangleId,:), :);
    lowerIndex = spOccupancyCellIndex(occupancy, min(triangle, [], 1) - width - halo);
    upperIndex = spOccupancyCellIndex(occupancy, max(triangle, [], 1) + width + halo);
    margin(lowerIndex(1):upperIndex(1), lowerIndex(2):upperIndex(2), ...
        lowerIndex(3):upperIndex(3)) = true;
end
labels(margin) = 2;

% One reusable FIFO labels 6-connected non-margin components without recursion.
fifo = zeros(cellTotal, 1, 'uint32');
for iy = 1:count(2)
    for iz = 1:count(3)
        floodBoundary(sub2ind(count, 1, iy, iz));
        floodBoundary(sub2ind(count, count(1), iy, iz));
    end
end
for ix = 1:count(1)
    for iz = 1:count(3)
        floodBoundary(sub2ind(count, ix, 1, iz));
        floodBoundary(sub2ind(count, ix, count(2), iz));
    end
end
for ix = 1:count(1)
    for iy = 1:count(2)
        floodBoundary(sub2ind(count, ix, iy, 1));
        floodBoundary(sub2ind(count, ix, iy, count(3)));
    end
end

for seed = 1:cellTotal
    if labels(seed) ~= 3, continue; end
    componentLength = floodComponent(seed, 0);
    representative = cellCentre(double(fifo(1)));
    if spExactPointInside(context, representative)
        labels(double(fifo(1:componentLength))) = 1;
    end
end
if any(labels(:) == 3)
    error('SpherePacking:OccupancyLabellingFailed', ...
        'Occupancy-grid construction left unclassified cells.');
end
occupancy.labels = labels;

    function floodBoundary(seed)
        if labels(seed) == 3
            floodComponent(seed, 0);
        end
    end

    function componentLength = floodComponent(seed, value)
        head = 1;
        componentLength = 1;
        fifo(1) = uint32(seed);
        labels(seed) = value;
        % Process small FIFO batches so large connected components do not
        % spend one MATLAB function call per neighbour.
        while head <= componentLength
            batchEnd = min(componentLength, head + 4095);
            current = double(fifo(head:batchEnd));
            head = batchEnd + 1;
            [cellX, cellY, cellZ] = ind2sub(count, current);
            neighbours = [current(cellX > 1) - 1; current(cellX < count(1)) + 1; ...
                current(cellY > 1) - count(1); current(cellY < count(2)) + count(1); ...
                current(cellZ > 1) - count(1)*count(2); ...
                current(cellZ < count(3)) + count(1)*count(2)];
            neighbours = unique(neighbours);
            neighbours = neighbours(labels(neighbours) == 3);
            if isempty(neighbours), continue; end
            labels(neighbours) = value;
            next = componentLength + (1:numel(neighbours));
            fifo(next) = uint32(neighbours);
            componentLength = next(end);
        end
    end

    function centre = cellCentre(linearIndex)
        [cellX, cellY, cellZ] = ind2sub(count, linearIndex);
        centre = occupancy.lower + ([cellX cellY cellZ] - 0.5) * occupancy.cellSize;
    end
end
