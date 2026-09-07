function occupancy = spBuildOccupancyGrid(context, maxRadius, options)
%SPBUILDOCCUPANCYGRID Conservatively label exact-inside shortcut cells.
% Labels are 0 outside, 1 inside, and 2 margin (which always falls back
% to the exact ray-parity test).  The temporary value 3 is never returned.

width = options.cellSize;
if isempty(width) || width==0
    width = max(maxRadius / 2, 32 * context.tolerance);
end
count = max(1, ceil((context.upper - context.lower) ./ width));
cellTotal = prod(double(count));
while cellTotal > options.maxCells
    width = width * ceil((cellTotal / options.maxCells)^(1/3));
    if ~isfinite(width), error('SpherePacking:OccupancyOverflow','Occupancy spacing overflow.'); end
    count = max(1, ceil((context.upper - context.lower) ./ width));
    cellTotal = prod(double(count));
end

occupancy = struct('enabled', true, 'lower', context.lower, 'cellSize', width, ...
    'cellCount', count, 'labels', []);
labels = 3 * ones(count, 'uint8');
margin = false(count);
% Keep the original AABB margin: tolerant ray parity can change at projected
% seams far from the physical surface. Surface-only components are unsafe.
halo = max(context.tolerance, 8 * eps(max(abs([context.lower context.upper]))));
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
%Visit linear seeds in X-fastest order, exactly as the C++ occupancy array.
for seed = 1:cellTotal
    if labels(seed)~=3, continue; end
    [ix,iy,iz]=ind2sub(count,seed);
    if any([ix iy iz]==1 | [ix iy iz]==count)
        floodComponent(seed,0);
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

    function componentLength = floodComponent(seed, value)
        head = 1;
        componentLength = 1;
        fifo(1) = uint32(seed);
        labels(seed) = value;
        %Expand one FIFO entry at a time: X-/X+, Y-/Y+, then Z-/Z+.
        %Do not sort neighbours; their discovery order fixes the traversal.
        strides=[1,count(1),count(1)*count(2)];
        while head <= componentLength
            current = double(fifo(head));
            head = head + 1;
            [cellX, cellY, cellZ] = ind2sub(count, current);
            index=[cellX cellY cellZ];
            for axis=1:3
                for direction=[-1 1]
                    coordinate=index(axis)+direction;
                    if coordinate<1 || coordinate>count(axis), continue; end
                    neighbour=current+direction*strides(axis);
                    if labels(neighbour)~=3, continue; end
                    labels(neighbour)=value;
                    componentLength=componentLength+1;
                    fifo(componentLength)=uint32(neighbour);
                end
            end
        end
    end

    function centre = cellCentre(linearIndex)
        [cellX, cellY, cellZ] = ind2sub(count, linearIndex);
        centre = occupancy.lower + ([cellX cellY cellZ] - 0.5) * occupancy.cellSize;
    end
end
