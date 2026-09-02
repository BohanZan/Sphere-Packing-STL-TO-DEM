function index = spOccupancyCellIndex(occupancy, point)
%SPOCCUPANCYCELLINDEX Map a point to one clamped one-based occupancy index.
index = floor((point - occupancy.lower) ./ occupancy.cellSize) + 1;
index = min(max(index, 1), occupancy.cellCount);
end
