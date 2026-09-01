function index = spCellIndex(context, point)
%SPCELLINDEX Return the clamped 3-D spatial-hash index of POINT.
%Convert physical coordinates to one-based cell indices and clamp boundaries.
index = floor((point - context.lower) ./ context.cellSize) + 1;
index = min(max(index, 1), context.cellCount);
end
