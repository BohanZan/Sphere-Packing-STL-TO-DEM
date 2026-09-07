function inside = spPointInside(context, point)
%SPPOINTINSIDE Classify known occupancy cells or preserve exact parity.
if numel(point)~=3 || ~isreal(point) || any(~isfinite(point(:)))
 error('SpherePacking:InvalidPoint','Point coordinates must be finite and real.');
end
point=double(point(:).');
if isfield(context, 'occupancy') && context.occupancy.enabled
    index = spOccupancyCellIndex(context.occupancy, point);
    label = context.occupancy.labels(index(1), index(2), index(3));
    if label == 0
        inside = false;
        return
    elseif label == 1
        inside = true;
        return
    end
end
inside = spExactPointInside(context, point);
end
