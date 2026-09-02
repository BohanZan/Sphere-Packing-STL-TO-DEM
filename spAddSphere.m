function state = spAddSphere(context, state, centre, radius)
%SPADDSPHERE Append one accepted sphere and index it in the sparse hash.
%Store centres and radii in insertion order for stable DEM output.
if state.count >= size(state.centres,1)
    error('SpherePacking:StateCapacityExceeded', ...
        'Cannot add a sphere beyond the preallocated state capacity.');
end
state.count=state.count+1; state.centres(state.count,:)=centre; state.radii(state.count,1)=radius;

%Register the new sphere under its current 3-D spatial cell.
idx=spCellIndex(context,centre);
state.cellIndices(state.count,:)=idx;
spHashInsert(state.sphereCells,idx,state.count);
end
