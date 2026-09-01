function state = spAddSphere(context, state, centre, radius)
%SPADDSPHERE Append one accepted sphere and index it in the sparse hash.
%Store centres and radii in insertion order for stable DEM output.
state.count=state.count+1; state.centres(state.count,:)=centre; state.radii(state.count,1)=radius;

%Register the new sphere under its current 3-D spatial cell.
idx=spCellIndex(context,centre);
spHashInsert(state.sphereCells,idx,state.count);
end
