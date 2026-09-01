function state = spAddSphere(context, state, centre, radius)
state.count=state.count+1; state.centres(state.count,:)=centre; state.radii(state.count,1)=radius;
idx=spCellIndex(context,centre);
spHashInsert(state.sphereCells,idx,state.count);
end
