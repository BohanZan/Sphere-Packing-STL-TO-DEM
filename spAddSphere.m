function state = spAddSphere(context, state, centre, radius)
state.count=state.count+1; state.centres(state.count,:)=centre; state.radii(state.count,1)=radius;
idx=spGeometry('cell',context,centre);
state.sphereCells{idx(1),idx(2),idx(3)}(end+1)=state.count;
end
