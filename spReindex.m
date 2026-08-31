function state = spReindex(context, state)
state.sphereCells=cell(context.cellCount(1),context.cellCount(2),context.cellCount(3));
for id=1:state.count
 idx=spGeometry('cell',context,state.centres(id,:));
 state.sphereCells{idx(1),idx(2),idx(3)}(end+1)=id;
end
end
