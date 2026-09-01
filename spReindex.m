function state = spReindex(context, state)
state.sphereCells=containers.Map('KeyType','char','ValueType','any');
for id=1:state.count
 idx=spCellIndex(context,state.centres(id,:));
 spHashInsert(state.sphereCells,idx,id);
end
end
