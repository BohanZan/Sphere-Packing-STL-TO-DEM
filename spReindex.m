function state = spReindex(context, state)
state.sphereCells=containers.Map('KeyType','char','ValueType','any');
for id=1:state.count
 idx=spGeometry('cell',context,state.centres(id,:));
 spGeometry('add',context,state.sphereCells,idx,id);
end
end
