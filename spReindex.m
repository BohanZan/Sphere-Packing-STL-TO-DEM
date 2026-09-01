function state = spReindex(context, state)
%SPREINDEX Rebuild the sphere spatial hash after a relaxation sweep.
%All sphere centres may have moved, so incremental updates are not reliable.
state.sphereCells=containers.Map('KeyType','char','ValueType','any');

%Insert every current sphere at its post-relaxation cell index.
for id=1:state.count
 idx=spCellIndex(context,state.centres(id,:));
 spHashInsert(state.sphereCells,idx,id);
end
end
