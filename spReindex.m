function state = spReindex(context, state, movedIds)
%SPREINDEX Rebuild all sphere cells or update only known moved sphere IDs.
if nargin < 3 || isempty(movedIds)
    state.sphereCells=containers.Map('KeyType','char','ValueType','any');
    for id=1:state.count
        idx=spCellIndex(context,state.centres(id,:));
        state.cellIndices(id,:)=idx;
        spHashInsert(state.sphereCells,idx,id);
    end
    return
end

for id=movedIds(:).'
    if id < 1 || id > state.count
        error('SpherePacking:InvalidMovedSphere', ...
            'Incremental reindex IDs must refer to accepted spheres.');
    end
    oldIndex=state.cellIndices(id,:);
    newIndex=spCellIndex(context,state.centres(id,:));
    if isequal(oldIndex,newIndex), continue; end

    oldKey=sprintf('%d,%d,%d',oldIndex(1),oldIndex(2),oldIndex(3));
    if ~isKey(state.sphereCells,oldKey)
        error('SpherePacking:MissingSphereCell', ...
            'The previous sparse cell is missing an accepted sphere.');
    end
    members=state.sphereCells(oldKey);
    position=find(members==id,1);
    if isempty(position)
        error('SpherePacking:MissingSphereIndex', ...
            'The previous sparse cell does not contain the accepted sphere.');
    end
    members(position)=[];
    if isempty(members)
        remove(state.sphereCells,oldKey);
    else
        state.sphereCells(oldKey)=members;
    end
    spHashInsert(state.sphereCells,newIndex,id);
    state.cellIndices(id,:)=newIndex;
end
end
