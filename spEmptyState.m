function state = spEmptyState(context)
state=struct('centres',zeros(0,3),'radii',zeros(0,1),'count',0,...
    'sphereCells',{cell(context.cellCount(1),context.cellCount(2),context.cellCount(3))});
end
