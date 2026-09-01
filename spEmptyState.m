function state = spEmptyState(~)
%SPEMPTYSTATE Initialise the mutable packing state before any insertion.
%The next progress threshold is advanced only after a successful placement.
state=struct('centres',zeros(0,3),'radii',zeros(0,1),'count',0,...
    'sphereCells',containers.Map('KeyType','char','ValueType','any'),...
    'nextProgressPercent',25);
end
