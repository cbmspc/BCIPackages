% function sensors = chan2lidx (ChanNames, Channels)
% Convert from channel names to logical index array for ChanNames
% Channels that are not found are not included
%
function lidx = chan2lidx (ChanNames, Channels)
idx = chan2idx(ChanNames, Channels, 1);
lidx = false(1, length(ChanNames));
lidx(idx) = 1;
