% function sensors = chan2idx (ChanNames, Channels)
% Convert from channel names to their index numbers in a list
%
function [sensors, allsensors] = chan2idx (ChanNames, Channels, SWnonzero)
if ~iscell(ChanNames)
    if ischar(ChanNames)
        ChanNames = string_to_cell(ChanNames,' ,');
    elseif isnumeric(ChanNames)
        ChanNames = string_to_cell(sprintf('%i ', ChanNames), ' ');
    end
end
if ~iscell(Channels)
    if ischar(Channels)
        Channels = string_to_cell(Channels,' ,');
    elseif isnumeric(Channels)
        Channels = string_to_cell(sprintf('%i ', Channels), ' ');
    end
end
N = length(Channels);
sensors = zeros(1,N);
allsensors = cell(1,N);
for i = 1:N
    I = find(strcmpi(ChanNames,Channels{i}));
    if ~isempty(I)
        sensors(i) = I(1);
        allsensors{i} = I;
    end
end

if exist('SWnonzero','var') && SWnonzero
    sensors = sensors(sensors>0);
end

