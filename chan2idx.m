% function sensors = chan2idx (ChanNames, Channels)
% Convert from channel names to their index numbers in a list
% sensors = the first matching channel in ChanNames for each element of Channels
% allsensors = all matching channels in ChanNames for each element of Channels
% SWnonzero = instead of leaving a 0 placeholder, omit the non-matched channels
% SWregexp = use regular expression instead of exact match
function [sensors, allsensors] = chan2idx (ChanNames, Channels, SWnonzero, SWregexp)
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

if ~exist('SWregexp','var') || isempty(SWregexp) || ~isscalar(SWregexp)
    SWregexp = false;
end

N = length(Channels);
sensors = zeros(1,N);
allsensors = cell(1,N);
for i = 1:N
    if SWregexp
        I = find(~cellfun(@isempty,regexpi(ChanNames, Channels{i}, 'once')));
    else
        I = find(strcmpi(ChanNames, Channels{i}));
    end
    if ~isempty(I)
        sensors(i) = I(1);
        allsensors{i} = I;
    end
end

if exist('SWnonzero','var') && ~isempty(SWregexp) && isscalar(SWregexp) && SWnonzero
    sensors = sensors(sensors>0);
end

