% function sensors = chan2idx (ChanNames, Channels)
% Convert from channel names to their index numbers in a list
% sensors = the first matching channel in ChanNames for each element of Channels
% allsensors = all matching channels in ChanNames for each element of Channels
% SWnonzero = instead of leaving a 0 placeholder, omit the non-matched channels
% SWregexp = use regular expression instead of exact match. 
%      Note: sensors will include all matching channels instead of only the
%            first matching channel. Nonmatching channels are not included
%            as placeholders.
% SWcasesensitive = case sensitive
% SWunique = remove duplicates while preserving order
%
function [sensors, allsensors] = chan2idx (ChanNames, Channels, SWnonzero, SWregexp, SWcasesensitive, SWunique)
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

if ~exist('SWnonzero','var') || ~isscalar(SWnonzero)
    SWnonzero = false;
end

if ~exist('SWregexp','var') || ~isscalar(SWregexp)
    SWregexp = false;
end

if ~exist('SWcasesensitive','var') || ~isscalar(SWcasesensitive)
    SWcasesensitive = false;
end

if ~exist('SWunique','var') || ~isscalar(SWunique)
    SWunique = false;
end

N = length(Channels);
sensors = zeros(1,N);
allsensors = cell(1,N);
for i = 1:N
    if SWregexp
        if SWcasesensitive
            I = find(~cellfun(@isempty,regexp(ChanNames, Channels{i}, 'once')));
        else
            I = find(~cellfun(@isempty,regexpi(ChanNames, Channels{i}, 'once')));
        end
    else
        if SWcasesensitive
            I = find(strcmp(ChanNames, Channels{i}));
        else
            I = find(strcmpi(ChanNames, Channels{i}));
        end
    end
    if ~isempty(I)
        sensors(i) = I(1);
        allsensors{i} = I;
    end
end

if SWnonzero
    sensors = sensors(sensors>0);
end

if SWregexp
    if isempty(allsensors)
        sensors = [];
    else
        sensors = [allsensors{:}];
    end
end

if SWunique
    sensors = unique(sensors, 'stable');
end

return

