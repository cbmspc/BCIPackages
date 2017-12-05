% This reads/imports the Nicolet .e file
% The format of the output Data is:
% Data{segment #}{ch} = actual signals
% SampleRates{segment #}(ch #) = sample rate for this channel
% ChannelNames{segment #}(ch #) = channel name for this channel
% Scales{segment #}(ch #) = scale for this channel
%
% Optional input arguments:
% SWskipdata == 1, no data is retrieved. Useful for only getting metadata.
% Segment: specify segments to retrieve. Example: [1 2]. Default = all
% StartSec: specify start time in seconds. Default = 0
% EndSec: specify end time. Default = end of each segment
% Channels: specify channel names to retrieve (not case sensitive) as cell array. Default = all
% Segment must be specified if any of the next 3 arguments are specified.
%
% written by Po T Wang
% uses Nicolet reader by jwagenaar at https://github.com/ieeg-portal/Nicolet-Reader

function [Data, SampleRates, ChannelNames, Scales, eventMarkers, MiscFields] = nicolet_read (File, SWskipdata, Segments, StartSec, EndSec, Channels)

if ~exist('SWskipdata', 'var')
    SWskipdata = 0;
end

if ~exist('Segments', 'var')
    Segments = [];
end

if ~exist('StartSec', 'var') || isempty(StartSec)
    StartSec = 0;
end

if ~exist('EndSec', 'var') || isempty(EndSec)
    EndSec = inf;
end

if ~exist('Channels', 'var') || isempty(Channels) || ~iscell(Channels)
    Channels = {};
end


OBJ = NicoletFile(File);
Nsegment = length(OBJ.segments);
for seg = Nsegment:-1:1
    % Note: Each channel can have a different sample rate.
    SampleRates{seg} = OBJ.segments(seg).samplingRate; % in Hz
    ChannelNames{seg} = OBJ.segments(seg).chName;
    Durations(seg) = OBJ.segments(seg).duration; % in seconds
    Scales{seg} = OBJ.segments(seg).scale;
end

if SWskipdata
    Data = {};
else
    for seg = Nsegment:-1:1
        for ch = 1:length(ChannelNames{seg})
            Data{seg}{ch} = [];
        end
    end
    for seg = 1:Nsegment
        % Note: Each channel can have a different sample rate.
        d = Durations(seg);
        Fs = SampleRates{seg};
        
        if ~SWskipdata
            if isempty(Segments)
                % Retrieve everything
                for ch = 1:length(ChannelNames{seg})
                    Data{seg}{ch} = getdata(OBJ, seg, [1 d*Fs(ch)], ch);
                end
            else
                % Retrieve only specified segment, time span, and channels
                if ismember(seg, Segments)
                    for ch = 1:length(ChannelNames{seg})
                        if ~isempty(Channels) && ~any(strcmpi(ChannelNames{seg}{ch}, Channels))
                            continue
                        end
                        SegmentStartIndex = round(StartSec*Fs(ch)+1);
                        if SegmentStartIndex < 1
                            SegmentStartIndex = 1;
                        end
                        SegmentEndIndex = round(EndSec*Fs(ch)+1);
                        if SegmentEndIndex > d*Fs(ch)
                            SegmentEndIndex = d*Fs(ch);
                        end
                        if SegmentStartIndex > d*Fs(ch) || SegmentEndIndex < 1 || SegmentEndIndex < SegmentStartIndex
                            % Nothing to retrieve
                            continue
                        end
                        Data{seg}{ch} = getdata(OBJ, seg, [SegmentStartIndex SegmentEndIndex], ch);
                    end
                end
                
            end
        end
    end
end

eventMarkers = OBJ.eventMarkers;
MiscFields.sigInfo = OBJ.sigInfo;
MiscFields.tsInfo = OBJ.tsInfo;
MiscFields.chInfo = OBJ.chInfo;
SegmentFieldsToCapture = {
    'dateOLE'
    'dateStr'
    'startDate'
    'startTime'
    'duration'
    'chName'
    'refName'
    'samplingRate'
    'scale'
};
for i = 1:Nsegment
    for j = 1:length(SegmentFieldsToCapture)
        MiscFields.segment(i).(SegmentFieldsToCapture{j}) = OBJ.segments(i).(SegmentFieldsToCapture{j});
    end
end
delete(OBJ);
return