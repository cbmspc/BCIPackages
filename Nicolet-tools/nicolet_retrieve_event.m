% Retrieves EEG data of an event from a given subject/date/time 
% Uses the Nicolet .e files as source
%
% Primary Syntax: 
% [data, channames, fs] = nicolet_retrieve_event (Subject, EventTiming, RetrieveTiming, ChanNameRegexp, NotChanNameRegexp, DatabaseDir)
%   Subject = Subject name (part of the folder name before the @ sign)
%   EventTiming = [YYYY MM DD HH nn ss]
%   RetrieveTiming = [SecondsBefore SecondsAfter]
%   (Optional) ChanNameRegexp = regexp pattern for including channels.
%   (Optional) NotChanNameRegexp = regexp pattern for excluding channels
%   (Optional) DatabaseDir = location of database root directory
%   data = (time points x channel)
%   channames = cell of channel names
%   fs = sample rate in Hz
%
% If DatabaseDir is not specified, the function will ask you for the location on the first run
%
% Example:
%  [data, channames, fs] = nicolet_retrieve_event ('DOE_JOHN', [2016 12 25 15 00 08], [0 300]);
%   retrieves 5 minutes of data starting at 2016-12-25 15:00:08.
%

function [subdata, channames, samplerate] = nicolet_retrieve_event (Subject, EventTiming, RetrieveTiming, ChanNameRegexp, NotChanNameRegexp, DatabaseDir)
persistent BaseDir

if ~exist('ChanNameRegexp','var') || isempty(ChanNameRegexp)
    ChanNameRegexp = '.*';
end

if ~exist('NotChanNameRegexp','var') || isempty(NotChanNameRegexp)
    NotChanNameRegexp = 'NOT CONNECTED';
end

if ~exist('DatabaseDir','var') || isempty(DatabaseDir)
    DatabaseDir = '';
end

if ~isempty(DatabaseDir) && exist(DatabaseDir, 'dir')
    BaseDir = DatabaseDir;
end

if isempty(BaseDir)
    BaseDir = uigetdir('.', 'Select the database root directory');
end

EventDate = sprintf('%04i%02i%02i', EventTiming(1), EventTiming(2), EventTiming(3));
EventDateNum = datenum(sprintf('%04i-%02i-%02i %02i:%02i:%02i', EventTiming(1), EventTiming(2), EventTiming(3), EventTiming(4), EventTiming(5), EventTiming(6)));
SecondsBefore = RetrieveTiming(1);
SecondsAfter = RetrieveTiming(2);

if nargout < 3
    error('3 output arguments are required');
end

% Locate the folder
a = dir([BaseDir filesep Subject '@*']);
if isempty(a)
    error('There is no folder that matches the given subject');
end

for i = length(a):-1:1
    b = regexp(a(i).name, '^(\w+)@(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})$', 'tokens', 'once');
    FolderDateNum = datenum(sprintf('%s-%s-%s %s:%s:%s', b{2}, b{3}, b{4}, b{5}, b{6}, b{7}));
    if FolderDateNum < EventDateNum
        Folder = [BaseDir filesep a(i).name];
        break
    end
end

% Locate the segment
a = dir([Folder filesep '*.e']);
if isempty(a)
    error('There is no Nicolet .e file in this folder.');
end
NicoletFile = [Folder filesep a(1).name];
[~, SampleRates, ChannelNames, ~, ~, MiscFields] = nicolet_read (NicoletFile, 1);

iSegment = 0;
for i = length(MiscFields.segment):-1:1
    SegmentStartDateNum = datenum(MiscFields.segment(i).dateStr);
    SegmentDuration = MiscFields.segment(i).duration;
    SegmentSampleRates = SampleRates{i};
    SegmentChannelNames = ChannelNames{i};
    if SegmentStartDateNum <= EventDateNum - SecondsBefore/86400 && EventDateNum + SecondsAfter/86400 <= SegmentStartDateNum + SegmentDuration/86400
        % Event is contained entirely in this segment
        iSegment = i;
        break
    end
end

if ~iSegment
    error('No segment contains this event.');
end

EventStartSec = (EventDateNum - SegmentStartDateNum)*86400;
RetrieveStartSec = EventStartSec - SecondsBefore;
RetrieveEndSec = EventStartSec + SecondsAfter;
MedianSampleRate = median(SegmentSampleRates);
cid = ~cellfun(@isempty,regexpi(SegmentChannelNames, ChanNameRegexp, 'match', 'once'));
cid = cid & cellfun(@isempty,regexpi(SegmentChannelNames, NotChanNameRegexp, 'match', 'once')) & SegmentSampleRates==MedianSampleRate;
Data = nicolet_read (NicoletFile, 0, iSegment, RetrieveStartSec, RetrieveEndSec, SegmentChannelNames(cid));
subdata = [Data{iSegment}{cid}];
channames = SegmentChannelNames(cid);
samplerate = MedianSampleRate;
fprintf('Event data retrieved\n');
return

