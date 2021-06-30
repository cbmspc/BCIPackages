
function [rawdata, SampleRate, ADCfactor, filteredvoltsignal, rawvoltsignal, suid] = nexus_scope (ChanNames, SampleRate, varargin)
mf = mfilename('fullpath');
a = dir([mf '.m']);
if ~isempty(a)
    md = a.date;
else
    md = 'date not available';
end
disp(['nexus_scope (' md ')']);
clear mf a md

%Warning: Turning off autospeed can cause data loss when program cannot
%process the incoming stream in real time
SWautospeed = 1;

% The following switches are only used for display:
SWcar = 0; % common average referencing
SWfilter = 0; % band pass frequency filter
HPF_Order = 2; % was 2.
LPF_Order = 8; % was 2.
SWnotch = 0; % power line notch filter
PowerLineFrequency = 60.00; % power line frequency in Hz
NotchFilter.order = 4; % Must be even
NotchFilter.qfactor = 10;
Funda = PowerLineFrequency;
for h = floor(SampleRate/2/Funda):-1:1
    % Create one for each harmonic
    d = fdesign.notch('N,F0,Q',NotchFilter.order,Funda*h/(SampleRate/2),NotchFilter.qfactor*h);
    notchHd{h} = design(d);
end
clear d h
SWnohardwaremode = 0; % No hardware mode (for viewing the scope window only)
fakecounter = 0;
Fcutoff = 30; % Filter cutoff (Hz) for low pass filter
Fcutoffhigh = 1; % Filter cutoff (Hz) for high pass filter
SWenvelope = 0; % Power envelope switch
Fenvelope = 1.5; % Enveloping filter
SWautoscale = 1; % auto-scale
SWovervoltagewarning = 1; % over-voltage warning
CSWmediandisp = [1 0];
SWzscore = 0;
AutoNumericChanNames = 0;  % When greater than 1, automatically swap between channel names and plug numbers every AutoNumericChanNames seconds.
SWmaximize = 1; % Maximize window on load.

AcqDuration = 10.0; % Do not change This is the acquisition buffer size
globalscale = 1;

% % Disallow wrap-around recording?
% SWnowrap = 0;

if nargin < 2
    SampleRate = 256;
end

if SampleRate < 128
    SampleRate = 128;
elseif SampleRate > 2048
    SampleRate = 2048;
end

if ~exist('SubjectName','var') || isempty(SubjectName)
    SubjectName = 'Subject Not Specified';
end

AnalysisDuration = 10.000;

% How far to look back to determine signal quality (seconds)
SigQualLookbackTime = 3.000;

% Peak-to-peak text readout look back time
TextPPLookbackTime = 1.500;

% Maximum sane voltage
MaxSaneVoltage = 0.30;

% Sometimes a channel's signal is too strong and spills into another
% channel's space. Set DispSpillAllow to 1.0 to cap it at each channel's
% ceiling/floor. Set it to higher to allow some spill.
DispSpillAllow = inf;

% How long to wait before stopping recording after trigger signal is lost. 
% This is only applicable when rawdata is specified as output.
%TrigLostStopWaitWarn = 5.000;
%TrigLostStopWait = 15.000;

TrigLostStopWaitWarn = 5.000;
TrigLostStopWait = inf;

% Turns yellow when peak-to-peak is over this voltage
OverVoltageLevel = 1400e-6; %100;

% Plot every N samples (saves time plotting less stuff)
PlotEveryNSample = 2;

% Allowable autoscale range (in volts)
%AutoScaleAllow = [1e-6 6.00];
AutoScaleAllow = [(6/2^22/20), 6.00];

% Allowable range for refresh period
RefreshPeriodAllow = [1/8, 8];

% Inactivity timeout auto exit
InactivityTimeoutSec = 3600;


%20100721
if nargin > 2
    disp('Command line arguments: ');
    for i = 1:2:nargin-2
        disp([varargin{i} ' = ' num2str(varargin{i+1})]);
        
        if ischar(varargin{i+1})
            eval([varargin{i} ' = ' '''' varargin{i+1} '''' ';']);
        else
            eval([varargin{i} ' = ' num2str(varargin{i+1}) ';']);
        end
    end
end

if ~exist('SWplot','var')
    SWplot = 1;
end

if ~SWplot
    SWautoscale = 0;
end

if ~exist('SWmemcheck','var')
    SWmemcheck = 0;
end

if ~exist('SWwiper','var')
    SWwiper = 1;
end

if RefreshPeriodAllow(1) < 1/SampleRate
    RefreshPeriodAllow = 1/SampleRate;
end

if TrigLostStopWaitWarn > TrigLostStopWait
    TrigLostStopWaitWarn = TrigLostStopWait;
end

% Notch filter parameters
%PowerlineFrequency = 60;
%[FiltcoefB,FiltcoefA] = butter(5,(PowerlineFrequency+[-3,3])/(SampleRate/2),'stop');
%[FiltcoefB,FiltcoefA] = butter(5,50/(SampleRate/2),'low');

%if SWfilter
[FiltcoefB, FiltcoefA, nlookback] = get_filter_coef(SampleRate, Fcutoffhigh, Fcutoff);
%end
%

[EnvelopeFiltcoefB, EnvelopeFiltcoefA, nlookback_env] = get_filter_coef(SampleRate, 0, Fenvelope);

BgdFigureColor = [0 0 0];
BgdFontColor = [0 0 0];
NormalFontColor = [0 1 0];
FailedFontColor = [1 0 0];
OverFontColor = [1 0.6 0];
HighestFontColor = [1 0.3 0];
AngryFontColor = [1 0 1];
BattLowBackgroundColor = [0 0 0];
BattLowEdgeColor = [1 1 0];
BattLowLineWidth = 2.5;
WaitFontColor = [0 1 1];
SavedFontColor = [0.3 0.3 1];
WarningNotRecordedString = 'NOT RECORDED';
WarningTriggerLostString = 'TRIGGER LOST';
WarningNoHardwareString = ['No connection to amplifiers.' 10 'Displayed signals are simulated.' 10 'NOT REAL DATA'];
WarningBatteryLowString = 'LOW BATTERY';
WarningDeviceErrorString = ['DEVICE ERROR' 10 'NO DATA RECEIVED'];
WarningNotRecordedColor1 = [1.0 0 0];
WarningNotRecordedColor2 = [0.5 0 0];
WarningNotRecordedColorChanged = 0;
% WarnedNotRecorded = 0;
% WarnedTriggerLost = 0;
% WarnedLowBattery = 0;
% WarnedDeviceError = 0;
% WarnedBufferFull = 0;
Warned.NotRecorded = 0;
Warned.TriggerLost = 0;
Warned.LowBattery = 0;
Warned.DeviceError = 0;
Warned.NoHardware = 0;
Warned.BufferFull = 0;
WarningBufferFullString = 'BUFFER ALMOST FULL';
WarnBufferFullBeforeSec = 120.000;
WarningNotRecordedFontSize = 36;
FontSize = 12; % was 10 for 32-chan
RestartNexusButton.normaltext = 'Restart Nexus';
RestartNexusButton.disabledtext = 'Cannot Restart';
RestartNexusButton.actiontext = 'Please Wait ..';
RestartNexusButton.confirmtext = 'Really Restart?';

RestartNexusConfirm = 0;

RefreshPeriod = 0.50; % starting refreshperiod.. will change automatically to best value
FetchExtraSec = 0.125; %
%MaxChanInRow = 34;
OriginalChanNames = ChanNames;
if length(ChanNames) > 34
    systemMode = 'synfi';
    AnalysisDuration = AnalysisDuration / 2;
    TwoColumnMode = 1;
else
    systemMode = 'fusbi';
    TwoColumnMode = 0;
end

% Line colors
kolors = 0.7*ones(67,1)*(1-BgdFigureColor);
NormalLineColor = 0.7*(1-BgdFigureColor);
OverLineColor = 0.7*OverFontColor;
AngryLineColor = 0.7*AngryFontColor;

% tmp = brighten(0.50);
% tmp2([1:2:64],:) = tmp([1:32],:);
% tmp2([2:2:64],:) = tmp([64:-1:33],:);

% Provide recording buffer using memory mapped file. 
% 63 minutes takes about 4 GiB at 2048 Hz
%BufferDurationSec = floor(63*60*2048/SampleRate);
BufferDurationSec = 63*60;


NexusChanNamesFileName = [getcccdatadir() filesep 'intermediate' filesep 'NexusChanNames.mat'];
tmp = dir(NexusChanNamesFileName);
if length(tmp) == 1
    NexusChanNamesFileDate = tmp.datenum;
else
    NexusChanNamesFileDate = 0;
end


if SWmemcheck
    rawdata = BufferDurationSec;
    return
end


if strcmpi(systemMode,'fusbi')
    MaxChan = 34;
    %METAchan = [33];
    METAchan = [33 34];
    STATUSchan = 33;
    DIAGchan = 34;
%     kolors(setdiff(1:64+length(METAchan),METAchan),:) = tmp2;
else
    MaxChan = 67;
    %METAchan = [33 66];
    METAchan = [33 66 67];
    STATUSchan = [33 66];
    DIAGchan = 67;
%     kolors(setdiff(1:64+length(METAchan),METAchan),:) = tmp2;
end

TRIGbitvals = uint16([2 64 4]);
TRIGbittext = {'READY', 'BATT LOW', 'TRIG'};
%

acqchan = [1:MaxChan];
Nchan = length(acqchan);
if strcmpi(systemMode,'fusbi')
    NchanPerRow = Nchan;
else
    NchanPerRow = floor(Nchan/2);
end
Ndatachan = Nchan - length(METAchan);
DataChans = setdiff(1:Nchan, METAchan);
RefreshPeriod = round(RefreshPeriod * SampleRate) / SampleRate;
AnalysisDuration = round(AnalysisDuration * SampleRate) / SampleRate;
Windowsize = AnalysisDuration * SampleRate;
Deltasize = RefreshPeriod * SampleRate;
FetchExtra = FetchExtraSec * SampleRate;
Qsize = SigQualLookbackTime * SampleRate;
Psize = min([TextPPLookbackTime SigQualLookbackTime]) * SampleRate;
WarnBufferFullBefore = WarnBufferFullBeforeSec * SampleRate;

% ADC parameters
for i = 1:67
    %ADCrange(i,:) = [-3 3];
    %ADCgain(i,:) = 20;
    %ADCbits(i,:) = 22;
    
    % 20130909: MindMedia says that SDK already converts to microvolts
    ADCrange(i,:) = [0 1];
    ADCgain(i,:) = 1e6;
    ADCbits(i,:) = 0;
end
for j = 1:length(STATUSchan)
    ADCrange(STATUSchan(j),:) = [0 1];
    ADCgain(STATUSchan(j),:) = 1;
    ADCbits(STATUSchan(j),:) = 0;
end
for j = 1:length(DIAGchan)
    ADCrange(DIAGchan(j),:) = [0 4];
    ADCgain(DIAGchan(j),:) = 63;
    ADCbits(DIAGchan(j),:) = 0;
end
ADCfactor = diff(ADCrange,[],2)./2.^ADCbits./ADCgain;
%

% Unit scaling parameters (applies to the text)
%GUNITfactor = 1e0;
%[~, NewPrefix] = log_scale_si_prefixes(1/UNITfactor);
%UNITname = [NewPrefix 'V'];
QTYname = '_{pp}';

% Display range scaling
for i = 1:67
    % What should the full swing (in volts) be?
    %DISPfull(i,:) = 2.5e-6;
    DISPfull(i,:) = 35e-6;
end
for j = 1:length(METAchan)
    DISPfull(METAchan(j),:) = 4;
end
DISPheight = 0.5; % The available height to display a full swing (in plot units)
DISPfullmedian = median(DISPfull(setdiff(1:size(DISPfull,1),METAchan),:));

globalscaleallow = AutoScaleAllow/DISPfullmedian;


try
    nexus_stop();
end

% Initialize device
%mp150_init_acq(SampleRate,acqchan);
try
    nexus_init(Nchan);
catch exception
    SWnohardwaremode = 1;
    FakeChanNames = {'3Hz', '7Hz', '11Hz', '17Hz', '23Hz', '31Hz', '75Hz', '90Hz', '100Hz', '110Hz', ...
        '3HzN1', '7HzN1', '11HzN1', '17HzN1', '23HzN1', '31HzN1', '75HzN1', '90HzN1', '100HzN1', '110HzN1', ...
        '3HzN2', '7HzN2', '11HzN2', '17HzN2', '23HzN2', '31HzN2', '75HzN2', '90HzN2', '100HzN2', '110HzN2', ...
        '3HzN3', '7HzN3', '11HzN3', '17HzN3', '23HzN3', '31HzN3', '75HzN3', '90HzN3', '100HzN3', '110HzN3', ...
        '3HzSquare', '7HzSquare', '11HzSquare', '17HzSquare', '23HzSquare', '31HzSquare', '75HzSquare', '90HzSquare', '100HzSquare', '110HzSquare'};
    OriginalChanNames(1:length(FakeChanNames)) = FakeChanNames;
    ChanNames = OriginalChanNames;
    disp('Error: Error initializing (the command "nexus_init" did not work).');
    disp(['The error message was: ' exception.message]);
    disp('The scope is continuing to load without a connection to the amplifiers.');
    %error('Error: Error initializing.');
    %return
end

try
    updatesystemstatus('nexus-init');
end

%Cursor = -Deltasize;

rawsignal = nan(Nchan,Windowsize);
rawvoltsignal = zeros(Nchan,Windowsize);
rawchunk = zeros(Nchan,Deltasize);
rawvoltchunk = rawchunk;
filteredvoltchunk = rawchunk;
Vhistory = zeros(Nchan,8);

missed = 0;
capaed = 0;
LastTrigTime = clock;
PrevTitleStr = '';

% close all
% close all hidden
fighand = figure(42974);
if SWmemcheck
    set(fighand, 'Visible', 'off');
else
    clf(fighand, 'reset');
end
set(fighand, 'NumberTitle', 'off', 'Name', 'nexus_scope');
set(fighand,'Menubar','none','Toolbar','none','DoubleBuffer','on');
%set(fighand,'Position',[1           31        1280         968]);
%set(fighand, 'Position', [49         264        1094         674]);
scrsz = get(0, 'ScreenSize');
scrsz = scrsz(1,:);
set(fighand, 'Position', scrsz + [0 157 -300 -190]);
set(fighand,'Color',BgdFigureColor);
figure(fighand);

%20200701
% https://www.mathworks.com/matlabcentral/answers/102219-how-do-i-make-a-figure-full-screen-programmatically-in-matlab
try
    if SWmaximize
        pause(0.00001);
        warning('off','MATLAB:ui:javacomponent:FunctionToBeRemoved');
        frame_h = get(handle(fighand),'JavaFrame');
        set(frame_h,'Maximized',1);
    end
end

thwarning = [];

WarnData.NotRecorded.String = WarningNotRecordedString;
WarnData.NotRecorded.FontSize = WarningNotRecordedFontSize;
WarnData.NotRecorded.Color1 = WarningNotRecordedColor1;
WarnData.NotRecorded.Color2 = WarningNotRecordedColor2;
WarnData.NotRecorded.BackgroundColor = 'none';
WarnData.NotRecorded.EdgeColor = 'none';
WarnData.NotRecorded.LineWidth = 0.5;
WarnData.NotRecorded.FontWeight = 'bold';
WarnData.NotRecorded.HorizontalAlignment = 'center';
WarnData.NotRecorded.VerticalAlignment = 'middle';
WarnData.NotRecorded.FontName = 'Arial Narrow';

WarnData.TriggerLost.String = WarningTriggerLostString;
WarnData.TriggerLost.FontSize = WarningNotRecordedFontSize;
WarnData.TriggerLost.Color1 = WarningNotRecordedColor1;
WarnData.TriggerLost.Color2 = WarningNotRecordedColor1;
WarnData.TriggerLost.BackgroundColor = 'none';
WarnData.TriggerLost.EdgeColor = 'none';
WarnData.TriggerLost.LineWidth = 0.5;
WarnData.TriggerLost.FontWeight = 'bold';
WarnData.TriggerLost.HorizontalAlignment = 'center';
WarnData.TriggerLost.VerticalAlignment = 'middle';
WarnData.TriggerLost.FontName = 'Arial Narrow';

WarnData.LowBattery.String = WarningBatteryLowString;
WarnData.LowBattery.FontSize = WarningNotRecordedFontSize;
WarnData.LowBattery.Color1 = WarningNotRecordedColor1;
WarnData.LowBattery.Color2 = WarningNotRecordedColor1;
WarnData.LowBattery.BackgroundColor = BattLowBackgroundColor;
WarnData.LowBattery.EdgeColor = BattLowEdgeColor;
WarnData.LowBattery.LineWidth = BattLowLineWidth;
WarnData.LowBattery.FontWeight = 'bold';
WarnData.LowBattery.HorizontalAlignment = 'center';
WarnData.LowBattery.VerticalAlignment = 'middle';
WarnData.LowBattery.FontName = 'Arial Narrow';

WarnData.DeviceError.String = WarningDeviceErrorString;
WarnData.DeviceError.FontSize = WarningNotRecordedFontSize;
WarnData.DeviceError.Color1 = WarningNotRecordedColor1;
WarnData.DeviceError.Color2 = WarningNotRecordedColor1;
WarnData.DeviceError.BackgroundColor = BattLowBackgroundColor;
WarnData.DeviceError.EdgeColor = BattLowEdgeColor;
WarnData.DeviceError.LineWidth = BattLowLineWidth;
WarnData.DeviceError.FontWeight = 'bold';
WarnData.DeviceError.HorizontalAlignment = 'center';
WarnData.DeviceError.VerticalAlignment = 'middle';
WarnData.DeviceError.FontName = 'Arial Narrow';

WarnData.NoHardware.String = WarningNoHardwareString;
WarnData.NoHardware.FontSize = WarningNotRecordedFontSize;
WarnData.NoHardware.Color1 = WarningNotRecordedColor1;
WarnData.NoHardware.Color2 = WarningNotRecordedColor1;
WarnData.NoHardware.BackgroundColor = BattLowBackgroundColor;
WarnData.NoHardware.EdgeColor = BattLowEdgeColor;
WarnData.NoHardware.LineWidth = BattLowLineWidth;
WarnData.NoHardware.FontWeight = 'bold';
WarnData.NoHardware.HorizontalAlignment = 'center';
WarnData.NoHardware.VerticalAlignment = 'middle';
WarnData.NoHardware.FontName = 'Arial Narrow';


WarnData.BufferFull.String = WarningBufferFullString;
WarnData.BufferFull.FontSize = WarningNotRecordedFontSize;
WarnData.BufferFull.Color1 = WarningNotRecordedColor1;
WarnData.BufferFull.Color2 = WarningNotRecordedColor1;
WarnData.BufferFull.BackgroundColor = BattLowBackgroundColor;
WarnData.BufferFull.EdgeColor = BattLowEdgeColor;
WarnData.BufferFull.LineWidth = BattLowLineWidth;
WarnData.BufferFull.FontWeight = 'bold';
WarnData.BufferFull.HorizontalAlignment = 'center';
WarnData.BufferFull.VerticalAlignment = 'middle';
WarnData.BufferFull.FontName = 'Arial Narrow';

% thwarningnotrecorded = [];
% thwarningtriglost = [];
% thwarninglowbattery = [];
% thwarningbufferfull = [];
% thwarningdeviceerror = [];
%pos2 = -MaxChanInRow/2;
%pos3 = -MaxChanInRow/2;
warningpos = -NchanPerRow/2;
pos3 = -NchanPerRow/2;
if Nchan <= NchanPerRow
    %set(gca,'Position',[0.04 0.03 0.86 0.94]);
    subplothand(1) = gca;
    set(gca,'Position',[0.04 0.03 0.86 0.94]);
    hold on
else
    subplothand(1) = subplot(1,2,1);
    set(gca,'Position',[0.04 0.03 0.36 0.94]);
    hold on
    subplothand(2) = subplot(1,2,2);
    set(gca,'Position',[0.54 0.03 0.36 0.94]);
    hold on
end

Choffsets = zeros(1,Nchan);

for ch = 1:Nchan
    if Nchan <= NchanPerRow
        choffset = ch;
    else
        if ch <= NchanPerRow
            subplot(subplothand(1));
            choffset = ch;
        else
            subplot(subplothand(2))
            choffset = ch - NchanPerRow;
        end
    end
    plot([1:PlotEveryNSample:size(rawsignal,2)]/SampleRate,(rawsignal(ch,1:PlotEveryNSample:end)-median(rawsignal(ch,1:PlotEveryNSample:end)))./(max(rawsignal(ch,1:PlotEveryNSample:end))-min(rawsignal(ch,1:PlotEveryNSample:end))+eps)/2-choffset,'Color',kolors(ch,:));
    tmp = get(gca,'Children');
    Kids(ch) = tmp(1);
    Choffsets(ch) = choffset;
    Grandkids(ch) = text(size(rawsignal,2)/SampleRate,-choffset,'INFODISPLY','BackgroundColor',BgdFontColor,'Color',1-BgdFontColor,'FontName','Consolas','FontSize',FontSize,'FontWeight','bold','LineWidth',0.50, 'BackgroundColor','none');
    %set(gca,'YLim',[-Nchan-1,0],'YTick',[-Nchan:-1],'YTickLabel',[Nchan:-1:1]);
    
    if ismember(ch, DIAGchan)
        set(Kids(ch), 'Visible', 'off');
        set(Grandkids(ch), 'Position', [size(rawsignal,2)/SampleRate,0.5,0]);
    end
    
end

if ~TwoColumnMode
    thwarning = text(mean([0 AnalysisDuration]), warningpos, '');
    %thwarningnotrecorded = text(mean([0 AnalysisDuration]), pos3, WarningNotRecordedString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow');
    %thwarningtriglost = text(mean([0 AnalysisDuration]), pos3, WarningTriggerLostString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow');
    %thwarninglowbattery = text(mean([0 AnalysisDuration]), warningpos, WarningBatteryLowString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    %thwarningbufferfull = text(mean([0 AnalysisDuration]), warningpos, WarningBufferFullString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    %thwarningdeviceerror = text(mean([0 AnalysisDuration]), warningpos, WarningDeviceErrorString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
else
    subplot(subplothand(1));
    thwarning(1) = text(mean([0 AnalysisDuration]), warningpos, '');
    %thwarningnotrecorded(1) = text(mean([0 AnalysisDuration]), pos3, WarningNotRecordedString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow');
    %thwarningtriglost(1) = text(mean([0 AnalysisDuration]), pos3, WarningTriggerLostString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow');
    %thwarninglowbattery(1) = text(mean([0 AnalysisDuration]), warningpos, WarningBatteryLowString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    %thwarningbufferfull(1) = text(mean([0 AnalysisDuration]), warningpos, WarningBufferFullString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    %thwarningdeviceerror(1) = text(mean([0 AnalysisDuration]), warningpos, WarningDeviceErrorString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    subplot(subplothand(2));
    thwarning(2) = text(mean([0 AnalysisDuration]), warningpos, '');
    %thwarningnotrecorded(2) = text(mean([0 AnalysisDuration]), pos3, WarningNotRecordedString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow');
    %thwarningtriglost(2) = text(mean([0 AnalysisDuration]), pos3, WarningTriggerLostString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow');
    %thwarninglowbattery(2) = text(mean([0 AnalysisDuration]), warningpos, WarningBatteryLowString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    %thwarningbufferfull(2) = text(mean([0 AnalysisDuration]), warningpos, WarningBufferFullString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
    %thwarningdeviceerror(2) = text(mean([0 AnalysisDuration]), warningpos, WarningDeviceErrorString, 'FontSize', WarningNotRecordedFontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', WarningNotRecordedColor1, 'FontName', 'Arial Narrow', 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
end

set(thwarning, 'Visible', 'off');

% set(thwarningnotrecorded, 'Visible', 'off');
% set(thwarningtriglost, 'Visible', 'off');
% set(thwarninglowbattery, 'Visible', 'off');
% set(thwarningbufferfull, 'Visible', 'off');
% set(thwarningdeviceerror, 'Visible', 'off');


if ~TwoColumnMode
    set(gca,'YLim',[-Nchan-1,0],'YTick',[-Nchan:-1],'YTickLabel',[Nchan:-1:1]);
else
    for i = 1:2
        set(subplothand(i),'YLim',[-NchanPerRow-1,0],'YTick',[-NchanPerRow:-1],'YTickLabel',(i-1)*NchanPerRow+[NchanPerRow:-1:1]);
    end
end

[ChanNames, Chans] = populate_channames (ChanNames, Nchan, acqchan, DIAGchan, NchanPerRow, TwoColumnMode, subplothand);

% if ~isempty(who('ChanNames')) && ~isempty(ChanNames) && iscell(ChanNames)
%     if length(ChanNames) < Nchan
%         tmp = length(ChanNames);
%         tmp = tmp + 1;
%         while tmp <= Nchan
%             ChanNames{tmp} = num2str(acqchan(tmp));
%             tmp = tmp + 1;
%         end
%     end
%     Chans = ChanNames(1:Nchan);
%     Chans(DIAGchan) = {''};
%     for i = 1:length(Chans)
%         YLabel{i} = Chans{i};
%     end
%     if ~TwoColumnMode
%         set(gca,'YTickLabel',fliplr(YLabel));
%     else
%         for i = 1:2
%             set(subplothand(i),'YTickLabel',fliplr(YLabel((i-1)*NchanPerRow+1:(i)*NchanPerRow)));
%         end
%     end
% end

if ~TwoColumnMode
    set(gca,'XLim',[0 AnalysisDuration]);
    %tmp = get(gca,'XTick');
    set(gca, 'XTick', [0 1]);
    %set(gca,'XTickLabel',fliplr(tmp));
    PlotYLim = get(gca,'YLim');
    text(0.5,PlotYLim(1)-0.4,'1 s','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
    EventIndicatorLight = text(1.5,PlotYLim(1)-0.4,'EVNT','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
    EventIndicatorText = text(2.5,PlotYLim(1)-0.4,'EVNT','HorizontalAlignment','left','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
    for t = 0:AnalysisDuration
        TimeKid(1,t+1) = text(t,PlotYLim(1)-0.4,'|','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
    end
    titlehand = title('EEG Scope','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
    %Kids = get(gca,'Children');
    axis off
else
    for i = 1:2
        set(subplothand(i),'XLim',[0 AnalysisDuration]);
        %tmp = get(subplothand(i),'XTick');
        %set(subplothand(i),'XTickLabel',-fliplr(tmp));
        set(gca, 'XTick', [0 1]);
        subplot(subplothand(i));
        PlotYLim = get(gca,'YLim');
        text(0.5,PlotYLim(1)-0.4,'1 s','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
        EventIndicatorLight(i) = text(1.5,PlotYLim(1)-0.4,'EVNT','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
        EventIndicatorText(i) = text(2.5,PlotYLim(1)-0.4,'EVNT','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
        for t = 0:AnalysisDuration
            TimeKid(i,t+1) = text(t,PlotYLim(1)-0.4,'|','HorizontalAlignment','center','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
        end
        titlehand(i) = title('EEG Scope','Color',1-BgdFontColor,'BackgroundColor',BgdFontColor);
        %Kids = get(gca,'Children');
        axis off
    end
end

scopecontrol.startrecording = uicontrol(fighand, 'Style','checkbox', 'String','REC','FontSize',12,'FontWeight','bold');
scopecontrol.startrecording_touched = clock;
set(scopecontrol.startrecording, 'ForegroundColor', [1 0.4 0.4], 'BackgroundColor', [0 0 0]);
set(scopecontrol.startrecording, 'Unit', 'normalized', 'Position', [0.93, 0.02, 0.10, 0.015]);

scopecontrol.paneltoggle = uicontrol(fighand, 'Style','togglebutton', 'String','CP','FontSize',8,'FontWeight','bold');
set(scopecontrol.paneltoggle, 'ForegroundColor', [0.6 0.6 0.6], 'BackgroundColor', [0 0 0]);
set(scopecontrol.paneltoggle, 'Unit', 'normalized', 'Position', [0.01, 0.02, 0.02, 0.02]);

scopecontrol.panel = uipanel(fighand, 'Title', 'Control Panel', 'FontWeight', 'bold', 'FontSize', 12);
set(scopecontrol.panel, 'ForegroundColor', [0.9 0.9 0.9], 'BackgroundColor', [0 0 0]);
%set(scopecontrol.panel, 'Unit', 'normalized', 'Position', [0.05, 0.05, 0.20, 0.40]);
set(scopecontrol.panel, 'Unit', 'pixels', 'Position', [65, 49.4, 256, 387.2]);
 

scopecontrol.sampleratetext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Sample rate:');
scopecontrol.samplerate = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', num2str(SampleRate));

scopecontrol.refreshratetext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Refresh rate:');
scopecontrol.refreshrate = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', num2str(1/RefreshPeriod));

scopecontrol.segmentbufferusedtext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Seg. buffer used:');
scopecontrol.segmentbufferused = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', '0');

scopecontrol.segmentbufferavailabletext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Seg. buffer avail.:');
scopecontrol.segmentbufferavailable = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', '0');

scopecontrol.contigbufferusedtext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Contig buffer used:');
scopecontrol.contigbufferused = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', '0');

scopecontrol.contigbufferavailabletext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Contig buffer avail.:');
scopecontrol.contigbufferavailable = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', '0');



scopecontrol.noplot = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Disable Plotting');
scopecontrol.car = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Common Average Reference');
scopecontrol.numericchannames = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Num Chan', 'Callback', @scopecontrol_numeric_channames);

scopecontrol.filter = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Bandpass');
scopecontrol.filterhighpasstext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'HPF:');
scopecontrol.filterhighpass = uicontrol(scopecontrol.panel, 'Style', 'edit', 'String', num2str(Fcutoffhigh));
scopecontrol.filterlowpasstext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'LPF:');
scopecontrol.filterlowpass = uicontrol(scopecontrol.panel, 'Style', 'edit', 'String', num2str(Fcutoff));
scopecontrol.filtertext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Filters only affect displayed traces.');
scopecontrol.notchfilter = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Notch');

scopecontrol.envelope = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Envelope');
scopecontrol.envelopelowpasstext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'ENV:');
scopecontrol.envelopelowpass = uicontrol(scopecontrol.panel, 'Style', 'edit', 'String', num2str(Fenvelope));

scopecontrol.autoscale = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Autoscale');
scopecontrol.zscore = uicontrol(scopecontrol.panel, 'Style', 'checkbox', 'String', 'Zscore');
scopecontrol.zoomtext = uicontrol(scopecontrol.panel, 'Style', 'text', 'String', 'Vertical Zoom:');
scopecontrol.zoom = uicontrol(scopecontrol.panel, 'Style', 'slider', 'String', 'Thingy');
set(scopecontrol.zoom, 'Min', log2(globalscaleallow(1)), 'Max', log2(globalscaleallow(end)), 'Value', log2(globalscaleallow(1)));
scopecontrol.reloadchannames = uicontrol(scopecontrol.panel, 'Style', 'pushbutton', 'String', 'Update ChanNames', 'Callback', @scopecontrol_reload_channames);
scopecontrol.restartnexus = uicontrol(scopecontrol.panel, 'Style', 'pushbutton', 'String', RestartNexusButton.normaltext, 'Callback', @scopecontrol_restart_nexus);

scopecontrol.zoom_in = uicontrol(scopecontrol.panel, 'Style', 'pushbutton', 'String', '+', 'Callback', @scopecontrol_zoom_in);
scopecontrol.zoom_out = uicontrol(scopecontrol.panel, 'Style', 'pushbutton', 'String', '-', 'Callback', @scopecontrol_zoom_out);
scopecontrol.pan_up = uicontrol(scopecontrol.panel, 'Style', 'pushbutton', 'String', '^', 'Callback', @scopecontrol_pan_up);
scopecontrol.pan_down = uicontrol(scopecontrol.panel, 'Style', 'pushbutton', 'String', 'v', 'Callback', @scopecontrol_pan_down);

clear tmp
tmp.gca = gca;
tmp.ylim = get(gca, 'YLim');
setappdata(scopecontrol.zoom_in, 'data', tmp);
setappdata(scopecontrol.zoom_out, 'data', tmp);
setappdata(scopecontrol.pan_up, 'data', tmp);
setappdata(scopecontrol.pan_down, 'data', tmp);
clear tmp

setappdata(scopecontrol.reloadchannames, 'ToReloadChanNames', 0);
setappdata(scopecontrol.restartnexus, 'ToRestartNexus', 0);
scopecontrol.restartnexus_touched = clock;
if globalscale >= globalscaleallow(1) && globalscale <= globalscaleallow(2)
    set(scopecontrol.zoom, 'Value', log2(globalscale));
end

panelitems = {
    scopecontrol.pan_up 1 0.9 0.07
    scopecontrol.zoom_in 2 0.9 0.07
    scopecontrol.zoom_out 3 0.9 0.07
    scopecontrol.pan_down 4 0.9 0.07
    scopecontrol.sampleratetext 1 0 0.40
    scopecontrol.samplerate 1 0.40 0.30
    scopecontrol.refreshratetext 2 0.0 0.40
    scopecontrol.refreshrate 2 0.40 0.30
    scopecontrol.segmentbufferusedtext 3 0 0.55
    scopecontrol.segmentbufferused 3 0.55 0.30
    scopecontrol.segmentbufferavailabletext 4 0 0.55
    scopecontrol.segmentbufferavailable 4 0.55 0.30
    scopecontrol.contigbufferusedtext 5 0 0.55
    scopecontrol.contigbufferused 5 0.55 0.30
    scopecontrol.contigbufferavailabletext 6 0 0.55
    scopecontrol.contigbufferavailable 6 0.55 0.30
    scopecontrol.noplot 7.5 0 0.45
    scopecontrol.numericchannames 7.5 0.45 0.45
    scopecontrol.car 9 0 0.9
    scopecontrol.filter 10.5 0 0.33
    scopecontrol.notchfilter 10.5 0.35 0.27
    scopecontrol.filterhighpasstext 11.5 0 0.13
    scopecontrol.filterhighpass 11.5 0.13 0.15
    scopecontrol.filterlowpasstext 11.5 0.3 0.13
    scopecontrol.filterlowpass 11.5 0.43 0.15
    scopecontrol.filtertext 12.5 0 0.9

    scopecontrol.envelope 10.5 0.65 0.33
    scopecontrol.envelopelowpasstext 11.5 0.65 0.13
    scopecontrol.envelopelowpass 11.5 0.78 0.15
    
    scopecontrol.autoscale 14 0 0.4
    scopecontrol.zscore 14 0.5 0.4
    scopecontrol.zoomtext 15 0 0.9
    scopecontrol.zoom 16 0 0.9
    scopecontrol.reloadchannames 17.5 0 0.51
    scopecontrol.restartnexus 17.5 0.56 0.37
    
    };

for i = 1:size(panelitems,1)
    if ishandle(panelitems{i,1})
        set(panelitems{i,1}, 'FontSize', 10, 'ForegroundColor', [0.9 0.9 0.9], 'Backgroundcolor', [0 0 0], 'HorizontalAlignment', 'left');
        set(panelitems{i,1}, 'Unit', 'normalized', 'Position', [0.02+panelitems{i,3}, 0.97-0.05*panelitems{i,2}, panelitems{i,4}, 0.05]);
    end
end

set(scopecontrol.notchfilter, 'Visible', 'off');

%keyboard

if ~SWplot
    set(scopecontrol.noplot, 'Value', 1);
end
if SWcar
    set(scopecontrol.car, 'Value', 1);
end
if SWfilter
    set(scopecontrol.filter, 'Value', 1);
end
if SWnotch
    set(scopecontrol.notchfilter, 'Value', 1);
end
if SWautoscale
    set(scopecontrol.autoscale, 'Value', 1);
end
if SWzscore
    set(scopecontrol.zscore, 'Value', 1);
end

filteredvoltsignal = zeros(Nchan,Windowsize);
%centeredsignal = zeros(Nchan,Windowsize);
rawshortsignal = nan(Nchan,Qsize);
%peaktopeak = nan(Nchan,1);
%fpeaktopeakvoltages = nan(Nchan,1);
%longpeaktopeak = nan(Nchan,1);
%peaktopeakvoltages = nan(Nchan,1);
FreeTime = zeros(1,Windowsize);
ows = ones(1,Windowsize);
%ocs = ones(1,Deltasize);

% Start acquisition (META channels are not filtered)
%mp150_start_acq_continuous(); % begin acquiring (push red button)

%witlim = 64;
%20100721
%if SWfilter
Zf_initial = zeros(length(setdiff(1:Nchan,METAchan)),nlookback).';
Zf = Zf_initial;
Zf_notch_initial = zeros(length(setdiff(1:Nchan,METAchan)),max(length(notchHd{1}.sosMatrix(1,:)), length(notchHd{1}.ScaleValues))-1).';
notchZf = Zf_notch_initial;
%end
Zf_env_initial = zeros(length(setdiff(1:Nchan,METAchan)),nlookback_env).';
Zf_env = Zf_env_initial;

%suid = datestr(now,'yyyymmdd-HHMMSS-FFF');
suid = dec2base(floor(unixtime*1000),36);
SavePathName = gettmpdir();
if ~exist('SaveDir','var') || isempty(SaveDir) || ~exist(SaveDir,'dir')
    SaveDir = SavePathName;
end

SaveFileName = ['nexusautosave-' suid '.mat'];
SaveUNC = [SavePathName filesep SaveFileName];
mmfilename = [SavePathName filesep 'nexus_scope_mmf_0.dat'];
set(thwarning, 'String', 'Preparing..', 'Visible', 'on');
drawnow
retry = 1;
while retry
    if retry >= 100
        error('Critical: Cannot open any memory mapped files for writing. Check that the temp directory is writable. Suggest restarting computer.');
    end
    [fid, errmsg] = fopen(mmfilename, 'w');
    if isempty(errmsg)
        retry = 0;
        break
    end
    warning('Cannot open %s for writing. Using the next file instead.\n', mmfilename);
    mmfilename = [SavePathName filesep 'nexus_scope_mmf_' num2str(retry) '.dat'];
    retry = retry + 1;
end
estimatedsize = 8*4 + 8*(Nchan+1)*BufferDurationSec*SampleRate + 2*1*BufferDurationSec*SampleRate;
fprintf('Creating a temp file for memory mapped file (%g MiB) to buffer %g seconds at %g Hz sample rate...\n', estimatedsize/1048576, BufferDurationSec, SampleRate);
fwrite(fid, zeros(1,2), 'double');
fwrite(fid, zeros(1,1), 'double');
fwrite(fid, zeros(1,1), 'double');
fwrite(fid, zeros((Nchan+1)*BufferDurationSec*SampleRate,1), 'double');
fwrite(fid, zeros(1*BufferDurationSec*SampleRate,1), 'uint16');
fclose(fid);
mmfobj = memmapfile(mmfilename, 'Format', ...
    {
    'double', [1,2], 'dimensions'
    'double', [1,1], 'samplerate'
    'double', [1,1], 'recordendposition'
    'double', [Nchan+1, BufferDurationSec*SampleRate], 'rawdata'
    'uint16', [1, BufferDurationSec*SampleRate], 'rawevent'
    }, ...
    'Repeat', 1, 'Writable', true); % The last channel is for time stamp
buffersize = size(mmfobj.Data.rawdata,2);
mmfobj.Data.dimensions = [Nchan+1, BufferDurationSec*SampleRate];
mmfobj.Data.samplerate = SampleRate;
rawdata = [];
rawevent = [];
fprintf('Memory mapped file: %s \n',mmfilename);

if nargout >= 2 && ~SWnohardwaremode
    SWrecord = 1;
    RecordStart = 1;
    segmentsuid = datestr(now,'yyyymmdd-HHMMSS-FFF');
    SegmentSaveFileName = ['nexusautosave-' suid '_segment-' segmentsuid '.mat'];
    SegmentSaveUNC = [SaveDir filesep SegmentSaveFileName];
else
    SWrecord = 0;
    Warned.NotRecorded = 1;
    RecordStart = 1;
end
set(thwarning, 'Visible', 'off');
drawnow


if ~SWnohardwaremode
    try
        nexus_start(SampleRate, AcqDuration); %20121126 changed from 2.0 sec
    catch exception
        SWnohardwaremode = 1;
        disp('Error: Error starting data acquisition (the command "nexus_start" did not work).');
        disp(['The error message from the command was: ' exception.message]);
        disp('The scope is continuing to load without a connection to the amplifiers.');
        %disp('     (1) Check if the amplifiers are turned on');
        %disp('     (2) Change batteries if they are running low (if you use batteries)');
        %disp('     (3) Check if SynFi or FUSBI is plugged in on both ends (USB and fiber)');
        %disp('     (4) Save your work if necessary, restart MATLAB, and try again');
        %if ishandle(fighand)
        %    close(fighand);
        %end
        rawdata(1) = NaN;
        %return
    end
end
CursorPosition = 0;
RecordedSamples = 0;
BattLowWarn = zeros(1,Nchan);
%DeviceErrorWarn = zeros(1,Nchan);
FailedToStart = zeros(1,Nchan);
WarningCount.DeviceError = 0;
DeviceHasBeenReady = 0;
errorfraction = 0;
HasChanges = zeros(1,Nchan);
ChannelStd = zeros(1,Nchan);
ChSpres = zeros(1,Nchan);


set(gcf, 'KeyPressFcn', @f_fig_keypress);

PrevTitleStr = '';


gcfpos = get(gcf, 'Position');
prevgcfpos = gcfpos;



% main loop
Timer.loopstart = clock;
Timer.recordstart = clock;
Timer.lastactivity = clock;
for i = 1:intmax
    tic
    %rawchunk = mp150_fetch_acq_different(SampleRate,RefreshPeriod,acqchan); % get and stop acquiring
    Timer.acqstart = clock;
    %20121025: Acquire more than requested amount, then trim
    %rawchunk = nexus_getdata(SampleRate,RefreshPeriod,Nchan);
    if ~SWnohardwaremode
        rawsafechunk = nexus_getdata(SampleRate,RefreshPeriod+FetchExtraSec,Nchan,[],1);
    else
        rawsafechunk = fakenexus_getsimulateddata(SampleRate,RefreshPeriod+FetchExtraSec,Nchan,[],1);
    end
    tmp = find(rawsafechunk(STATUSchan(1),:),1,'last');
    if isempty(tmp)
        tmp = 0;
    end
    gapsize = size(rawsafechunk,2) - tmp;
    if gapsize <= 0 && FetchExtraSec <= RefreshPeriod / 4
        FetchExtraSec = FetchExtraSec * 2;
        fprintf('[%10.0f ms] Increased FetchExtraSec to %g sec.\n',etime(clock,Timer.loopstart)*1000, FetchExtraSec);
    end
    rawsafechunk = rawsafechunk(:,rawsafechunk(STATUSchan(1),:)>0);
    
    Deltasize = size(rawsafechunk,2);
    rawchunk = rawsafechunk;
    
%     tmp = size(rawsafechunk,2);
%     if tmp >= Deltasize
%         if tmp > Deltasize
%             Deltasize = tmp;
%             %ocs = ones(1,Deltasize);
%             %fprintf('warning: buffer overflow %i >= %i. excess data points are lost.\n', tmp, Deltasize);
%         end
%         rawchunk = rawsafechunk(:, 1:Deltasize);
%     else
%         rawchunk(:) = 0;
%         rawchunk(:,1:tmp) = rawsafechunk;
%     end
    %rawchunk = rawsafechunk(:, end-Deltasize+1:end);
    
    Timer.acqtotal = etime(clock,Timer.acqstart);
    Timer.procstart = clock;

    % immediately apply voltage conversion here:
    rawvoltchunk = rawchunk .* (ADCfactor(1:Nchan,:)*ones(1,Deltasize));
    
    ValidChunkIdx = rawchunk(STATUSchan(1),:) > 0;
    
    %20120308 : Always filter.
    % % Wait until device is ready before doing filtering
    %if min(min(rawchunk(METAchan,:))) == 0 || max(peaktopeakvoltages) > MaxSaneVoltage
    %    filteredvoltchunk = rawvoltchunk;
    %    %20100721 %Zf = zeros(length(setdiff(1:Nchan,METAchan)),nlookback).';
    %else
    
    filteredvoltchunk = rawvoltchunk;
    
    % notch 20210609 - does not work yet
    if SWnotch
        %for h = 1:floor(SampleRate/2/Funda)
            %[filteredvoltchunk(setdiff(1:Nchan,METAchan),:), notchZf] = filter(notchHd{h}.sosMatrix(1,:), notchHd{h}.ScaleValues, filteredvoltchunk(setdiff(1:Nchan,METAchan),:), notchZf, 2);
        %end
    end
    
    
    
    % car / filter
    
    %20100805,20120308
    if SWcar
        filteredvoltchunk(setdiff(1:Nchan,METAchan),:) = filteredvoltchunk(setdiff(1:Nchan,METAchan),:) - ones(Ndatachan,1) * mean(filteredvoltchunk(setdiff(1:Nchan,METAchan),:),1);
    end
    
    
    %20100721,20120308
    if SWfilter && ~isempty(filteredvoltchunk)
        [filteredvoltchunk(setdiff(1:Nchan,METAchan),:),Zf] = filter(FiltcoefB,FiltcoefA,filteredvoltchunk(setdiff(1:Nchan,METAchan),:),Zf,2);
    else
        %filteredvoltchunk(setdiff(1:Nchan,METAchan),:) = filteredvoltchunk(setdiff(1:Nchan,METAchan),:);
    end
    %end
    
    %20130930
    if SWenvelope && ~isempty(filteredvoltchunk)
        filteredvoltchunk(setdiff(1:Nchan,METAchan),:) = filteredvoltchunk(setdiff(1:Nchan,METAchan),:).^2;
        [filteredvoltchunk(setdiff(1:Nchan,METAchan),:),Zf_env] = filter(EnvelopeFiltcoefB,EnvelopeFiltcoefA,filteredvoltchunk(setdiff(1:Nchan,METAchan),:),Zf_env,2);
        %filteredvoltchunk = filteredvoltchunk.^0.5;
    end
    
    rawvoltsignal = rawvoltsignal(:,[Deltasize+1:end,1:Deltasize]);
    filteredvoltsignal = filteredvoltsignal(:,[Deltasize+1:end,1:Deltasize]);
    rawshortsignal = rawshortsignal(:,[Deltasize+1:end,1:Deltasize]);
    
    if SWrecord
        Timer.lastactivity = clock;
        if RecordedSamples+Deltasize > size(mmfobj.Data.rawdata,2)-WarnBufferFullBefore
            if ~Warned.BufferFull
                %set(thwarningbufferfull, 'Visible', 'on');
                Warned.BufferFull = 1;
            end
        end
        
        if RecordedSamples+Deltasize > size(mmfobj.Data.rawdata,2)
            T = 'STOPPED';
            CC = FailedFontColor;
            while length(T) < 10
                T = [' ' T];
            end
            for ch = 1:Nchan
                T2 = Chans{ch};
                while length(T2) < 4
                    T2 = [' ' T2];
                end
                set(Grandkids(ch),'String',[T2 ' ' T],'Color',CC);
            end
            fprintf('Recording stopped. Data exceeded rawdata buffer.\n');
            break
        end
        if etime(clock,LastTrigTime) > TrigLostStopWaitWarn && etime(clock,Timer.recordstart) > 60+TrigLostStopWaitWarn-TrigLostStopWait
            if ~Warned.TriggerLost
                %set(thwarningtriglost, 'Visible', 'on');
                Warned.TriggerLost = 1;
            end
            if etime(clock,LastTrigTime) > TrigLostStopWait && etime(clock,Timer.recordstart) > 60
                T = 'STOPPED';
                CC = FailedFontColor;
                while length(T) < 10
                    T = [' ' T];
                end
                for ch = 1:Nchan
                    T2 = Chans{ch};
                    while length(T2) < 4
                        T2 = [' ' T2];
                    end
                    set(Grandkids(ch),'String',[T2 ' ' T],'Color',CC);
                end
                fprintf('Recording stopped. Trigger signal lost for %g sec.\n',TrigLostStopWait);
                break
            end
        elseif Warned.TriggerLost
            %set(thwarningtriglost, 'Visible', 'off');
            Warned.TriggerLost = 0;
        end
    else
        Warned.TriggerLost = 0;
    end
    
    a = CursorPosition+1;
    b = CursorPosition+Deltasize;
    if ~isempty(a:b)
        ma = mod(a-1, buffersize)+1;
        mb = mod(b-1, buffersize)+1;
        if mb >= ma
            mmfobj.Data.rawdata(1:end-1,ma:mb) = rawchunk;
            mmfobj.Data.rawdata(end,ma:mb) = a:b;
        else
            % Wrap around required
            mmfobj.Data.rawdata(1:end-1,[ma:buffersize,1:mb]) = rawchunk;
            mmfobj.Data.rawdata(end,[ma:buffersize,1:mb]) = a:b;
        end
    end

    rawvoltsignal(:,[end-Deltasize+1:end]) = rawvoltchunk;
    filteredvoltsignal(:,[end-Deltasize+1:end]) = filteredvoltchunk;
    %filtereddata(:,CursorPosition+1:CursorPosition+Deltasize) = filteredvoltchunk;
    rawshortsignal(:,[end-Deltasize+1:end]) = rawchunk;
    peaktopeak = max(rawshortsignal(:,end-Psize+1:end),[],2) - min(rawshortsignal(:,end-Psize+1:end),[],2);
    fpeaktopeakvoltages = max(filteredvoltsignal(:,end-Psize+1:end),[],2) - min(filteredvoltsignal(:,end-Psize+1:end),[],2);
    longpeaktopeak = max(rawvoltsignal,[],2) - min(rawvoltsignal,[],2);
    
    if mod(i-1,0.5/RefreshPeriod)==0
        Vhistory = [Vhistory(:,2:end) peaktopeak(:,1)];
    end
    CursorPosition = CursorPosition + Deltasize;
    mmfobj.Data.recordendposition = CursorPosition;
    if SWrecord
        RecordedSamples = RecordedSamples + Deltasize;
    end
    
    if SWzscore
        centeredsignal = zscore(filteredvoltsignal,[],2);
    %elseif any(CSWmediandisp)
        %centeredsignal = filteredvoltsignal - median(filteredvoltsignal,2) * ows;
    else
        %centeredsignal = filteredvoltsignal;
        centeredsignal = filteredvoltsignal - median(filteredvoltsignal,2) * ows;
    end
    
    peaktopeakvoltages = peaktopeak .* ADCfactor(1:Nchan,:);
    
    if SWautoscale
        if SWzscore
            newglobalscale = min(4/DISPfullmedian, globalscaleallow(2));
        else
            tmp = fpeaktopeakvoltages(intersect(DataChans,find(longpeaktopeak~=0)));
            newglobalscale = 4^round( log2(median(tmp)/DISPfullmedian)/2 );
            if isnan(newglobalscale)
                newglobalscale = 0;
            end
        end
        if newglobalscale < globalscaleallow(1)
            newglobalscale = globalscaleallow(1);
        elseif newglobalscale > globalscaleallow(2)
            newglobalscale = globalscaleallow(2);
        end
        if SWzscore
            globalscale = newglobalscale;
        else
            if newglobalscale > globalscale
                globalscale = min(newglobalscale, globalscale * 2);
            elseif newglobalscale < globalscale
                globalscale = max(newglobalscale, globalscale / 2);
            end
        end
        set(scopecontrol.zoom, 'Value', log2(globalscale));
    else
        newglobalscale = 2^get(scopecontrol.zoom, 'Value');
        if isnan(newglobalscale)
            newglobalscale = 0;
        end
        if newglobalscale < globalscaleallow(1)
            newglobalscale = globalscaleallow(1);
        elseif newglobalscale > globalscaleallow(2)
            newglobalscale = globalscaleallow(2);
        end
        globalscale = newglobalscale;
    end
    
    NumNormalChannels = 0;
    
    for ch = 1:Nchan
        %if Nchan <= NchanPerRow || ch <= NchanPerRow
        %    choffset = ch;
        %else
        %    choffset = ch - NchanPerRow;
        %end
        choffset = Choffsets(ch);
        
        if SWplot
            sigpres = 1;
            if longpeaktopeak(ch) == 0
                sigpres = NaN;
            end
            S = centeredsignal(ch,1:PlotEveryNSample:end) * sigpres;
            if intersect(METAchan,ch)
                S = S/DISPfull(ch,:)*DISPheight;
            else
                S = S/DISPfull(ch,:)/globalscale*DISPheight;
            end
            
            % Cap the ceiling and floor
            S(S>DISPheight/2*DispSpillAllow) = DISPheight/2*DispSpillAllow;
            S(S<-DISPheight/2*DispSpillAllow) = -DISPheight/2*DispSpillAllow;
            
            S = S  - Choffsets(ch);
            
            if SWwiper
                WiperPos = mod(CursorPosition,Windowsize)/Windowsize*size(S,2);
                WiperWidth = max(1,floor(AnalysisDuration/100*SampleRate/PlotEveryNSample));
                % This makes the oldest samples unplotted.
                S(1,1:WiperWidth) = NaN;
                % This makes the plot appear non-moving.
                %20121026: Round.
                S = circshift(S, [0 round(WiperPos)]);
            end
            set(Kids(ch),'YData',S);
            ChSpres(ch) = sigpres;
        end
        
        
        [UNITfactor, UNITname] = get_best_unitfactor(mean(peaktopeakvoltages(ch,1)));
        V = peaktopeakvoltages(ch,1) * UNITfactor;
        
        if intersect(STATUSchan,ch)
            %bittest = bitand(uint16(max(filteredvoltchunk(ch,:))),TRIGbitvals)~=0;
            tmp = max(rawchunk(ch,:));
            if isempty(tmp)
                tmp = 0;
            end
            bittest = bitand(uint16(tmp),TRIGbitvals)~=0;
            if bittest(strcmpi(TRIGbittext,'BATT LOW'))
                BattLowWarn(ch) = 1;
            else
                BattLowWarn(ch) = 0;
            end
            if bittest(strcmpi(TRIGbittext,'READY'))
                DeviceHasBeenReady = 1;
                FailedToStart(ch) = 0;
            elseif DeviceHasBeenReady || etime(clock,Timer.loopstart) > 20
                FailedToStart(ch) = 1;
            end
        end
        
        if (isempty(rawchunk) || min(min(rawchunk(STATUSchan,:))) == 0) && ~DeviceHasBeenReady
            T = 'Starting';
            CC = WaitFontColor;
            LastTrigTime = clock;
        %elseif any(isnan(rawvoltsignal(ch,:)))
        %    T = 'Buffering';
        %    CC = WaitFontColor;
        %    LastTrigTime = clock;
        elseif abs(V) > MaxSaneVoltage*UNITfactor && SWovervoltagewarning && isempty(intersect(METAchan,ch))
            T = 'Over Voltage';
            CC = AngryFontColor;
            set(Kids(ch),'Color',AngryLineColor);
            NumNormalChannels = NumNormalChannels - Ndatachan;
        else
            
%             if etime(clock,Timer.loopstart) > 20 && ~WarningNotRecordedColorChanged
%                 %set(thwarningnotrecorded(1), 'Color', WarningNotRecordedColor2);
%                 if length(thwarningnotrecorded) == 2
%                     set(thwarningnotrecorded(2), 'Color', WarningNotRecordedColor2);
%                 end
%                 WarningNotRecordedColorChanged = 1;
%             end
            
            if intersect(STATUSchan,ch)
                %if max(filteredvoltchunk(ch,:)) > TRIGthres
                if ~isempty(filteredvoltchunk)
                    bittest = bitand(uint16(max(filteredvoltchunk(ch,:))),TRIGbitvals)~=0;
                    if bittest(strcmpi(TRIGbittext,'TRIG'))
                        LastTrigTime = clock;
                    end
                end
                %if bittest(strcmpi(TRIGbittext,'BATT LOW'))
                %    BattLowWarn(ch) = 1;
                %else
                %    BattLowWarn(ch) = 0;
                %end
                %if bittest(strcmpi(TRIGbittext,'READY'))
                %    DeviceHasBeenReady = 1;
                %    DeviceErrorWarn(ch) = 0;
                %elseif DeviceHasBeenReady || etime(clock,Timer.loopstart) > 20
                %    DeviceErrorWarn(ch) = 1;
                %end
                for bt = length(bittest):-1:1
                    if bittest(bt)
                        T = TRIGbittext{bt};
                        break
                    end
                end
                HasChanges(ch) = 1;
            elseif intersect(DIAGchan,ch)
                ud = diff(rawshortsignal(ch,rawshortsignal(ch,:)>0));
                % normally it goes 1:2:63 then back to 1
                %mud = mod(ud,64);
                %mud( mud>32 ) = 64 - mud( mud>32 );
                %lostpackets = sum(mud(mud>2)/2-1);
                %errorfraction = lostpackets / length(mud);
                if strcmpi(systemMode, 'synfi')
                    lostpackets = sum(ud~=1);
                else
                    lostpackets = sum((mod( ud(ud ~= 2 & ud ~= -62), 64) - 2) / 2);
                end
                errorfraction = lostpackets / length(ud);
                %tquality = (nnz(ud==2) + nnz(ud==-62))/length(ud);
                %errorfraction = 1 - tquality;
                T = regexprep(sprintf('%4.1f%% ERR', min(100,errorfraction*100)),'100.0%','100%');
            else
                T = [num2str(V,'%6.2f') ' ' UNITname QTYname];
                ChannelStd(ch) = std(rawshortsignal(ch,:));
                HasChanges(ch) = any(diff(Vhistory(ch,:),[],2));
            end
            if any(diff(rawshortsignal(ch,:))) && HasChanges(ch)
                if abs(peaktopeakvoltages(ch,1)) < OverVoltageLevel || ismember(ch,METAchan) % Normal signals should be less than X microvolts
                    CC = NormalFontColor;
                    set(Kids(ch),'Color',NormalLineColor);
                    if ~ismember(ch, METAchan)
                        NumNormalChannels = NumNormalChannels + 1;
                    end
                else
                    if peaktopeakvoltages(ch,1) == max(peaktopeakvoltages)
                        CC = HighestFontColor;
                        set(Kids(ch),'Color',HighestFontColor);
                    else
                        CC = OverFontColor;
                        set(Kids(ch),'Color',OverLineColor);
                    end
                    if ~ismember(ch, METAchan)
                        NumNormalChannels = NumNormalChannels - Ndatachan*0.8;
                    end
                end
            else
                CC = FailedFontColor;
            end
        end
        while length(T) < 10
            T = [' ' T];
        end
        T2 = Chans{ch};
        while length(T2) < 4
            T2 = [' ' T2];
        end
        set(Grandkids(ch),'String',[T2 ' ' T],'Color',CC);
        if BattLowWarn(ch) && etime(clock,Timer.loopstart) >= AnalysisDuration+2.5
            set(Grandkids(ch), 'BackgroundColor', BattLowBackgroundColor, 'EdgeColor', BattLowEdgeColor, 'LineWidth', BattLowLineWidth);
        else
            set(Grandkids(ch), 'BackgroundColor', 'none', 'EdgeColor', 'none', 'LineWidth', 0.50);
        end
        
        if ismember(ch, DIAGchan) && etime(clock,Timer.loopstart) >= AnalysisDuration
            if errorfraction <= 0.05
                set(Grandkids(ch), 'Color', NormalFontColor);
            elseif errorfraction <= 0.25
                set(Grandkids(ch), 'Color', OverFontColor);
            else
                set(Grandkids(ch), 'Color', AngryFontColor);
            end
        end
    end
    
    if SWnohardwaremode
        if ~Warned.NoHardware
            Warned.NoHardware = 1;
        end
    end
    
    
    if any(BattLowWarn)
        if ~Warned.LowBattery
            %set(thwarninglowbattery, 'Visible', 'on');
            Warned.LowBattery = 1;
        end
    else
        if Warned.LowBattery
            %set(thwarninglowbattery, 'Visible', 'off');
            Warned.LowBattery = 0;
        end
    end
    
    FailedToChange = ChannelStd > 0 & HasChanges == 0;
    FailedToChange(METAchan) = 0;
    if nnz(FailedToChange) >= 2
        DeviceErrorWarn = FailedToStart | FailedToChange;
    else
        DeviceErrorWarn = FailedToStart;
    end
    
    if any(DeviceErrorWarn) && ~SWnohardwaremode
        if ~Warned.DeviceError
            WarningCount.DeviceError = WarningCount.DeviceError + 1;
            if WarningCount.DeviceError >= Deltasize
                %set(thwarningdeviceerror, 'Visible', 'on');
                Warned.DeviceError = 1;
            end
        end
    else
        if Warned.DeviceError
            WarningCount.DeviceError = 0;
            %set(thwarningdeviceerror, 'Visible', 'off');
            Warned.DeviceError = 0;
        end
    end
    
    drawnow;
%     rawenergy = signalpower(filteredvoltsignal.',SampleRate,[8 24]).';
%     for ch = 1:Nchan
%         C = min(rawenergy(ch,:),1);
%         set(Kids(ch),'Color',[0 C 0]);
%     end
    
    
    Timer.proctotal = etime(clock,Timer.procstart);
    elap = toc;
    if Deltasize
        FreeTime(mod(i-1,floor(Windowsize/Deltasize))+1) = RefreshPeriod-elap;
        FTsorted = sort(FreeTime(1:floor(Windowsize/Deltasize)));
        if length(FTsorted) > 1
            MinFreeTime = FTsorted(2);
            %MaxFreeTime = FTsorted(end-1);
        else
            MinFreeTime = FTsorted(1);
            %MaxFreeTime = FTsorted(end);
        end
        %AvgFreeTime = mean(FreeTime(1:Windowsize/Deltasize));
    else
        MinFreeTime = 0;
    end
    
    if SWautospeed
        if elap > RefreshPeriod && RefreshPeriod*2 <= RefreshPeriodAllow(2)
            %missed = missed+1;
            %if missed*RefreshPeriod > 2 || MinFreeTime < 0
            if MinFreeTime < 0
                fprintf('[%10.0f ms] Processing too slow. RefreshPeriod doubled.\n',etime(clock,Timer.loopstart)*1000);
                RefreshPeriod = RefreshPeriod*2;
                Deltasize = RefreshPeriod * SampleRate;
                %ocs = ones(1,Deltasize);
                %missed = 0;
                FreeTime = FreeTime*0;
            end
        elseif elap < RefreshPeriod/2 && RefreshPeriod/2 >= RefreshPeriodAllow(1)
            %missed = max(missed-1,0);
            %capaed = capaed+1;
            %if capaed*RefreshPeriod > 2 && MinFreeTime >= RefreshPeriod/2
            if MinFreeTime >= RefreshPeriod/2
                fprintf('[%10.0f ms] Refresh rate too slow. RefreshPeriod halved.\n',etime(clock,Timer.loopstart)*1000);
                RefreshPeriod = RefreshPeriod/2;
                Deltasize = RefreshPeriod * SampleRate;
                %ocs = ones(1,Deltasize);
                %capaed = 0;
                FreeTime = FreeTime*0;
            end
        end
    end
    if mod(i-1,10/RefreshPeriod) == 0
        fprintf('[%10.0f ms] Acq=%5.0f ms; Proc=%5.0f ms; Refresh=%5.0f ms; Freetime=%5.0f ms; Idle=%5.0f s\n',etime(clock,Timer.loopstart)*1000,Timer.acqtotal*1000,Timer.proctotal*1000,RefreshPeriod*1000,MinFreeTime*1000, etime(clock,Timer.lastactivity));
    end
    
    %20200701
    try
        if AutoNumericChanNames && mod(i-1,AutoNumericChanNames/RefreshPeriod) == 0
            set(scopecontrol.numericchannames, 'Enable', 'off');
            if mod(i-1,2*AutoNumericChanNames/RefreshPeriod) == 0
                setappdata(scopecontrol.numericchannames, 'ToEnableNumericChanNames', 1);
                setappdata(scopecontrol.numericchannames, 'ToDisableNumericChanNames', 0);
            else
                setappdata(scopecontrol.numericchannames, 'ToEnableNumericChanNames', 0);
                setappdata(scopecontrol.numericchannames, 'ToDisableNumericChanNames', 1);
            end
        end
    end
    
    BlinkOffTimeSec = 0.00; % max(0.50, RefreshPeriod);
    BlinkOnTimeSec = max(3.00-BlinkOffTimeSec, BlinkOffTimeSec*5);
    %blinkblinkwarning([WarnedNotRecorded, WarnedTriggerLost, WarnedLowBattery, WarnedBufferFull WarnedDeviceError], [thwarningnotrecorded; thwarningtriglost; thwarninglowbattery; thwarningbufferfull; thwarningdeviceerror], BlinkOnTimeSec, BlinkOffTimeSec);
    ERRORS = warningsystem(thwarning, Warned, WarnData, BlinkOnTimeSec, BlinkOffTimeSec, etime(clock,Timer.loopstart));
    PreAllocatedSize = size(mmfobj.Data.rawdata,2);
    updatetitle(Timer, CursorPosition, Deltasize, SampleRate, PreAllocatedSize, SWrecord, SWfilter, SWcar, NumNormalChannels, Ndatachan, DISPfullmedian, globalscale, [], [], titlehand, fighand, ERRORS, Fcutoff, Fcutoffhigh, SubjectName, [num2str(SampleRate) ' Hz'], SWzscore);
    if ~ishandle(fighand)
        Timer.lastactivity = clock;
        fprintf('[%10.0f ms] Recording stopped. Display window disappeared.\n',etime(clock,Timer.loopstart)*1000);
        break
    end
    if ~SWrecord && get(scopecontrol.startrecording, 'Value') && ~SWnohardwaremode
        Timer.lastactivity = clock;
        SWrecord = 1;
        Warned.NotRecorded = 0;
        RecordStart = CursorPosition + 1;
        RecordedSamples = 0;
        Timer.recordstart = clock;
        scopecontrol.startrecording_touched = clock;
        set(scopecontrol.startrecording, 'ForegroundColor', [0.4 1 0.4]);
        set(scopecontrol.startrecording, 'Enable', 'off');
        set(scopecontrol.restartnexus, 'Enable', 'off', 'String', RestartNexusButton.disabledtext);
        drawnow
        segmentsuid = datestr(now,'yyyymmdd-HHMMSS-FFF');
        SegmentSaveFileName = ['nexusautosave-' suid '_segment-' segmentsuid '.mat'];
        SegmentSaveUNC = [SaveDir filesep SegmentSaveFileName];
    elseif SWrecord && ~get(scopecontrol.startrecording, 'Value')
        Timer.lastactivity = clock;
        SWrecord = 0;
        Warned.NotRecorded = 1;
        RecordedSamples = 0;
        scopecontrol.startrecording_touched = clock;
        scopecontrol.restartnexus_touched = clock;
        set(scopecontrol.startrecording, 'ForegroundColor', [1 0.4 0.4]);
        set(scopecontrol.startrecording, 'Enable', 'off');
        drawnow
        nexus_savesegment(mmfilename, mmfobj, RecordStart, STATUSchan, SegmentSaveUNC, SampleRate, ChanNames);
    elseif get(scopecontrol.startrecording, 'Value') && SWnohardwaremode
        set(scopecontrol.startrecording, 'Value', 0);
    end
    
    if get(scopecontrol.paneltoggle, 'Value') && strcmp(get(scopecontrol.panel, 'Visible'), 'off')
        Timer.lastactivity = clock;
        set(scopecontrol.panel, 'Visible', 'on');
        set(scopecontrol.paneltoggle, 'BackgroundColor', [1 1 0]);
    elseif ~get(scopecontrol.paneltoggle, 'Value') && strcmp(get(scopecontrol.panel, 'Visible'), 'on')
        Timer.lastactivity = clock;
        set(scopecontrol.panel, 'Visible', 'off');
        set(scopecontrol.paneltoggle, 'BackgroundColor', [0 0 0]);
    end
    
    if ~SWcar && get(scopecontrol.car, 'Value')
        Timer.lastactivity = clock;
        SWcar = 1;
    elseif SWcar && ~get(scopecontrol.car, 'Value')
        Timer.lastactivity = clock;
        SWcar = 0;
    end
    
    if ~SWfilter && get(scopecontrol.filter, 'Value')
        Timer.lastactivity = clock;
        SWfilter = 1;
        if Fcutoffhigh > 0 && SWfilter && CSWmediandisp(1) == 1
            CSWmediandisp(1) = 0;
        else
            CSWmediandisp(1) = 1;
        end
    elseif SWfilter && ~get(scopecontrol.filter, 'Value')
        Timer.lastactivity = clock;
        SWfilter = 0;
        if Fcutoffhigh > 0 && SWfilter && CSWmediandisp(1) == 1
            CSWmediandisp(1) = 0;
        else
            CSWmediandisp(1) = 1;
        end
    end
    
    if ~SWnotch && get(scopecontrol.notchfilter, 'Value')
        Timer.lastactivity = clock;
        SWnotch = 1;
    elseif SWnotch && ~get(scopecontrol.notchfilter, 'Value')
        Timer.lastactivity = clock;
        SWnotch = 0;
    end

    
    if SWfilter
        tmp_1 = str2double(get(scopecontrol.filterhighpass, 'String'));
        tmp_2 = str2double(get(scopecontrol.filterlowpass, 'String'));
        if ~isfinite(tmp_1) || ~isfinite(tmp_2) || tmp_1 < 0 || tmp_1 >= tmp_2 || tmp_2 >= SampleRate/2
            set(scopecontrol.filterhighpass, 'String', num2str(Fcutoffhigh));
            set(scopecontrol.filterlowpass, 'String', num2str(Fcutoff));
        else
            Fcutoffhigh = tmp_1;
            Fcutoff = tmp_2;
            [FiltcoefB, FiltcoefA] = get_filter_coef(SampleRate, Fcutoffhigh, Fcutoff);
        end
    end
    
    if ~SWenvelope && get(scopecontrol.envelope, 'Value')
        Timer.lastactivity = clock;
        SWenvelope = 1;
        CSWmediandisp(2) = 1;
        SWzscore = 1;
        set(scopecontrol.zscore, 'Value', 1);
        globalscale = min(4/DISPfullmedian, globalscaleallow(2));
        set(scopecontrol.zoom, 'Value', log2(globalscale));
    elseif SWenvelope && ~get(scopecontrol.envelope, 'Value')
        Timer.lastactivity = clock;
        SWenvelope = 0;
        CSWmediandisp(2) = 0;
    end
    
    if SWenvelope
        tmp_2 = str2double(get(scopecontrol.envelopelowpass, 'String'));
        if ~isfinite(tmp_2) || tmp_2 <= 0 || tmp_2 >= SampleRate/2
            set(scopecontrol.envelopelowpass, 'String', num2str(Fenvelope));
        else
            Fenvelope = tmp_2;
            [EnvelopeFiltcoefB, EnvelopeFiltcoefA] = get_filter_coef(SampleRate, 0, Fenvelope);
        end
    end
    
    if ~SWautoscale && get(scopecontrol.autoscale, 'Value')
        Timer.lastactivity = clock;
        SWautoscale = 1;
    elseif SWautoscale && ~get(scopecontrol.autoscale, 'Value')
        Timer.lastactivity = clock;
        SWautoscale = 0;
    end
        
    if ~SWzscore && get(scopecontrol.zscore, 'Value')
        Timer.lastactivity = clock;
        SWzscore = 1;
        globalscale = min(4/DISPfullmedian, globalscaleallow(2));
        set(scopecontrol.zoom, 'Value', log2(globalscale));
    elseif SWzscore && ~get(scopecontrol.zscore, 'Value')
        Timer.lastactivity = clock;
        SWzscore = 0;
    end
    
    if SWplot && get(scopecontrol.noplot, 'Value')
        Timer.lastactivity = clock;
        SWplot = 0;
    elseif ~SWplot && ~get(scopecontrol.noplot, 'Value')
        Timer.lastactivity = clock;
        SWplot = 1;
        Zf = Zf_initial;
        notchZf = Zf_notch_initial;
        Zf_env = Zf_env_initial;
    end
    
    %set(scopecontrol.segmentbufferused, 'String', [num2str(RecordedSamples/SampleRate/60, '%.2f') ' min.']);
    %set(scopecontrol.segmentbufferavailable, 'String', [num2str((buffersize-RecordedSamples)/SampleRate/60, '%.2f') ' min.']);
    %set(scopecontrol.contigbufferused, 'String', [num2str(min(buffersize,CursorPosition)/SampleRate/60, '%.2f') ' min.']);
    %set(scopecontrol.contigbufferavailable, 'String', [num2str((buffersize-min(buffersize,CursorPosition))/SampleRate/60, '%.2f') ' min.']);
    set(scopecontrol.segmentbufferused, 'String', sec2hmsstr(RecordedSamples/SampleRate,2));
    set(scopecontrol.segmentbufferavailable, 'String', sec2hmsstr((buffersize-RecordedSamples)/SampleRate,2));
    set(scopecontrol.contigbufferused, 'String', sec2hmsstr(min(buffersize,CursorPosition)/SampleRate,2));
    set(scopecontrol.contigbufferavailable, 'String', sec2hmsstr((buffersize-min(buffersize,CursorPosition))/SampleRate,2));
    set(scopecontrol.refreshrate, 'String', [num2str(1/RefreshPeriod, '%.2f') ' Hz']);
    set(scopecontrol.samplerate, 'String', [num2str(SampleRate, '%.2f') ' Hz']);
    
    if etime(clock, scopecontrol.startrecording_touched) > 3.00
        set(scopecontrol.startrecording, 'Enable', 'on');
    end
    
    if getappdata(scopecontrol.reloadchannames, 'ToReloadChanNames')
        setappdata(scopecontrol.reloadchannames, 'ToReloadChanNames', 0);
        %fprintf('DEBUG: Button pressed\n');
        if exist(NexusChanNamesFileName, 'file')
            %fprintf('DEBUG: File exists\n');
            tmp = dir(NexusChanNamesFileName);
            if tmp.datenum > NexusChanNamesFileDate
                %fprintf('DEBUG: File newer\n');
                % This is a newer version of the channames file
                NexusChanNamesFileDate = tmp.datenum;
                tmp2 = load(NexusChanNamesFileName, 'NexusChanNames');
                %fprintf('DEBUG: File loaded\n');
                if length(ChanNames) == length(tmp2.NexusChanNames.synfi)
                    ChanNames = tmp2.NexusChanNames.synfi;
                    %OriginalChanNames = ChanNames;
                elseif length(ChanNames) == length(tmp2.NexusChanNames.nexus1)
                    ChanNames = tmp2.NexusChanNames.nexus1;
                    %OriginalChanNames = ChanNames;
                end
                [ChanNames, Chans] = populate_channames (ChanNames, Nchan, acqchan, DIAGchan, NchanPerRow, TwoColumnMode, subplothand);
            end
        end
    end
    
    if getappdata(scopecontrol.numericchannames, 'ToEnableNumericChanNames')
        setappdata(scopecontrol.numericchannames, 'ToEnableNumericChanNames', 0);
        NumericChanNames = ChanNames;
        for j = 1:length(NumericChanNames)
            if j < 34
                NumericChanNames{j} = sprintf('1-%02i',j);
            else
                NumericChanNames{j} = sprintf('2-%02i',j-33);
            end
        end
        [NumericChanNames, Chans] = populate_channames (NumericChanNames, Nchan, acqchan, DIAGchan, NchanPerRow, TwoColumnMode, subplothand);
    elseif getappdata(scopecontrol.numericchannames, 'ToDisableNumericChanNames')
        setappdata(scopecontrol.numericchannames, 'ToDisableNumericChanNames', 0);
        [ChanNames, Chans] = populate_channames (ChanNames, Nchan, acqchan, DIAGchan, NchanPerRow, TwoColumnMode, subplothand);
    end

    
    if getappdata(scopecontrol.restartnexus, 'ToRestartNexus')
        setappdata(scopecontrol.restartnexus, 'ToRestartNexus', 0);
        scopecontrol.restartnexus_touched = clock;
        if RestartNexusConfirm == 0
            RestartNexusConfirm = 1;
            set(scopecontrol.restartnexus, 'Enable', 'off', 'String', RestartNexusButton.confirmtext);
        elseif RestartNexusConfirm == 1
            RestartNexusConfirm = 0;
            set(scopecontrol.restartnexus, 'Enable', 'off', 'String', RestartNexusButton.actiontext);
            if ~SWrecord && ~SWnohardwaremode
                restart_nexus(Nchan, SampleRate, AcqDuration);
                scopecontrol.restartnexus_touched = clock;
            end
        end
    end
    
    if ~SWrecord
        if etime(clock, scopecontrol.restartnexus_touched) > 5.00
            set(scopecontrol.restartnexus, 'Enable', 'on', 'String', RestartNexusButton.normaltext);
            RestartNexusConfirm = 0;
        elseif RestartNexusConfirm == 1 && etime(clock, scopecontrol.restartnexus_touched) > 0.50
            set(scopecontrol.restartnexus, 'Enable', 'on');
        end
    end
    
    pause(RefreshPeriod - elap);

    
    if etime(clock,Timer.loopstart) > 120 && etime(clock,Timer.lastactivity) > InactivityTimeoutSec
        % If no user activity over this long, automatically stop the acquisition to allow Nexus to turn off.
        fprintf('[%10.0f ms] Stopping. No operator activity for %g seconds.\n',etime(clock,Timer.loopstart)*1000, InactivityTimeoutSec);
        close(fighand);
    end
    
    if ~ishandle(fighand)
        Timer.lastactivity = clock;
        fprintf('[%10.0f ms] Recording stopped. Display window disappeared.\n',etime(clock,Timer.loopstart)*1000);
        break
    end
    
    prevgcfpos = gcfpos;
    gcfpos = get(gcf, 'Position');
    if any(prevgcfpos ~= gcfpos)
        Timer.lastactivity = clock;
    end
    
    
    
end

if SWrecord
    if ishandle(fighand)
        set(thwarning, 'Visible', 'off');
        T = 'STOPPED';
        CC = SavedFontColor;
        while length(T) < 10
            T = [' ' T];
        end
        for ch = 1:Nchan
            T2 = Chans{ch};
            while length(T2) < 4
                T2 = [' ' T2];
            end
            set(Grandkids(ch),'String',[T2 ' ' T],'Color',CC);
        end
        drawnow
    end
end

if ~SWnohardwaremode
    nexus_exit();
end

fprintf('Memory mapped file: %s \n',mmfilename);

if SWrecord
    nexus_savesegment(mmfilename, mmfobj, RecordStart, STATUSchan, SegmentSaveUNC, SampleRate, ChanNames);
end

rawdata = mmfobj.Data.rawdata(:,mmfobj.Data.rawdata(end,:)>0);
rawevent = mmfobj.Data.rawevent(:,mmfobj.Data.rawdata(end,:)>0);
[~, i1] = min(rawdata(end,:));
[~, i2] = max(rawdata(end,:));
if i2 < i1
    rawdata = rawdata(1:end,[i1:end, 1:i2]);
    rawevent = rawevent(1:end,[i1:end, 1:i2]);
else
    rawdata = rawdata(1:end,i1:i2);
    rawevent = rawevent(1:end,i1:i2);
end
fprintf('Autosaving the memory mapped file. This is separate from the segment save.\n');
save(SaveUNC, 'rawdata', 'rawevent', 'SampleRate', 'ChanNames', 'OriginalChanNames', '-v7.3');
fprintf('Autosave file: %s \n',SaveUNC);
fprintf('done.\n');
rawdata = rawdata(:,rawdata(STATUSchan(1),:)>0);
rawdata = rawdata(:,rawdata(end,:) >= RecordStart);
rawdata = rawdata(1:end-1,:);


    function f_fig_keypress (hObject, eventdata)
        Timer.lastactivity = clock;
        Key = eventdata.Key;
        if ~isempty(regexp(Key, '^numpad\d', 'match', 'once'))
            Key = Key(7:end);
        end
        if ~isempty(regexp(Key, '([0-9a-z])|(space)|(return)', 'match', 'once'))
            
            % Check if there are existing keys recorded in the same time frame
            a = CursorPosition+1;
            b = CursorPosition+Deltasize;
            L = length(a:b);
            if ~isempty(a:b)
                ma = mod(a-1, buffersize)+1;
                mb = mod(b-1, buffersize)+1;
                if mb >= ma
                    Event = mmfobj.Data.rawevent(1,ma:mb);
                else
                    % Wrap around required
                    Event = mmfobj.Data.rawevent(1,[ma:buffersize,1:mb]);
                end
                Event = Event(Event>0);
            else
                Event = [];
            end
            
            if length(Key) == 1
                Event = [Event uint8(Key)];
            elseif strcmp(Key, 'space')
                Event = [Event uint8(32)];
            elseif strcmp(Key, 'return')
                Event = [Event uint8(13)];
            end
            
            if length(Event) > L
                Event = Event(1:L);
            else
                Event = [Event zeros(1,L - length(Event))];
            end
            
            % Write the event back
            if ~isempty(a:b)
                if mb >= ma
                    mmfobj.Data.rawevent(1,ma:mb) = Event;
                else
                    % Wrap around required
                    mmfobj.Data.rawevent(1,[ma:buffersize,1:mb]) = Event;
                end
            end
        end

        
    end


    function nexus_exit ()
        try
            nexus_unload;
            try
                assignin('base', 'naa_tnexusactive', tic);
            end
            fprintf('DAQ unloaded.\n');
        catch
            fprintf('Failed to unload DAQ.\n');
        end
    end


    function blinkblinkwarning (WarnStates, thwarninghandles, OnTimeSec, OffTimeSec)
        persistent blinktic
        if isempty(blinktic)
            blinktic = tic;
        end
        TotalSec = OnTimeSec + OffTimeSec;
        Tock = mod(toc(blinktic), TotalSec);
        if Tock <= OnTimeSec
            TurnOn = 1;
        else
            TurnOn = 0;
        end
        
        
        for i = 1:length(WarnStates)
            if WarnStates(i)
                for j = 1:size(thwarninghandles,2)
                    try
                        if strcmp(get(thwarninghandles(i,j), 'Visible'), 'on')
                            if ~TurnOn
                                set(thwarninghandles(i,j), 'Visible', 'off');
                            end
                        elseif strcmp(get(thwarninghandles(i,j), 'Visible'), 'off')
                            if TurnOn
                                set(thwarninghandles(i,j), 'Visible', 'on');
                            end
                        end
                    end
                end
            end
        end
        
        drawnow
    end

    function ERRORS = warningsystem (thwarning, Warned, WarnData, OnTimeSec, OffTimeSec, ElapsedTime)
        
        % persistent blinktic
        % if isempty(blinktic)
        %     blinktic = tic;
        % end
        % TotalSec = OnTimeSec*sum(a) + OffTimeSec;
        % Tock = mod(toc(blinktic), TotalSec);
        % if Tock <= OnTimeSec
        %     TurnOn = 1;
        % else
        %     TurnOn = 0;
        % end
        
        ERRORS = '';
        
        if ~ishandle(thwarning)
            return
        end
        
        a = struct2cell(Warned);
        a = cat(2,a{:});
        if any(a)
            set(thwarning, 'Visible', 'on');
        else
            set(thwarning, 'Visible', 'off');
            return
        end
        
        
        TotalSec = OnTimeSec*sum(a) + OffTimeSec;
        ind = floor(mod(ElapsedTime, TotalSec) / OnTimeSec) + 1;
        if mod(ElapsedTime, TotalSec) > OnTimeSec*sum(a)
            set(thwarning, 'Visible', 'off');
            return
        end
        
        fn = fieldnames(Warned);
        k = 0;
        for i = 1:length(fn)
            if Warned.(fn{i})
                k = k + 1;
                if k == ind
                    set(thwarning, 'String', WarnData.(fn{i}).String);
                    set(thwarning, 'FontSize', WarnData.(fn{i}).FontSize);
                    set(thwarning, 'Color', WarnData.(fn{i}).Color1);
                    if ElapsedTime > 20
                        set(thwarning, 'Color', WarnData.(fn{i}).Color2);
                    end
                    set(thwarning, 'BackgroundColor', WarnData.(fn{i}).BackgroundColor);
                    set(thwarning, 'EdgeColor', WarnData.(fn{i}).EdgeColor);
                    set(thwarning, 'LineWidth', WarnData.(fn{i}).LineWidth);
                    set(thwarning, 'FontWeight', WarnData.(fn{i}).FontWeight);
                    set(thwarning, 'HorizontalAlignment', WarnData.(fn{i}).HorizontalAlignment);
                    set(thwarning, 'VerticalAlignment', WarnData.(fn{i}).VerticalAlignment);
                    set(thwarning, 'FontName', WarnData.(fn{i}).FontName);
                    ERRORS = WarnData.(fn{i}).String;
                    break
                end
            end
        end
        drawnow
    end


    function updatetitle (Timer, CursorPosition, Deltasize, SampleRate, PreAllocatedSize, SWrecord, SWfilter, SWcar, NumNormalChannels, Ndatachan, DISPfullmedian, globalscale, ~, ~, titlehand, fighand, ERRORS, Fcutoff, Fcutoffhigh, SUBJECT, SAMPRATE, SWzscore)
        LoopElapTimeStr = [num2str(etime(clock,Timer.loopstart)*1,'%10.0f') ' sec'];
        NRwarning = '';
        TRIGstatus = '';
        FIwarning = '';
        NOTCHwarning = '';
        CAwarning = '';
        GCcount = '';
        if ~SWrecord
            NRwarning = '(Not Recorded)';
        else
            NRwarning = sprintf('(Buffer=%.0f/%.0f min)', (CursorPosition+Deltasize)/SampleRate/60, PreAllocatedSize/SampleRate/60);
        end
        if SWfilter
            FIwarning = ['(BPF ' num2str(Fcutoffhigh,'%.2f') '-' num2str(Fcutoff,'%.2f') 'Hz)'];
        end
        if SWnotch
            NOTCHwarning = ['(NOTCH ' num2str(PowerLineFrequency,'%.2f') 'Hz)'];
        end
        if SWcar
            CAwarning = '(CAR)';
        end
        if (NumNormalChannels/Ndatachan) >= 0.25
            GCcount = ['(Stable: ' num2str(NumNormalChannels) '/' num2str(Ndatachan) ')'];
        else
            GCcount = '(Unstable)';
        end
        
        % trigs = bitand(uint16(rawshortsignaltrig),TRIGbitvals(strcmpi(TRIGbittext,'TRIG')));
        % trigs = sum(trigs,1);
        % bb = getdigitalbounds(trigs.');
        % if ~isempty(bb)
        %     pulsewidth = mean((diff(bb,[],2)+1))/SampleRate
        %     iti = mean((diff(bb(:,1),[],1)))/SampleRate
        % end
        
        [UNITfactor, UNITname] = get_best_unitfactor (DISPfullmedian*globalscale);
        if SWzscore
            UNITfactor = 1;
            UNITname = 'Zscores';
        end
        
        TitleStr = ['Expm time: ' LoopElapTimeStr '; ' 'Scale: ' num2str(DISPfullmedian*globalscale*UNITfactor) ' ' UNITname ' ' GCcount FIwarning NOTCHwarning CAwarning NRwarning TRIGstatus];
        if ~strcmp(TitleStr,PrevTitleStr)
            if min(ishandle(titlehand)) && ishandle(fighand)
                set(titlehand,'String',TitleStr);
                if mod(floor(now*86400/2),2)
                    % 2016-04-11 Po: Alternate between title and errors
                    set(fighand, 'Name', ['nexus_scope ' SUBJECT ' ' SAMPRATE ' ' ERRORS ' ' GCcount]);
                else
                    set(fighand, 'Name', [ERRORS ' ' 'nexus_scope ' SUBJECT ' ' SAMPRATE ' ' GCcount]);
                end
                
            end
        end
        PrevTitleStr = TitleStr;
        
        
        % Check for event too
        if ishandle(EventIndicatorLight)
            a = CursorPosition+1-3*SampleRate;
            b = CursorPosition+1;
            if ~isempty(a:b)
                ma = mod(a-1, buffersize)+1;
                mb = mod(b-1, buffersize)+1;
                if mb >= ma
                    Event = mmfobj.Data.rawevent(1,ma:mb);
                else
                    % Wrap around required
                    Event = mmfobj.Data.rawevent(1,[ma:buffersize,1:mb]);
                end
                Event = Event(Event>0);
            else
                Event = [];
            end
            if isempty(Event)
                set(EventIndicatorLight, 'Visible', 'off');
                set(EventIndicatorText, 'String', '');
            else
                set(EventIndicatorLight, 'Visible', 'on');
                set(EventIndicatorText, 'String', char(Event));
            end
        end
        
    end



    function nexus_savesegment (mmfilename, mmfobj, RecordStart, STATUSchan, SegmentSaveUNC, SampleRate, ChanNames)
        fprintf('Memory mapped file: %s \n',mmfilename);
        fprintf('Saving this segment to MATLAB format.\n');
        % rawdata = mmfobj.Data.rawdata(:,mmfobj.Data.rawdata(end,:)>=RecordStart);
        % [~, i1] = min(rawdata(end,:));
        % [~, i2] = max(rawdata(end,:));
        % if i2 < i1
        %     rawdata = rawdata(1:end,[i1:end, 1:i2]);
        % else
        %     rawdata = rawdata(1:end,i1:i2);
        % end
        % rawdata = rawdata(1:end-1,rawdata(STATUSchan(1),:)>0);
        % save(SegmentSaveUNC, 'rawdata', 'SampleRate', 'ChanNames', '-v7.3');
        
        % 2016-04-11 Po: New save format is [time points, channels] to facilitate
        % faster saving and loading
        rawsignaldata = permute(mmfobj.Data.rawdata(:,mmfobj.Data.rawdata(end,:)>=RecordStart),[2 1]);
        rawevent = permute(mmfobj.Data.rawevent(:,mmfobj.Data.rawdata(end,:)>=RecordStart),[2 1]);
        [~, i1] = min(rawsignaldata(:,end));
        [~, i2] = max(rawsignaldata(:,end));
        if i2 < i1
            rawsignaldata = rawsignaldata([i1:end, 1:i2],:);
            rawevent = rawevent([i1:end, 1:i2],:);
        else
            rawsignaldata = rawsignaldata(i1:i2,:);
            rawevent = rawevent(i1:i2,:);
        end
        % 20180503: Bug fix: Swapped these two lines
        rawevent = rawevent(rawsignaldata(:,STATUSchan(1))>0,:); 
        rawsignaldata = rawsignaldata(rawsignaldata(:,STATUSchan(1))>0,1:end-1); 
        save(SegmentSaveUNC, 'rawsignaldata', 'rawevent', 'SampleRate', 'ChanNames', '-v7.3');
        
        fprintf('Saved this segment to file: %s \n',SegmentSaveUNC);
    end



    function [FiltcoefB, FiltcoefA, nlookback] = get_filter_coef(SampleRate, Fcutoffhigh, Fcutoff)
        FiltcoefB = 1;
        FiltcoefA = 1;
        
        if Fcutoffhigh > 0
            [FB2,FA2] = butter(HPF_Order,Fcutoffhigh/(SampleRate/2),'high');
            FiltcoefB = conv(FiltcoefB,FB2);
            FiltcoefA = conv(FiltcoefA,FA2);
        end
        
        tmp = Fcutoff/(SampleRate/2);
        if tmp < 1
            [FB2,FA2] = butter(LPF_Order,tmp,'low');
            FiltcoefB = conv(FiltcoefB,FB2);
            FiltcoefA = conv(FiltcoefA,FA2);
        end
        
        nlookback = max(length(FiltcoefA),length(FiltcoefB)) - 1;
    end




    function [ChanNames, Chans] = populate_channames (ChanNames, Nchan, acqchan, DIAGchan, NchanPerRow, TwoColumnMode, subplothand)
        if ~isempty(who('ChanNames')) && ~isempty(ChanNames) && iscell(ChanNames)
            if length(ChanNames) < Nchan
                tmp = length(ChanNames);
                tmp = tmp + 1;
                while tmp <= Nchan
                    ChanNames{tmp} = num2str(acqchan(tmp));
                    tmp = tmp + 1;
                end
            end
            Chans = ChanNames(1:Nchan);
            Chans(DIAGchan) = {''};
            for i = 1:length(Chans)
                YLabel{i} = Chans{i};
            end
            if ~TwoColumnMode
                set(gca,'YTickLabel',fliplr(YLabel));
            else
                for i = 1:2
                    set(subplothand(i),'YTickLabel',fliplr(YLabel((i-1)*NchanPerRow+1:(i)*NchanPerRow)));
                end
            end
        end
    end



    function scopecontrol_zoom_in (hObject, eventdata)
        Timer.lastactivity = clock;
        data = getappdata(hObject, 'data');
        ax = data.gca;
        ylmax = data.ylim;
        yl = get(ax, 'YLim');
        c = (yl(1)+yl(2))/2;
        r = (yl(2)-yl(1))/2;
        r = r / sqrt(2);
        yl = c + [-r, r];
        if yl(1) < ylmax(1), y1(1) = ylmax(1); end
        if yl(2) > ylmax(2), yl(2) = ylmax(2); end
        set(ax, 'YLim', yl);
    end


    function scopecontrol_zoom_out (hObject, eventdata)
        Timer.lastactivity = clock;
        data = getappdata(hObject, 'data');
        ax = data.gca;
        ylmax = data.ylim;
        yl = get(ax, 'YLim');
        c = (yl(1)+yl(2))/2;
        r = (yl(2)-yl(1))/2;
        r = r * sqrt(2);
        yl = c + [-r, r];
        if yl(1) < ylmax(1), yl(1) = ylmax(1); end
        if yl(2) > ylmax(2), yl(2) = ylmax(2); end
        set(ax, 'YLim', yl);
    end


    function scopecontrol_pan_up (hObject, eventdata)
        Timer.lastactivity = clock;
        data = getappdata(hObject, 'data');
        ax = data.gca;
        ylmax = data.ylim;
        yl = get(ax, 'YLim');
        yl = yl + 1;
        set(ax, 'YLim', yl);
    end


    function scopecontrol_pan_down (hObject, eventdata)
        Timer.lastactivity = clock;
        data = getappdata(hObject, 'data');
        ax = data.gca;
        ylmax = data.ylim;
        yl = get(ax, 'YLim');
        yl = yl - 1;
        set(ax, 'YLim', yl);
    end


    function scopecontrol_reload_channames (hObject, eventdata)
        Timer.lastactivity = clock;
        setappdata(hObject, 'ToReloadChanNames', 1);
    end


    function scopecontrol_numeric_channames (hObject, eventdata)
        Timer.lastactivity = clock;
        if get(hObject, 'value')
            setappdata(hObject, 'ToEnableNumericChanNames', 1);
            setappdata(hObject, 'ToDisableNumericChanNames', 0);
        else
            setappdata(hObject, 'ToEnableNumericChanNames', 0);
            setappdata(hObject, 'ToDisableNumericChanNames', 1);
        end
    end


    function scopecontrol_restart_nexus (hObject, eventdata)
        Timer.lastactivity = clock;
        setappdata(hObject, 'ToRestartNexus', 1);
    end


    function restart_nexus (Nchan, SampleRate, AcqDuration)
        disp('Attempting to restart Nexus');
        
        try
            nexus_stop();
        catch
            disp('Error stopping');
        end
        
        try
            nexus_init(Nchan);
        catch exception
            disp(['Error: Error initialzing. ' exception.message]);
            return
        end
        
        try
            nexus_start(SampleRate, AcqDuration);
        catch exception
            disp(['Error: Error starting data acquisition. ' exception.message]);
            return
        end
    end

    function fakedata = fakenexus_getsimulateddata (AcqSamplRateHz, AcqTimeSeconds, NumChannels, libnames, ForwardFetch)
        fakedata = zeros(NumChannels, AcqSamplRateHz*AcqTimeSeconds);
        tt = (0:AcqSamplRateHz*AcqTimeSeconds-1)/AcqSamplRateHz;
        %mainAmp = 10;
        mainAmp = mod(now*86400,60);
        reo = regexp(FakeChanNames, '^(\d+)Hz$', 'tokens', 'once');
        reoN1 = regexp(FakeChanNames, '^(\d+)HzN1$', 'tokens', 'once');
        reoN2 = regexp(FakeChanNames, '^(\d+)HzN2$', 'tokens', 'once');
        reoN3 = regexp(FakeChanNames, '^(\d+)HzN3$', 'tokens', 'once');
        reoN4 = regexp(FakeChanNames, '^(\d+)HzN4$', 'tokens', 'once');
        reoN5 = regexp(FakeChanNames, '^(\d+)HzN5$', 'tokens', 'once');
        reoN6 = regexp(FakeChanNames, '^(\d+)HzN6$', 'tokens', 'once');
        reoN7 = regexp(FakeChanNames, '^(\d+)HzN7$', 'tokens', 'once');
        reoSquare = regexp(FakeChanNames, '^(\d+)HzSquare$', 'tokens', 'once');
        targetHz = zeros(length(FakeChanNames),1);
        for ich = 1:length(FakeChanNames)
            if ~isempty(reo{ich})
                targetHz(ich) = str2double(reo{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
            elseif ~isempty(reoN1{ich})
                targetHz(ich) = str2double(reoN1{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 1 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoN2{ich})
                targetHz(ich) = str2double(reoN2{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 10 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoN3{ich})
                targetHz(ich) = str2double(reoN3{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 100 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoN4{ich})
                targetHz(ich) = str2double(reoN4{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 1000 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoN5{ich})
                targetHz(ich) = str2double(reoN5{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 10000 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoN6{ich})
                targetHz(ich) = str2double(reoN6{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 100000 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoN7{ich})
                targetHz(ich) = str2double(reoN7{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = mainAmp * sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                realHzN = PowerLineFrequency * (1+randn/1000);
                fakedata(ich,:) = fakedata(ich,:) + 1000000 * sin(2*pi*tt*realHzN) .* (1+randn(1,length(tt))/100); % add power line noise
            elseif ~isempty(reoSquare{ich})
                targetHz(ich) = str2double(reoSquare{ich}{1});
                realHz = targetHz(ich) * (1+randn/100);
                fakedata(ich,:) = sin(2*pi*tt*realHz) .* (1+randn(1,length(tt))/100);
                fakedata(ich,fakedata(ich,:)<0) = -mainAmp;
                fakedata(ich,fakedata(ich,:)>0) = mainAmp;
            end
        end
        
        fakedata(STATUSchan,:) = 2;
        fakedata(DIAGchan,:) = fakecounter + (1:AcqSamplRateHz*AcqTimeSeconds);
        if ~ForwardFetch
            pause(AcqTimeSeconds);
        else
            pause(AcqTimeSeconds/10);
        end
    end

end