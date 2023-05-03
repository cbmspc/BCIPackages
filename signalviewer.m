% This program expects Signal to be in microvolts and SampleRate to be in
% samples per second
%
% Optional input arguments:
% First, create a struct called opts. Every field is optional.
%
%   opts.EventTimeStamps should be a Nx2 cell {1.2345, 'Abc'; 4.5, 'Def'} with first column in seconds, and 2nd column a string
%
%   opts.ica_W = ICA separating matrix with orientation (Nsource x Nchan), i.e. ica_W * Signal.' = Source.'
%   opts.ica_A = ICA mixing matrix with orientation (Nchan x Nsource), i.e. ica_A * Source.' = Signal.'
%    Example: [ica_sig, ica_A, ica_W] = fastica(Signal.', 'stabilization', 'on', 'maxNumIterations', 200);
%    If ica matrices are not specified, FastICA is used
%
%   opts.FooterMessage = 'footer messaage'
%
%   opts.pwelch_nfft = 4096
%   NFFT for pwelch
%
%   opts.stitch_mult = 2
%   Stitching algorithm does a high pass on each segment before joining them
%
%   opts.blankaround_stitch_samples = 9
%   This many samples to the left and to the right of the signal are NaN'd to
%   hide the stitching artifact for display only; does not affect filtering. 
%
%   opts.notch = 60
%   If specified, turns on the notch filter at the specified powerline frequency (and harmonics)
%
%   opts.bandpass = [1 35]
%   If specified, turns on both high pass and low pass filters at the specified corner frequencies
%
%   opts.highpass = 1
%   If specified, turns on high pass filter at the specified corner frequency. This field is ignored if bandpass is specified.
%
%   opts.lowpass = 30
%   If specified, turns on low pass filter at the specified corner frequency. This field is ignored if bandpass is specified.
%
%   opts.sensitivity = 100
%   If specified, sets vertical sensitivity (in microvolts) after viewer is opened
%   By default, the viewer autoscales after applying filters. This can be disabled by setting opts.sensitivity = -1
%

function argout1 = signalviewer(Signal, SampleRate, ChanNames, opts, argin5, argin6, argin7)
%t_program_start = tic;
try
    mfnp = mfilename('fullpath');
    if isempty(regexp(mfnp, '\\.m$', 'match'))
        mfnp = [mfnp '.m'];
    end
    tmp = dir(mfnp);
    lmoddate = datestr(tmp.datenum, 'YYmmdd.HHMM');
    clear tmp
catch
    lmoddate = '';
end


if ~iscell(Signal) && size(Signal,3) > 1
    % Assume ch x time x trial
    tmp = cell(size(Signal,3));
    for i = 1:size(Signal,3)
        tmp{i} = Signal(:,:,i).';
    end
    Signal = tmp;
end

pwelch_nfft = [];

stitch_mult = 2;
% Stitching algorithm does a high pass on each segment before joining them
% To turn the stitching algorithm off, pass in opts.stitch_mult = 0

blankaround_stitch_samples = round(0.035 * SampleRate);
% This many samples to the left and to the right of the signal are NaN'd to
% hide the stitching artifact for display only; does not affect filtering. 
% Example: blankaround_stitch_samples = 15
% To turn the blanking off, pass in opts.blankaround_stitch_samples = 0

set_chansep_to = 0;
% After the viewer is opened, scale to this vertical sensitivity



if nargin >= 4 && exist('opts', 'var') && isstruct(opts)
    if isfield(opts, 'EventTimeStamps')
        EventTimeStamps = opts.EventTimeStamps;
    elseif isfield(opts, 'eventtimestamps')
        EventTimeStamps = opts.eventtimestamps;
    else
        clear EventTimeStamps;
    end
    if isfield(opts, 'ica_W')
        ica_W = opts.ica_W;
    elseif isfield(opts, 'ica_w')
        ica_W = opts.ica_w;
    else
        ica_W = [];
    end
    if isfield(opts, 'ica_A')
        ica_A = opts.ica_A;
    elseif isfield(opts, 'ica_a')
        ica_A = opts.ica_a;
    else
        ica_A = [];
    end
    if isfield(opts, 'FooterMessage')
        FooterMessage = opts.FooterMessage;
    elseif isfield(opts, 'footermessage')
        FooterMessage = opts.footermessage;
    else
        FooterMessage = '';
    end
    if isfield(opts, 'pwelch_nfft')
        pwelch_nfft = opts.pwelch_nfft;
    else
        pwelch_nfft = [];
    end
    if isfield(opts, 'stitch_mult') && isnumeric(opts.stitch_mult) && numel(opts.stitch_mult) == 1 && isfinite(opts.stitch_mult) && opts.stitch_mult >= 0
        stitch_mult = opts.stitch_mult;
    end
    if isfield(opts, 'blankaround_stitch_samples') && isnumeric(opts.blankaround_stitch_samples) && numel(opts.blankaround_stitch_samples) == 1 && isfinite(opts.blankaround_stitch_samples) && opts.blankaround_stitch_samples >= 0
        blankaround_stitch_samples = opts.blankaround_stitch_samples;
    end
    if isfield(opts, 'sensitivity') && isnumeric(opts.sensitivity) && numel(opts.sensitivity) == 1 && isfinite(opts.sensitivity) && opts.sensitivity >= 0
        set_chansep_to = opts.sensitivity;
    end
elseif nargin >= 4
    % Old callers still use the 4th input argument as EventTimeStamps
    
    if exist('opts', 'var') && ~isstruct(opts)
        EventTimeStamps = opts;
    end
    if exist('argin5', 'var')
        ica_W = argin5;
    end
    if exist('argin6', 'var')
        ica_A = argin6;
    end
    if exist('argin7', 'var')
        FooterMessage = argin7;
    end
end

SignalIsStitched = 0;
if iscell(Signal)
    UserSuppliedEventNames = 0;
    if exist('EventTimeStamps','var') && isnumeric(EventTimeStamps)
        CellEvents = EventTimeStamps;
        UserSuppliedEventNames = 1;
    else
        CellEvents = 1:length(Signal);
    end
    EventTimeStamps = zeros(length(Signal),3);
    tmp = 0;
    for i = 1:size(EventTimeStamps,1)
        EventTimeStamps(i,:) = [tmp, tmp + size(Signal{i},1) / SampleRate, CellEvents(i)];
        tmp = tmp + size(Signal{i},1) / SampleRate;
    end
    
    % 2022-05-24: Keep the time points in each cell before merging
    SignalIsStitched = 1;
    EventTimePoints = EventTimeStamps(:,1:2) * SampleRate;
    EventTimePoints(:,1) = EventTimePoints(:,1) + 1;
    
    if ~UserSuppliedEventNames
        EventTimeStamps = EventTimeStamps(:,1);
    end
    
    
    Signal = cat(1,Signal{:});
else
    EventTimePoints = [];
end


Signal = double(Signal);


if ~exist('ChanNames', 'var') || isempty(ChanNames)
    if size(Signal,1) == 1 && size(Signal,2) > 1
        %Signal = Signal.';
        Signal = permute(Signal,[2 1 3 4 5 6 7 8 9]);
    elseif size(Signal,2) >= SampleRate && size(Signal,1) < SampleRate / 2
        %Signal = Signal.';
        Signal = permute(Signal,[2 1 3 4 5 6 7 8 9]);
    end
    
    nchan = size(Signal,2);
    npad = floor(log10(nchan))+1;
    ChanNames = string_to_cell(num2str(1:nchan,['ch%0' num2str(npad) 'i,']),',');
end
ChanNames = ChanNames(:).';
if size(Signal,1) == length(ChanNames) && size(Signal,2) ~= length(ChanNames)
    Signal = permute(Signal,[2 1 3 4 5 6 7 8 9]);
end
Fs = SampleRate;

% 2022-05-24
Fcut_minfreq = Fs/10000;

averagesegmentduration = [];

StitchSignalCell();




if ~exist('ica_W', 'var') || ~exist('ica_A', 'var')
    ica_W = [];
    ica_A = [];
end

if ~exist('FooterMessage', 'var')
    FooterMessage = '';
end

if ~exist('opts', 'var') || ~isstruct(opts)
    opts = struct;
end

viewhand_ica_sig = ceil(rand*1000000000);
viewhand_ica_A = ceil(rand*1000000000);
viewhand_ica_W = ceil(rand*1000000000);
viewhand_psd = ceil(rand*1000000000);
viewhand_psd_axe = -1;
selected_plothand = -1;
selected_timepoint = -1;
previously_selected_timepoint = -2;
FineSnapScale = 100;
selected_cursortype = 0;
CursorEnable = 0;

% Remove channels that are completely flat
%n = size(Signal,1);
% if n > 491520
%     % Sample only these many points in 4 areas
%     chnc = true(1,size(Signal,2));
%     b = floor(n/4)*(0:3)' * [1 1] + ones(4,1)*[1 122880];
%     for ch = 1:size(Signal,2) %#ok<*FXUP>
%         if chnc(ch)
%             for i = 1:size(b,1)
%                 chnc(ch) = chnc(ch) & nanstd(Signal(b(i,1):b(i,2),ch),[],1) == 0;
%                 if ~chnc(ch)
%                     break
%                 end
%             end
%         end
%     end
% else
%     chnc = nanstd(Signal,[],1) == 0;
% end


% Identify not-connected channels by signal variance
chnc = nanmax(Signal) - nanmin(Signal) == 0;

% Identify not-connected channels by zero-length channel names
if length(chnc) == length(ChanNames)
    chnc = chnc | cellfun(@isempty,ChanNames);
end

if exist('ica_W', 'var') && exist('ica_A', 'var') && ~isempty(ica_W) && ~isempty(ica_A)
    if size(ica_W,2) == size(Signal,2) && size(ica_A,1) == size(Signal,2)
        ica_W = ica_W(:,~chnc);
        ica_A = ica_A(~chnc,:);
    end
else
    ica_W = [];
    ica_A = [];
end
Signal = Signal(:,~chnc);
ChanNames = ChanNames(~chnc);


% Fix NaNs if any
Signal = fix_nans_for_filtering(Signal);



selchan = 1:length(ChanNames);
selica = [];


if isempty(ChanNames)
    error('There is nothing to plot. Either the input is an empty matrix or all channels are completely flat.');
end


%[ChanNames,ix] = sort(ChanNames);
%Signal = Signal(:,ix);
%ica_W = ica_W(:,ix);
%ica_A = ica_A(ix,:);
Signal_postica = Signal;
Signal_postnotch = Signal;
Signal_postfilter = Signal;
Signal4 = Signal(:,1);
Signal_psd_source = Signal(:,1) * 0;

% If a signal is known to be bandlimited due to low-pass and/or envelope
% filter, reduce the number of points to plot to speed up
SigBandwidth = Fs/2;
%BandLimitedInterpolation = 1;
BLIM = 1;

% If there are more signal points than horizontal screen resolution, reduce
% the number of points to plot to speed up
ScreenLimitedDownsampling = 1;
SLD_H = 32768;
DefaultYColor = [0.15 0.15 0.15];
BusyYColor = [1 0 0];
InactiveYColor = [0.65 0.65 0.65];


PlotHold = 0;


PowerLineFrequency = 60;
MaxFilterOrder = 4;
FilterOrder = MaxFilterOrder;
FilterReflectMult = 0;
FilterChunkSec = 600;
FilterBusy = 0;
MovementBusy = 0;
ICA_Initialized = 0;
ica_sig = [];
icachans = {};

Nsch = length(ChanNames(selchan));
Ntp = size(Signal,1);

PermittedXZoomRanges   = [0.001 0.005 0.01 0.05 0.1 0.5 1 2 5 10 20 30 60 120 300 600 1200 1800 3600:3600:6*3600 8*3600 12*3600 24*3600 7*24*3600];
PermittedXZoomRanges = PermittedXZoomRanges(PermittedXZoomRanges > 1/Fs);
PermittedChanSepRanges = [1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000];
XTickSpacingsAndUnits = {
    1e-6    1e-6    'us'
    2e-6    1e-6    'us'
    10e-6   1e-6    'us'
    20e-6   1e-6    'us'
    100e-6  1e-6    'us'
    200e-6  1e-6    'us'
    1e-3    0.001   'ms'
    2e-3    0.001   'ms'
    10e-3   0.001   'ms'
    20e-3   0.001   'ms'
    0.1     1   's'
    0.2     1   's'
    0.5     1   's'
    1       1   's'
    2       1   's'
    5       1   's'
    10      1   's'
    20      1   's'
    30      1   's'
    60      60  'm'
    300     60  'm'
    600     60  'm'
    1200    60  'm'
    1800    60  'm'
    3600    3600    'h'
    7200    3600    'h'
    3600*4  3600    'h'
    3600*8  3600    'h'
    3600*24 86400    'd'
};
XTickSpacings = cat(1,XTickSpacingsAndUnits{:,1});

tmp = ceil(Nsch/3)*3;
Kolor = jet(tmp)*0.50;
tmp = reshape(reshape(rearrange_top_bottom((1:tmp).'),tmp/3,[]).',1,[]);
Kolor = Kolor(tmp,:);

Nkolor = size(Kolor,1);

EventKolor = [0.75 0.75 0];
EventFontSize = 14;
EventFontWeight = 'bold';

AxesFontName = 'Consolas';
AxesFontSize = 16;
CursorTextFontSize = 8;

AxesPosition = [0.0300    0.0600    0.88    0.93];

chansep = 100;

screensize = get(groot,'Screensize');
tmp = mod(ceil(unixtimemillis()),604800000)+100000001;
while ishandle(tmp)
    tmp = tmp + 1;
end
fighand = figure(tmp);
if nargout > 0
    argout1 = fighand;
end
clf
set(fighand, 'ToolBar', 'none', 'MenuBar', 'none');
set(fighand, 'CloseRequestFcn', @f_main_close);
if SignalIsStitched
    tmp = sprintf('Signal is stitched from %g segments with %.2g s average length', size(EventTimePoints,1), averagesegmentduration);
else
    tmp = '';
end
set(fighand, 'Name', ['Signal Viewer: ' num2str(size(Signal,2)) ' channels, ' num2str(size(Signal,1)/Fs) ' seconds record duration, ' num2str(Fs) ' Hz sample rate. ' tmp]);
hold on
set(fighand,'Position',screensize);
YLim = [-chansep*Nsch-0.5*chansep, -chansep+0.5*chansep];
XLim = [0 10];
set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim, 'XLim', XLim, 'Position', AxesPosition, 'FontWeight', 'bold', 'FontName', AxesFontName, 'FontUnits', 'pixel', 'FontSize', AxesFontSize); %#ok<*NBRAK>
axehand = get(fighand,'CurrentAxes');

t1 = 1;
t2 = 2;
Time = (0:Ntp-1)/Fs;
Time_min = 0;
Time_max = (Ntp-1)/Fs;

cursorlinehand = plot([0 0], [0 0], '-', 'LineWidth', 2, 'Color', [0.8 0.8 0.8]);
for ch = Nsch:-1:1
    plotpip1hand(ch) = plot(0, 0, '+', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 24, 'LineWidth', 2, 'Visible', 'off');
    plotpip2hand(ch) = plot(0, 0, '+', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 24, 'LineWidth', 2, 'Visible', 'off');
    plotpip3hand(ch) = plot(0, 0, '+', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 24, 'LineWidth', 2, 'Visible', 'off');
    plothand(ch) = plot(Time(t1:t2), Signal(t1:t2,selchan(ch)) - nanmean(Signal(t1:t2,selchan(ch))) - chansep*ch);
    plottext1hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', CursorTextFontSize);
    plottext2hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', CursorTextFontSize);
    plottext3hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', CursorTextFontSize);
    set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:));
    set(plothand(ch), 'ButtonDownFcn', @f_plothand_buttondown);
    setappdata(plothand(ch), 'chanind', selchan(ch));
    setappdata(plothand(ch), 'channame', ChanNames{selchan(ch)});
end
if ~isempty(plothand)
    selected_plothand = plothand(1);
end

set(axehand, 'FontUnits', 'normalized');
AxesFontSize = get(axehand, 'FontSize');

tmax = max(Time);
SigChunkRendered = false(ceil(tmax / FilterChunkSec),1);
ChunkIndexMax = length(SigChunkRendered);

% Allow different kinds of event time stamp formats
if exist('EventTimeStamps','var') && ~isempty(EventTimeStamps)
    if isstruct(EventTimeStamps)
        % event_time format
        tmp = EventTimeStamps;
        EventTimeStamps = {};
        tmpf = fieldnames(tmp);
        k = 0;
        for i = 1:length(tmpf)
            if isnumeric(tmp.(tmpf{i}))
                for j = 1:length(tmp.(tmpf{i}))
                    k = k + 1;
                    EventTimeStamps(k,:) = {tmp.(tmpf{i})(j), tmpf{i}}; %#ok<*AGROW>
                end
            end
        end
    elseif isnumeric(EventTimeStamps) && size(EventTimeStamps,2) == 3
        % boundary format with augmented label
        tmp = EventTimeStamps;
        EventTimeStamps = {};
        k = 0;
        for i = 1:size(tmp,1)
            k = k + 1;
            EventTimeStamps(k,:) = {tmp(i,1), sprintf('Event %g Start', tmp(i,3))};
            k = k + 1;
            EventTimeStamps(k,:) = {tmp(i,2), sprintf('Event %g End', tmp(i,3))};
        end
    elseif isnumeric(EventTimeStamps) && size(EventTimeStamps,2) == 2
        % boundary format without augmented label
        tmp = EventTimeStamps;
        EventTimeStamps = {};
        k = 0;
        for i = 1:size(tmp,1)
            k = k + 1;
            EventTimeStamps(k,:) = {tmp(i,1), sprintf('Event %i Start', i)};
            k = k + 1;
            EventTimeStamps(k,:) = {tmp(i,2), sprintf('Event %i End', i)};
        end
    elseif isnumeric(EventTimeStamps) && size(EventTimeStamps,2) == 1
        % simple time stamps with no labels
        tmp = EventTimeStamps;
        EventTimeStamps = {};
        k = 0;
        for i = 1:size(tmp,1)
            k = k + 1;
            EventTimeStamps(k,:) = {tmp(i,1), sprintf('Event %i', i)};
        end
    end
    
    
end

if exist('EventTimeStamps','var') && ~isempty(EventTimeStamps) && iscell(EventTimeStamps) && size(EventTimeStamps,2) == 2
    EventEnable = 1;
else
    EventEnable = 0;
end


if EventEnable
    try
        EventTimeStamps = sortrows(EventTimeStamps);
    catch
        [B,I] = sortrows(EventTimeStamps(:,1));
        EventTimeStamps = [B EventTimeStamps(I,2)];
        clear B I
    end
    EventTimes = cell2mat(EventTimeStamps(:,1));
    YPos = sort(mean(YLim)+diff(YLim)/20*([-8:2:8]), 'descend');
    NYPos = length(YPos);
    for i = size(EventTimeStamps,1):-1:1
        eventplothand(i) = plot( EventTimeStamps{i,1}*[1 1], [-10000000*(Nsch+1), 10000000], '-', 'Color', EventKolor ); %#ok<NASGU>
        eventtexthand(i) = text( EventTimeStamps{i,1}, YPos(mod(i-1,NYPos)+1), EventTimeStamps{i,2}, 'FontUnits', 'pixel', 'FontSize', EventFontSize, 'FontWeight', EventFontWeight);
    end
    set(eventtexthand(i), 'Visible', 'off');
end


HighPassFilter.state = 0;
%HighPassFilter.cutoff = ceil(Fs/99)/100; % Fs/10000 likely the stability limit
HighPassFilter.cutoff = ceil(Fcut_minfreq*101)/100;
LowPassFilter.state = 0;
LowPassFilter.cutoff = floor(Fs/2-1);
EnvelopeFilter.state = 0;
EnvelopeFilter.cutoff = 2;
NotchFilter.state = 0;
NotchFilter.order = 4; % Must be even
NotchFilter.qfactor = 10;
ZscoreFilter.state = 0;
ZscoreFilter.multiplier = 100;
ZscoreFilter.off_chansep = chansep;
ZscoreFilter.on_chansep = [];

%text_off = 'OFF';
%text_on = 'ON';
fontweight_on1 = 'bold';
fontweight_on2 = 'bold';
fontweight_off = 'normal';
fontcolor_on1 = [0 0.5 0];
fontcolor_on2 = [1 0 0];
fontcolor_off = [0 0 0];



h_hugetext = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.35 0.5 0.25 0.08], 'String', 'BIG TEXT', 'FontUnits', 'normalized', 'FontSize', 0.8, 'Visible', 'off', 'BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1]);
h_bigtext = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.25 0.5 0.45 0.05], 'String', 'BIG TEXT', 'FontUnits', 'normalized', 'FontSize', 0.4, 'Visible', 'off');

NormalizedControlFontSize = 0.8;

h_xoomtitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.96 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'Horizontal Zoom', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_xzoomout = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.94 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_xzoomlevel = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.94 0.04 0.02], 'BackgroundColor', [0.7 0.7 0.9], 'String', '10 s', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_xzoomin = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.975 0.94 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_yoomtitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.91 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', '# chans on screen', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_yzoomout = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.89 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_yzoomlevel = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.89 0.04 0.02], 'BackgroundColor', [0.7 0.7 0.9], 'String', '64 ch', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_yzoomin = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.975 0.89 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_sensititle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.86 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'Y Sensitivity', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_sepup = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.84 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_sensitivity = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.84 0.04 0.02], 'BackgroundColor', [0.7 0.7 0.9], 'String', '100 uV', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_sepdown = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.975 0.84 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_pantitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.81 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'Pan', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_panleft = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.95 0.775 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '<', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_panright = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.775 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '>', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_panup = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.80 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '^', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_pandown = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.775 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'v', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_hold_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.92 0.75 0.045 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Pause Plotting', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8);
%h_hold_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.965 0.74 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_psd_plot = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.70 0.03 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'PSD', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8);
h_autofit = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.70 0.03 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Auto Fit', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.5);


h_passfilttitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.68 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Filters', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_notch_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.92 0.66 0.034 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', ['Notch ' num2str(PowerLineFrequency) 'Hz'], 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
%h_notch_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.937 0.66 0.017 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_notch_order = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.954 0.66 0.010 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(NotchFilter.order), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_notch_orderunit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.964 0.66 0.010 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ord', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_notch_qfactor = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.975 0.66 0.010 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(NotchFilter.qfactor), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_notch_qfactorunit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.985 0.66 0.005 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Q', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>


h_hpf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.92 0.64 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'HighPass', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_hpf_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.64 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_hpf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.64 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(HighPassFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_hpf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.64 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_lpf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.92 0.62 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'LowPass', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_lpf_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.62 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_lpf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.62 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(LowPassFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_lpf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.62 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_evf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.92 0.60 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Envelope', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_evf_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.60 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_evf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.60 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(EnvelopeFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_evf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.60 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_zscore_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.92 0.58 0.030 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Zscore', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
%h_zscore_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.94 0.58 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_fastdraw_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.96 0.58 0.030 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Fast', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Value', 1, 'FontWeight', fontweight_on1, 'ForegroundColor', fontcolor_on1);
%h_fastdraw_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.98 0.58 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_on, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_chansel_title = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.525 0.030 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Plotted Channels', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.5); %#ok<NASGU>
h_chansel_list = uicontrol(fighand, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.92 0.12 0.040, 0.400], 'FontUnits', 'pixel');
h_chansel_confirm = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.10 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Plot', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_chansel_auto = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.94 0.10 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Auto', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_chansel_reset = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.945 0.10 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'R', 'FontUnits', 'normalized', 'Visible', 'off');

h_cursor_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.94 0.08 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Cursor', 'FontUnits', 'normalized', 'FontName', 'Times New Roman');
%h_cursor_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.95 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'off', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_signal_export = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.06 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Export Sig', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);


h_icasel_title = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.96 0.525 0.030 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ICA Comps', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.5); %#ok<NASGU>
h_icasel_list = uicontrol(fighand, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.96 0.12 0.030, 0.400], 'FontUnits', 'pixel');
h_icasel_confirm = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.10 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Start', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_reset = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.10 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'R', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_sources = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'S', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_mixmat = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_sepmat = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_export = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.06 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Export ICA', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);

h_axesfont_inc = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.04 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_axesfont_dec = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.93 0.04 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_windowhsize_inc = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.95 0.04 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_windowhsize_dec = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.04 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);

h_xspan_text = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.020 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'full t range [0,12345]', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);

h_xspan_edittext1intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.020 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'XLim(1) =', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_xspan_edit1 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.04 0.020 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_xspan_edittext1unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.07 0.020 0.025 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'seconds', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>

h_lmoddate = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.10 0.020 0.14 0.012], 'BackgroundColor', [0.7 0.7 0.7], 'String', ['SignalViewer version: ' lmoddate], 'FontUnits', 'normalized', 'FontSize', 1.0, 'FontName', 'FixedWidth', 'HorizontalAlignment', 'left'); %#ok<NASGU>
h_hintbar = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.25 0.020 0.57 0.018], 'BackgroundColor', [0.7 0.9 0.7], 'String', '', 'FontUnits', 'normalized', 'FontSize', 0.8, 'FontName', 'Verdana', 'HorizontalAlignment', 'left'); %#ok<NASGU>

h_xspan_edittext2intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.83 0.020 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'XLim(2) =', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_xspan_edit2 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.86 0.020 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_xspan_edittext2unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.89 0.020 0.025 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'seconds', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>

% h_yspan_edittext1intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.41 0.025 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Selected channel YMin =', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edit1 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.44 0.025 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edittext1unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.47 0.025 0.025 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'µV', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% 
% h_yspan_edittext2intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.50 0.025 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Selected channel YMax =', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edit2 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.53 0.025 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edittext2unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.56 0.025 0.025 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'µV', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);


h_footer_message = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.10 0.020 0.70 0.020], 'BackgroundColor', [0.8 0.8 0.8], 'String', FooterMessage, 'FontUnits', 'normalized', 'FontSize', 16, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
if isempty(FooterMessage)
    set(h_footer_message, 'Visible', 'off');
end

set(h_icasel_reset, 'Enable', 'off');
set(h_icasel_view_sources, 'Enable', 'off');
set(h_icasel_view_mixmat, 'Enable', 'off');
set(h_icasel_view_sepmat, 'Enable', 'off');
set(h_icasel_view_export, 'Enable', 'off');
set(h_chansel_list, 'String', ChanNames);
set(h_chansel_list, 'Value', selchan);

set(h_xzoomin, 'Callback', @f_xzoomin);
set(h_xzoomout, 'Callback', @f_xzoomout);
set(h_yzoomin, 'Callback', @f_yzoomin);
set(h_yzoomout, 'Callback', @f_yzoomout);
set(h_panleft, 'Callback', @f_panleft);
set(h_panright, 'Callback', @f_panright);
set(h_panup, 'Callback', @f_panup);
set(h_pandown, 'Callback', @f_pandown);
set(h_sepup, 'Callback', @f_sepup);
set(h_sepdown, 'Callback', @f_sepdown);

set(h_xspan_edit1, 'Callback', @f_xspan_edit1);
set(h_xspan_edit2, 'Callback', @f_xspan_edit2);

set(h_hold_switch, 'Callback', @f_hold_switch);

set(h_hpf_switch, 'Callback', @f_hpf_switch);
set(h_lpf_switch, 'Callback', @f_lpf_switch);
set(h_hpf_cutoff, 'Callback', @f_hpf_cutoff);
set(h_lpf_cutoff, 'Callback', @f_lpf_cutoff);
set(h_evf_switch, 'Callback', @f_evf_switch);
set(h_evf_cutoff, 'Callback', @f_evf_cutoff);
set(h_notch_switch, 'Callback', @f_notch_switch);
set(h_notch_order, 'Callback', @f_notch_order);
set(h_notch_qfactor, 'Callback', @f_notch_qfactor);
set(h_zscore_switch, 'Callback', @f_zscore_switch);
set(h_fastdraw_switch, 'Callback', @f_fastdraw_switch);
set(h_chansel_list, 'Callback', @f_chansel_list);
set(h_chansel_confirm, 'Callback', @f_chansel_confirm);
set(h_chansel_reset, 'Callback', @f_chansel_reset);
set(h_psd_plot, 'Callback', @f_psd_plot);
set(h_autofit, 'Callback', @f_autofit);
set(h_icasel_confirm, 'Callback', @f_icasel_confirm);
set(h_icasel_reset, 'Callback', @f_icasel_reset);
set(h_icasel_view_sources, 'Callback', @f_icasel_view_sources);
set(h_icasel_view_mixmat, 'Callback', @f_icasel_view_mixmat);
set(h_icasel_view_sepmat, 'Callback', @f_icasel_view_sepmat);
set(h_icasel_view_export, 'Callback', @f_icasel_view_export);
set(h_chansel_auto, 'Callback', @f_chansel_auto);
set(h_cursor_switch, 'Callback', @f_cursor_enable);
set(h_signal_export, 'Callback', @f_signal_export);
set(h_axesfont_inc, 'Callback', @f_axesfont_inc);
set(h_axesfont_dec, 'Callback', @f_axesfont_dec);
set(h_windowhsize_inc, 'Callback', @f_windowhsize_inc);
set(h_windowhsize_dec, 'Callback', @f_windowhsize_dec);
set(fighand, 'KeyPressFcn', @f_fig_keypress);
set(fighand,'Position',screensize);

notch_update();
filter_update();
%set(h_bigtext, 'Visible', 'off', 'String', '');

set(h_xspan_text, 'String', ['full t range [' num2str(round(min(Time))) ', ' num2str(round(max(Time))) '] s']);

try
    pause(0.00001);
    oldWarningState = warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');
    frame_h = get(handle(fighand),'JavaFrame');
    set(frame_h,'Maximized',1);
    warning(oldWarningState);
end


if isfield(opts, 'notch')
    if ~isempty(opts.notch) && opts.notch > 0
        if opts.notch >= 40 && opts.notch <= 300
            PowerLineFrequency = opts.notch;
        end
        f_notch_switch(fighand, []);
    end
end
if isfield(opts, 'bandpass')
    if numel(opts.bandpass) == 2 && opts.bandpass(1) > 0 && opts.bandpass(2) < Fs/2 && opts.bandpass(1) < opts.bandpass(2)
        set(h_hpf_cutoff, 'String', num2str(opts.bandpass(1)));
        set(h_lpf_cutoff, 'String', num2str(opts.bandpass(2)));
        f_bpf_switch(fighand, []);
    end
else
    if isfield(opts, 'highpass')
        if numel(opts.highpass) == 1 && opts.highpass(1) > 0 && opts.highpass(1) < Fs/2
            set(h_hpf_cutoff, 'String', num2str(opts.highpass(1)));
            f_hpf_switch(fighand, []);
        end
    end
    if isfield(opts, 'lowpass')
        if numel(opts.lowpass) == 1 && opts.lowpass(1) > 0 && opts.lowpass(1) < Fs/2
            set(h_lpf_cutoff, 'String', num2str(opts.lowpass(1)));
            f_lpf_switch(fighand, []);
        end
    end
end

%drawnow

set(fighand, 'HandleVisibility', 'callback');

if set_chansep_to >= 0
    autofit();
end
if set_chansep_to > 0
    refit(set_chansep_to);
end






    function f_fig_keypress(hObject, eventdata)
        Key = eventdata.Key;
        if ismember('control', eventdata.Modifier)
            Ctrl = 1;
        else
            Ctrl = 0;
        end
        if ismember('alt', eventdata.Modifier)
            Alt = 1;
        else
            Alt = 0;
        end
        if ismember('shift', eventdata.Modifier)
            Shift = 1;
        else
            Shift = 0;
        end
        
        switch Key
            case 'g'
                uin = inputdlg('Center on specific time in seconds:');
                gotosec = str2double(uin);
                if isfinite(gotosec)
                    if FilterBusy
                        return;
                    end
                    XRange = XLim(2)-XLim(1);
                    XLim = gotosec + [-1 1]*XRange/2;
                    set(axehand, 'XLim', XLim);
                    resnap_pan();
                end
            case 'leftarrow'
                if Ctrl
                    f_xzoomout(hObject, []);
                elseif Alt && EventEnable
                    u = find(EventTimes < XLim(1), 1, 'last');
                    if ~isempty(u)
                        XRange = XLim(2) - XLim(1);
                        XLim(1) = EventTimes(u) - XRange/2;
                        XLim(2) = EventTimes(u) + XRange/2;
                        resnap_pan();
                    elseif XLim(1) < min(EventTimes) && XLim(1) > 0
                        f_panleft(hObject, 5.0);
                    end
                elseif Alt && ~EventEnable
                    if XLim(1) > 0
                        f_panleft(hObject, 5.0);
                    end
                elseif Shift
                    f_panleft(hObject, 0.10);
                else
                    f_panleft(hObject, []);
                end
            case 'rightarrow'
                if Ctrl
                    f_xzoomin(hObject, []);
                elseif Alt && EventEnable
                    u = find(EventTimes > XLim(2), 1, 'first');
                    if ~isempty(u)
                        XRange = XLim(2) - XLim(1);
                        XLim(1) = EventTimes(u) - XRange/2;
                        XLim(2) = EventTimes(u) + XRange/2;
                        resnap_pan();
                    elseif XLim(2) > max(EventTimes) && XLim(2) < Time_max
                        f_panright(hObject, 5.0);
                        %fprintf('a');
                    end
                elseif Alt && ~EventEnable
                    if XLim(2) < Time_max
                        f_panright(hObject, 5.0);
                    end
                elseif Shift
                    f_panright(hObject, 0.10);
                else
                    f_panright(hObject, []);
                end
            case 'uparrow'
                if Ctrl
                    f_sepdown(hObject, []);
                else
                    f_panup(hObject, []);
                end
            case 'downarrow'
                if Ctrl
                    f_sepup(hObject, []);
                else
                    f_pandown(hObject, []);
                end
            case 'pageup'
                f_sepdown(hObject, []);
            case 'pagedown'
                f_sepup(hObject, []);
            case 'home'
                f_yzoomin(hObject, []);
            case 'end'
                f_yzoomout(hObject, []);
            case 'insert'
                f_xzoomin(hObject, []);
            case 'delete'
                f_xzoomout(hObject, []);
                
        end
        set(h_hintbar, 'String', 'Ctrl left/right: change time scale. Ctrl up/down: change sensitivity. Shift left/right: scroll slowly. Alt left/right: go to events.');
    end


    function f_hold_switch(hObject, eventdata) %#ok<*INUSD>
        if PlotHold
            %set(h_hold_state, 'String', text_off);
            set(h_hold_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            PlotHold = 0;
            set(h_hugetext, 'Visible', 'off');
            redraw();
            set(h_hintbar, 'String', '');
        else
            %set(h_hold_state, 'String', text_on);
            set(h_hold_switch, 'Value', 1, 'ForegroundColor', fontcolor_on2, 'FontWeight', fontweight_on2);
            set(h_hugetext, 'String', 'Plotting Paused', 'Visible', 'on');
            PlotHold = 1;
            set(h_hintbar, 'String', 'Plotting paused');
        end
    end



    function f_hpf_switch(hObject, eventdata)
        disable_filter_switches();
        if ~HighPassFilter.state
            HighPassFilter.state = 1;
        else
            HighPassFilter.state = 0;
        end
        
        if HighPassFilter.state
            %set(h_hpf_state, 'String', text_on);
            set(h_hpf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
        else
            %set(h_hpf_state, 'String', text_off);
            set(h_hpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        
        filter_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_lpf_switch(hObject, eventdata)
        disable_filter_switches();
        if ~LowPassFilter.state
            LowPassFilter.state = 1;
        else
            LowPassFilter.state = 0;
        end
        
        if LowPassFilter.state
            %set(h_lpf_state, 'String', text_on);
            set(h_lpf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
        else
            %set(h_lpf_state, 'String', text_off);
            set(h_lpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        
        filter_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_bpf_switch(hObject, eventdata)
        % This is not a real switch but rather lets a batch script update
        % two filters at the same time
        disable_filter_switches();
        HighPassFilter.state = 1;
        %set(h_hpf_state, 'String', text_on);
        set(h_hpf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
        LowPassFilter.state = 1;
        %set(h_lpf_state, 'String', text_on);
        set(h_lpf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
        filter_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_evf_switch(hObject, eventdata)
        disable_filter_switches();
        if ~EnvelopeFilter.state
            EnvelopeFilter.state = 1;
        else
            EnvelopeFilter.state = 0;
        end
        
        if EnvelopeFilter.state
            %set(h_evf_state, 'String', text_on);
            set(h_evf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
        else
            %set(h_evf_state, 'String', text_off);
            set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        
        filter_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_hpf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if HighPassFilter.state
            filter_update();
        end
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_lpf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if LowPassFilter.state
            filter_update();
        end
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_evf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if EnvelopeFilter.state
            filter_update();
        end
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function disable_filter_switches()
        set([h_hpf_switch, h_lpf_switch, h_evf_switch, h_hpf_cutoff, h_lpf_cutoff, h_evf_cutoff, h_notch_switch, h_notch_order, h_notch_qfactor], 'Enable', 'off');
        drawnow
    end

    function enable_filter_switches()
        set([h_hpf_switch, h_lpf_switch, h_evf_switch, h_hpf_cutoff, h_lpf_cutoff, h_evf_cutoff, h_notch_switch, h_notch_order, h_notch_qfactor], 'Enable', 'on');
    end

    function disable_movement_switches()
        set([h_panleft h_panright h_panup h_pandown h_xzoomout h_xzoomin h_yzoomout h_yzoomin h_sepup h_sepdown], 'Enable' , 'off');
    end

    function enable_movement_switches()
        set([h_panleft h_panright h_panup h_pandown h_xzoomout h_xzoomin h_yzoomout h_yzoomin h_sepup h_sepdown], 'Enable' , 'on');
    end

    function f_notch_switch(hObject, eventdata)
        disable_filter_switches();
        %set(h_notch_state, 'String', 'Wait'); drawnow;
        if NotchFilter.state
            NotchFilter.state = 0;
        else
            NotchFilter.state = 1;
        end
        notch_update();
        filter_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_notch_order(hObject, eventdata)
        disable_filter_switches();
        ord = str2double(get(h_notch_order, 'String'));
        if isfinite(ord) && mod(ord,2) == 0 && ord >= 2 && ord <= 10
            %set(h_notch_state, 'String', 'Wait'); 
            drawnow;
            NotchFilter.order = ord;
            notch_update();
            filter_update();
        else
            set(h_notch_order, 'String', num2str(NotchFilter.order));
        end
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_notch_qfactor(hObject, eventdata)
        disable_filter_switches();
        q = str2double(get(h_notch_qfactor, 'String'));
        if isfinite(q) && q > 0
            %set(h_notch_state, 'String', 'Wait'); 
            drawnow;
            NotchFilter.qfactor = q;
            notch_update();
            filter_update();
        else
            set(h_notch_qfactor, 'String', num2str(NotchFilter.qfactor));
        end
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    
    function notch_update()
        if NotchFilter.state
            FilterBusy = 1;
            set(h_bigtext, 'Visible', 'on', 'String', ['Preparing notch filter...']); drawnow;
            clear d Hd
            Funda = PowerLineFrequency;
            for h = 1:floor(Fs/2/Funda)
                % Create one for each harmonic
                d = fdesign.notch('N,F0,Q',NotchFilter.order,Funda*h/(Fs/2),NotchFilter.qfactor*h);
                Hd{h} = design(d);
            end
            warning('off', 'signal:filtfilt:ParseSOS');
            warning('off', 'signal:filtfilt:ParseB');
            for ch = size(Signal_postica,2):-1:1
                set(h_bigtext, 'Visible', 'on', 'String', ['Applying notch filter (' num2str(ch) ' chans to go)']); drawnow;
                %Signal_postnotch(:,ch) = freqfilter(Signal_postica(:,ch), Fs, [PowerLineFrequency+[-2 2], NotchOrder], 'stop', 'butter', 1*Fs);
                Signal_postnotch(:,ch) = Signal_postica(:,ch);
                for h = 1:length(Hd)
                    Signal_postnotch(:,ch) = filtfilt(Hd{h}.sosMatrix,Hd{h}.ScaleValues,Signal_postnotch(:,ch));
                end
            end
            set(h_bigtext, 'Visible', 'on', 'String', ['Finishing notch filter...']); drawnow;
            %set(h_notch_state, 'String', text_on);
            set(h_notch_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', '');
        else
            Signal_postnotch = Signal_postica;
            %set(h_notch_state, 'String', text_off);
            set(h_notch_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
    end


    function f_zscore_switch(hObject, eventdata)
        set(h_zscore_switch, 'Enable', 'off');
        %set(h_zscore_state, 'String', 'Wait'); 
        drawnow;
        if ZscoreFilter.state
            ZscoreFilter.on_chansep = chansep;
            ZscoreFilter.state = 0;
            %set(h_zscore_state, 'String', text_off);
            set(h_zscore_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            refit(ZscoreFilter.off_chansep);
        else
            ZscoreFilter.off_chansep = chansep;
            ZscoreFilter.state = 1;
            %set(h_zscore_state, 'String', text_on);
            set(h_zscore_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            if isempty(ZscoreFilter.on_chansep)
                autofit();
            else
                refit(ZscoreFilter.on_chansep);
            end
        end
        redraw();
        set(h_zscore_switch, 'Enable', 'on');
    end


    function f_fastdraw_switch(hObject, eventdata)
        set(h_fastdraw_switch, 'Enable', 'off');
        %set(h_fastdraw_state, 'String', 'Wait'); 
        drawnow;
        if ScreenLimitedDownsampling
            ScreenLimitedDownsampling = 0;
            %set(h_fastdraw_state, 'String', text_off);
            set(h_fastdraw_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        else
            ScreenLimitedDownsampling = 1;
            %set(h_fastdraw_state, 'String', text_on);
            set(h_fastdraw_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
        end
        redraw();
        set(h_fastdraw_switch, 'Enable', 'on');
    end


    function filter_update()
        hpf = str2double(get(h_hpf_cutoff, 'String'));
        lpf = str2double(get(h_lpf_cutoff, 'String'));
        evf = str2double(get(h_evf_cutoff, 'String'));
        
        % 20220711: Set them to impossible values but still valid
        % numerics to avoid crashing.
        if ~isfinite(hpf)
            hpf = 0;
        end
        if ~isfinite(lpf)
            lpf = Fs;
        end
        
        if ~isfinite(lpf)
            evf = Fs;
        end
        
        if hpf <= Fcut_minfreq % Not likely to be stable, or not even available
            HighPassFilter.state = 0;
            set(h_hpf_cutoff, 'String', num2str(ceil(Fcut_minfreq*101)/100));
            %set(h_hpf_state, 'String', text_off);
            set(h_hpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        if lpf >= Fs/2
            LowPassFilter.state = 0;
            set(h_lpf_cutoff, 'String', num2str(floor(Fs/2*99)/100));
            %set(h_lpf_state, 'String', text_off);
            set(h_lpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        if hpf >= lpf
            HighPassFilter.state = 0;
            %set(h_hpf_state, 'String', text_off);
            set(h_hpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            LowPassFilter.state = 0;
            %set(h_lpf_state, 'String', text_off);
            set(h_lpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        if evf >= Fs/2
            EnvelopeFilter.state = 0;
            set(h_evf_cutoff, 'String', num2str(floor(Fs/2*99)/100));
            %set(h_evf_state, 'String', text_off);
            set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        
        if HighPassFilter.cutoff ~= hpf || LowPassFilter.cutoff ~= lpf
            HighPassFilter.cutoff = hpf;
            LowPassFilter.cutoff = lpf;
        end
        if EnvelopeFilter.cutoff ~= evf
            EnvelopeFilter.cutoff = evf;
        end
        if HighPassFilter.state || LowPassFilter.state || EnvelopeFilter.state
            SigChunkRendered(:) = 0;
            FilterOrder = MaxFilterOrder;
            filter_render();
        else
            Signal_postfilter = Signal_postnotch;
            SigChunkRendered(:) = 1;
        end
        
        redraw();
    end


    function chunkindexrange = xlim_to_sigchunk()
        chunkindexrange = [max(1,floor(XLim(1)/FilterChunkSec)+1), min(ChunkIndexMax,ceil(XLim(2) / FilterChunkSec))];
    end


    function filter_render()
        hpf = str2double(get(h_hpf_cutoff, 'String'));
        lpf = str2double(get(h_lpf_cutoff, 'String'));
        evf = str2double(get(h_evf_cutoff, 'String'));
        chunkindexrange = xlim_to_sigchunk();
        if SigChunkRendered(chunkindexrange(1):chunkindexrange(2))
            %redraw();
            return;
        end
        
        if HighPassFilter.state || LowPassFilter.state || EnvelopeFilter.state
            FilterBusy = 1;
            disable_movement_switches();
            %set([h_hpf_state h_lpf_state h_evf_state], 'String', 'Wait');
            set(h_bigtext, 'Visible', 'on', 'String', 'Preparing filters...'); 
            %redraw();
            drawnow;
            
            % Only need to render the chunks and 1/8 chunk to the left and
            % to the right
            
            
            %e1 = max(Time_min, (chunkindexrange(1)-1-0.125)*FilterChunkSec);
            %e2 = min(Time_max, (chunkindexrange(2)+0.125)*FilterChunkSec);
            
            %ti1 = find(Time>=e1,1);
            %ti2 = find(Time<=e2,1,'last');
            
            if HighPassFilter.state
                for ch = size(Signal_postnotch,2):-1:1
                    set(h_bigtext, 'Visible', 'on', 'String', ['Applying high-pass filter (' num2str(ch) ' chans to go)']); drawnow;
                    [Signal3a(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [hpf FilterOrder], 'high', 'butter', FilterReflectMult*FilterOrder/hpf*Fs);
                    while any(FilterInfo.ButterUnstable)
                        if FilterOrder <= 1
                            Signal3a(:,ch) = Signal_postnotch(:,ch);
                            HighPassFilter.state = 0;
                            %set(h_hpf_state, 'String', text_off);
                            set(h_hpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
                            SigChunkRendered(:) = 0;
                            break
                        end
                        FilterOrder = ceil(FilterOrder / 2);
                        [Signal3a(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [hpf FilterOrder], 'high', 'butter', FilterReflectMult*FilterOrder/hpf*Fs);
                    end
                    if ~HighPassFilter.state
                        Signal3a = Signal_postnotch;
                        break
                    end
                end
            else
                Signal3a = Signal_postnotch;
            end
            
            SigBandwidth = Fs/2;
            
            if LowPassFilter.state
                for ch = size(Signal_postnotch,2):-1:1
                    set(h_bigtext, 'Visible', 'on', 'String', ['Applying low-pass filter (' num2str(ch) ' chans to go)']); drawnow;
                    [Signal3b(:,ch), FilterInfo] = freqfilter(Signal3a(:,ch), Fs, [lpf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/lpf*Fs);
                    while any(FilterInfo.ButterUnstable)
                        if FilterOrder <= 1
                            Signal3b(:,ch) = Signal3a(:,ch);
                            LowPassFilter.state = 0;
                            %set(h_lpf_state, 'String', text_off);
                            set(h_lpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
                            SigChunkRendered(:) = 0;
                            break
                        end
                        FilterOrder = ceil(FilterOrder / 2);
                        [Signal3b(:,ch), FilterInfo] = freqfilter(Signal3a(:,ch), Fs, [lpf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/lpf*Fs);
                    end
                    if ~LowPassFilter.state
                        Signal3b = Signal3a;
                        break
                    end
                end
                if LowPassFilter.state
                    SigBandwidth = lpf;
                end
            else
                Signal3b = Signal3a;
            end
            
            if EnvelopeFilter.state
                for ch = size(Signal_postnotch,2):-1:1
                    set(h_bigtext, 'Visible', 'on', 'String', ['Applying envelope filter (' num2str(ch) ' chans to go)']); drawnow;
                    [Signal_postfilter(:,ch), FilterInfo] = freqfilter(Signal3b(:,ch).^2, Fs, [evf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/evf*Fs);
                    while any(FilterInfo.ButterUnstable)
                        if FilterOrder <= 1
                            Signal_postfilter(:,ch) = Signal3b(:,ch);
                            EnvelopeFilter.state = 0;
                            %set(h_evf_state, 'String', text_off);
                            set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
                            SigChunkRendered(:) = 0;
                            break
                        end
                        FilterOrder = ceil(FilterOrder / 2);
                        [Signal_postfilter(:,ch), FilterInfo] = freqfilter(Signal3b(:,ch).^2, Fs, [evf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/evf*Fs);
                    end
                    if ~EnvelopeFilter.state
                        Signal_postfilter = Signal3b;
                        break
                    end
                end
                if EnvelopeFilter.state
                    SigBandwidth = min(SigBandwidth,evf);
                end
            else
                Signal_postfilter = Signal3b;
            end            
            
            
%             set(h_bigtext, 'Visible', 'on', 'String', 'Interpolating to the plot grid...'); drawnow;
%             if BandLimitedInterpolation
%                 BLIM = max(1,floor(Fs/2^nextpow2(SigBandwidth*4)));
%                 Signal_postfilter(ti1:BLIM:ti2,:) = interp1(ETime, Signal3c, Time(ti1:BLIM:ti2));
%             else
%                 Signal_postfilter(ti1:ti2,:) = interp1(ETime, Signal3c, Time(ti1:ti2));
%             end
            
            clear Signal3a Signal3b
            
            SigChunkRendered(chunkindexrange(1):chunkindexrange(2)) = 1;
            FilterBusy = 0;
            enable_movement_switches();
            
            if HighPassFilter.state
                %set(h_hpf_state, 'String', text_on);
                set(h_hpf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            else
                %set(h_hpf_state, 'String', text_off);
                set(h_hpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            end
            if LowPassFilter.state
                %set(h_lpf_state, 'String', text_on);
                set(h_lpf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            else
                %set(h_lpf_state, 'String', text_off);
                set(h_lpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            end
            if EnvelopeFilter.state
                %set(h_evf_state, 'String', text_on);
                set(h_evf_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            else
                %set(h_evf_state, 'String', text_off);
                set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            end
        end
        
        set(h_bigtext, 'Visible', 'on', 'String', 'Rendering...'); drawnow;
        redraw();
        set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
        
        
    end











    function f_xzoomin(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        XRange = XLim(2)-XLim(1);
        [~,in] = min(abs(XRange - PermittedXZoomRanges));
        if in > 1
            XRange = PermittedXZoomRanges(in-1);
        else
            XRange = PermittedXZoomRanges(in);
        end
        XLim(2) = XLim(1) + XRange;
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['Zoomed in to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end

    function f_xzoomout(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        XRange = XLim(2)-XLim(1);
        [~,in] = min(abs(XRange - PermittedXZoomRanges));
        if in < length(PermittedXZoomRanges)
            XRange = PermittedXZoomRanges(in+1);
        else
            XRange = PermittedXZoomRanges(in);
        end
        
        if XRange > Time_max
            XRange = Time_max;
        end
        
        XLim(2) = XLim(1) + XRange;
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['Zoomed out to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end

    function f_yzoomin(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YRange = YRange / 2;
        YLim(1) = YLim(2) - YRange;
        set(axehand, 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['Zoomed in to ' get(h_yzoomlevel, 'String') ' on screen']);
    end

    function f_yzoomout(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        %YCenter = (YLim(1)+YLim(2))/2;
        YRange = YRange * 2;
        YLim(1) = YLim(2) - YRange;
        set(axehand, 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['Zoomed out to ' get(h_yzoomlevel, 'String') ' on screen']);
    end

    function f_panleft(hObject, eventdata) %#ok<*INUSL>
        if FilterBusy || MovementBusy || XLim(1) == Time_min
            return;
        end
        if ~isempty(eventdata) && isnumeric(eventdata) && eventdata > 0 && eventdata < 10
            panfrac = eventdata;
        else
            panfrac = 1;
        end
        XRange = XLim(2)-XLim(1);
        XLim = XLim - XRange*panfrac;
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end

    function f_panright(hObject, eventdata)
        if FilterBusy || MovementBusy || XLim(2) == Time_max
            return;
        end
        if ~isempty(eventdata) && isnumeric(eventdata) && eventdata > 0 && eventdata < 10
            panfrac = eventdata;
        else
            panfrac = 1;
        end        
        XRange = XLim(2)-XLim(1);
        XLim = XLim + XRange*panfrac;
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end

    function f_panup(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YLim = YLim + YRange;
        set(axehand, 'YLim', YLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        yt = get(axehand, 'YTick');
        ind = find(yt >= YLim(1) & yt <= YLim(2));
        ytl = get(axehand, 'YTickLabel');
        if ~isempty(ind)
            set(h_hintbar, 'String', ['Panned up to ' ytl{ind(end)} ' -- ' ytl{ind(1)} ' (' num2str(length(ind)) ' channels)']);
        else
            set(h_hintbar, 'String', '');
        end
    end

    function f_pandown(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YLim = YLim - YRange;
        set(axehand, 'YLim', YLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        yt = get(axehand, 'YTick');
        ind = find(yt >= YLim(1) & yt <= YLim(2));
        ytl = get(axehand, 'YTickLabel');
        if ~isempty(ind)
            set(h_hintbar, 'String', ['Panned down to ' ytl{ind(end)} ' -- ' ytl{ind(1)} ' (' num2str(length(ind)) ' channels)']);
        else
            set(h_hintbar, 'String', '');
        end
    end


    function f_xspan_edit1(hObject, eventdata)
        if FilterBusy || MovementBusy
            set(h_xspan_edit1, 'String', XLim(1));
            return;
        end
        v = str2double(get(h_xspan_edit1, 'String'));
        if isfinite(v) && imag(v) == 0 && v < XLim(2)
            XLim(1) = v;
        end
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end



    function f_xspan_edit2(hObject, eventdata)
        if FilterBusy || MovementBusy
            set(h_xspan_edit2, 'String', XLim(2));
            return;
        end
        v = str2double(get(h_xspan_edit2, 'String'));
        if isfinite(v) && imag(v) == 0 && v > XLim(1)
            XLim(2) = v;
        end
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end


    function f_sepup(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        yl2 = YLim(2)-0.5*chansep;
        yl1 = YLim(1)+0.5*chansep;
        FirstChViewable = ceil(-(yl2+chansep/2)/chansep);
        Nchviewable = ceil((yl2-yl1)/chansep);
        [~,in] = min(abs(chansep - PermittedChanSepRanges));
        if in < length(PermittedChanSepRanges)
            chansep = PermittedChanSepRanges(in+1);
        else
            chansep = PermittedChanSepRanges(in);
        end
        %YLim = [-chansep*Nsch-0.5*chansep, -chansep+0.5*chansep];
        %set(axehand, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames), 'YLim', YLim);
        YLim(2) = 0.5*chansep - FirstChViewable*chansep;
        YLim(1) = YLim(2) - chansep*Nchviewable - 0.5*chansep;
        set(axehand, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        %redraw();
        set(h_hintbar, 'String', ['Decreased sensitivity to ' get(h_sensitivity, 'String') ' between channel name tickmarks']);
    end

    function f_sepdown(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        yl2 = YLim(2)-0.5*chansep;
        yl1 = YLim(1)+0.5*chansep;
        FirstChViewable = ceil(-(yl2+chansep/2)/chansep);
        Nchviewable = ceil((yl2-yl1)/chansep);
        [~,in] = min(abs(chansep - PermittedChanSepRanges));
        if in > 1
            chansep = PermittedChanSepRanges(in-1);
        else
            chansep = PermittedChanSepRanges(in);
        end
        %YLim = [-chansep*Nsch-0.5*chansep, -chansep+0.5*chansep];
        %set(axehand, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames), 'YLim', YLim);
        YLim(2) = 0.5*chansep - FirstChViewable*chansep;
        YLim(1) = YLim(2) - chansep*Nchviewable - 0.5*chansep;
        set(axehand, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        %redraw();
        set(h_hintbar, 'String', ['Increased sensitivity to ' get(h_sensitivity, 'String') ' between channel name tickmarks']);
    end

    function resnap_pan()
        XRange = XLim(2)-XLim(1);
        YRange = YLim(2)-YLim(1);
        if XLim(2) > Time_max
            XLim(2) = Time_max;
            XLim(1) = XLim(2) - XRange;
        end
        if XLim(1) < Time_min
            XLim(1) = Time_min;
            XLim(2) = XLim(1) + XRange;
        end
        if YLim(1) < -chansep*Nsch-chansep/2
            YLim(1) = -chansep*Nsch-chansep/2;
            YLim(2) = YLim(1) + YRange;
        end
        if YLim(2) > -chansep/2
            YLim(2) = -chansep/2;
            YLim(1) = YLim(2) - YRange;
        end
        
        
        XLim = round_xlim(XLim, XRange);
        XLim(2) = XLim(1) + XRange;
        set(axehand, 'XLim', XLim, 'YLim', YLim);
        set(h_xspan_edit1, 'String', num2str(XLim(1)));
        set(h_xspan_edit2, 'String', num2str(XLim(2)));
        filter_render();
        redraw();
        
    end


    function XLim = round_xlim(XLim, XRange)
        if XRange >= 1
            XLim(1) = floor(XLim(1)*FineSnapScale)/FineSnapScale;
        elseif XRange >= 0.1
            XLim(1) = floor(XLim(1)*10*FineSnapScale)/10/FineSnapScale;
        elseif XRange >= 0.01
            XLim(1) = floor(XLim(1)*100*FineSnapScale)/100/FineSnapScale;
        elseif XRange >= 0.001
            XLim(1) = floor(XLim(1)*1000*FineSnapScale)/1000/FineSnapScale;
        elseif XRange >= 0.0001
            XLim(1) = floor(XLim(1)*10000*FineSnapScale)/10000/FineSnapScale;
        elseif XRange >= 0.00001
            XLim(1) = floor(XLim(1)*100000*FineSnapScale)/100000/FineSnapScale;
        elseif XRange >= 0.000001
            XLim(1) = floor(XLim(1)*1000000*FineSnapScale)/1000000/FineSnapScale;
        end
    end

    function resnap_zoom()
        XRange = XLim(2)-XLim(1);
        YRange = YLim(2)-YLim(1);
        
        if YRange > chansep/2 - (-chansep*Nsch-chansep/2)
            YRange = chansep/2 - (-chansep*Nsch-chansep/2);
            YLim(1) = YLim(2) - YRange;
        end
        if YRange < chansep
            YRange = chansep;
            YLim(1) = YLim(2) - YRange;
        end
        
        YRange = ceil(YRange/chansep)*chansep;
        YLim(1) = YLim(2) - YRange;

        XLim = round_xlim(XLim, XRange);
        XLim(2) = XLim(1) + XRange;
        
        set(axehand, 'XLim', XLim, 'YLim', YLim);
        set(h_xspan_edit1, 'String', num2str(XLim(1)));
        set(h_xspan_edit2, 'String', num2str(XLim(2)));
        filter_render();
        redraw();

    end

    function redraw()
        %fprintf('0001 axehand xlim %g %g\n', get(axehand, 'XLim'));
        t1 = find(Time<=XLim(1),1,'last');
        t2 = find(Time>=XLim(2),1,'first');
        if isempty(t1)
            t1 = 1;
        end
        if isempty(t2)
            t2 = Ntp;
        end
        
        Nsch = length(selchan);
        tdrawupdate = tic;
        for ch = randperm(Nsch)
            
            if -chansep/2-chansep*ch < YLim(1) || chansep/2-chansep*ch > YLim(2)
                % out of plotting range
                %fprintf('%i is out of range\n', ch);
                set(plothand(ch), 'Visible', 'off');
                continue;
            end
            
            if ~PlotHold
                
                if ~isempty(EventTimePoints)
                    % NaN around blankaround_stitch_samples
                    l1l = EventTimePoints(:,1) >= t1 & EventTimePoints(:,1) <= t2;
                    l2l = EventTimePoints(:,2) >= t1 & EventTimePoints(:,2) <= t2;
                    lcomlist = find(l1l | l2l);
                else
                    lcomlist = [];
                end
                if ~isempty(lcomlist)
                    tmp = Signal_postfilter(:,selchan(ch));
                    for i2 = 1:length(lcomlist)
                        i = lcomlist(i2);
                        if blankaround_stitch_samples > 0 && blankaround_stitch_samples <= EventTimePoints(i,2)-EventTimePoints(i,1)+1
                            tmp(round(EventTimePoints(i,1) + [0:blankaround_stitch_samples-1]),:) = NaN;
                            tmp(round(EventTimePoints(i,2) - [blankaround_stitch_samples-1:-1:0]),:) = NaN;
                        end
                    end
                    Signal4 = tmp(t1:BLIM:t2,:);
                else
                    Signal4 = Signal_postfilter(t1:BLIM:t2,selchan(ch));
                end
                
                Time4 = Time(t1:BLIM:t2);
                if ScreenLimitedDownsampling && length(Time4) > 2*SLD_H
                    tdiff = Time(t2) - Time(t1);
                    Fs_pref = 2^nextpow2(SLD_H / tdiff);
                    if Fs_pref < Fs/BLIM
                        [Signal4, Time4] = downsamplecustom(Signal4, Time4, SLD_H);
                    end
                end
                
                YDATA = Signal4;
                
                if ZscoreFilter.state
                    YDATA = nanzscore(YDATA)*ZscoreFilter.multiplier;
                else
                    YDATA = YDATA - nanmedian(YDATA);
                end
                
                set(plothand(ch), 'XData', Time4, 'YData', YDATA - chansep*ch, 'Visible', 'on');
                set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:));
                setappdata(plothand(ch), 'chanind', selchan(ch));
                setappdata(plothand(ch), 'channame', ChanNames{selchan(ch)});

                update_psd();
                
            else
                set(plothand(ch), 'Visible', 'off');
            end

            
            if toc(tdrawupdate) > 0.2
                set(axehand, 'YColor', BusyYColor);
                drawnow
                tdrawupdate = tic;
            end
            
        end
        for ch = length(plothand):-1:Nsch+1
            set(plothand(ch), 'Visible', 'off');
        end
        
        if PlotHold
            set(axehand, 'YColor', InactiveYColor);
        else
            set(axehand, 'YColor', DefaultYColor);
        end
        
        [~,in] = min(abs((Time(t2)-Time(t1))./XTickSpacings - 20));
        tm1 = floor(Time(t1)/XTickSpacings(in))*XTickSpacings(in);
        tm2 = ceil(Time(t2)/XTickSpacings(in))*XTickSpacings(in);
        tms = tm1:XTickSpacings(in):tm2;
        tmslabel = cell(1,length(tms));
        for i = 1:length(tms)
            tmslabel{i} = sprintf('%g%s', tms(i)/XTickSpacingsAndUnits{in,2}, XTickSpacingsAndUnits{in,3});
        end
        set(axehand, 'XTick', tms, 'XTickLabel', tmslabel);
        
        if EventEnable
            %YPos = [mean(YLim)+diff(YLim)/8     mean(YLim)    mean(YLim)-diff(YLim)/8];
            YPos = sort(mean(YLim)+diff(YLim)/20*([-8:2:8]), 'descend');
            NYPos = length(YPos);
            for i = size(EventTimeStamps,1):-1:1
                tmp = get(eventtexthand(i), 'Position');
                tmp(2) = YPos(mod(i-1,NYPos)+1);
                set(eventtexthand(i), 'Position', tmp);
            end
            clear tmp
        end

        
        XRange = XLim(2)-XLim(1);
        if XRange < 1e-6
            set(h_xzoomlevel, 'String', sprintf('%.3g ns', XRange*1000000000));
        elseif XRange < 1e-3
            set(h_xzoomlevel, 'String', sprintf('%.3g µs', XRange*1000000));
        elseif XRange < 1
            set(h_xzoomlevel, 'String', sprintf('%.3g ms', XRange*1000));
        elseif XRange < 60
            set(h_xzoomlevel, 'String', sprintf('%.3g s', XRange));
        elseif XRange < 3600
            set(h_xzoomlevel, 'String', sprintf('%.3g min', XRange/60));
        else
            set(h_xzoomlevel, 'String', sprintf('%.3g hr', XRange/3600));
        end
        
        YRange = YLim(2)-YLim(1);
        set(h_yzoomlevel, 'String', sprintf('%g ch', YRange/chansep));
        
        tmp_sq = '';
        if EnvelopeFilter.state
            tmp_sq = '²';
        end
        if ZscoreFilter.state
            set(h_sensitivity, 'String', sprintf('%g sdev', chansep/ZscoreFilter.multiplier));
        elseif chansep < 1
            set(h_sensitivity, 'String', sprintf('%g nV%s', chansep*1000, tmp_sq));
        elseif chansep < 1000
            set(h_sensitivity, 'String', sprintf('%g µV%s', chansep, tmp_sq));
        elseif chansep < 1e6
            set(h_sensitivity, 'String', sprintf('%g mV%s', chansep/1000, tmp_sq));
        else
            set(h_sensitivity, 'String', sprintf('%g V%s', chansep/1000000, tmp_sq));
        end
        
        if EventEnable
            l = XLim(1) <= EventTimes & EventTimes <= XLim(2);
            set(eventtexthand(l), 'Visible', 'on');
            set(eventtexthand(~l), 'Visible', 'off');
        end
        
        
        update_cursorline();
        
    end


    function f_chansel_auto(hObject, eventdata)
        try
            cids = GetChannelsWithSimilarStdev(Signal_postfilter);
            set(h_chansel_list, 'Value', cids);
            f_chansel_confirm([], []);
        catch
        end
    end

    function f_cursor_enable(hObject, eventdata)
        if CursorEnable
            CursorEnable = 0;
            %set(h_cursor_state, 'String', text_off);
            set(h_cursor_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            update_cursorline();
            set(h_hintbar, 'String', '');
        else
            CursorEnable = 1;
            %set(h_cursor_state, 'String', text_on);
            set(h_cursor_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            set(h_hintbar, 'String', 'Click on signal to get time/value. Hold Shift to get local max. Hold Ctrl to get local min.');
        end
    end

    function f_chansel_list(hObject, eventdata)
        set(h_hintbar, 'String', 'Hold Ctrl to select/deselect multiple. Hold Shift to select contiguous range. Click "Plot" to finalize selection.');
    end


    function f_chansel_confirm(hObject, eventdata)
        selchan = get(h_chansel_list, 'Value');
        Nsch = length(selchan);
        set(axehand, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)));
        redraw();
        %Automatically zoom out if less-than-full channels fill the screen
        if Nsch < (YLim(2)-YLim(1))/chansep
            f_yzoomout(hObject, eventdata);
        end
    end

    function f_chansel_reset(hObject, eventdata)
        set(h_chansel_list, 'Value', selchan);
    end

    function f_psd_plot(hObject, eventdata)
        if ~ishandle(viewhand_psd)
            viewhand_psd = figure;
            viewhand_psd_axe = axes('parent',viewhand_psd);
            set(viewhand_psd, 'KeyPressFcn', @f_fig_keypress);
            Signal_psd_source = Signal_psd_source*0;
            update_psd();
            set(viewhand_psd, 'MenuBar', 'none', 'Toolbar', 'figure', 'HandleVisibility', 'callback');
        else
            figure(viewhand_psd);
        end
        set(h_hintbar, 'String', 'PSD updated');
    end


    function update_cursorline(sel_type)
        
        if ~exist('sel_type','var') || isempty(sel_type)
            sel_type = selected_cursortype;
        end
        
        set([cursorlinehand plottext1hand plottext2hand plottext3hand plotpip1hand plotpip2hand plotpip3hand], 'Visible', 'off');
        
        if ~CursorEnable
            return
        end
        
        
        if sel_type == 0
            return
        end
        
        
        selected_cursortype = sel_type;
        
        
        yl = get(axehand, 'YLim');
        set(cursorlinehand, 'YData', yl, 'XData', selected_timepoint*[1 1], 'Visible', 'on');
        
        ypoint = [];
        ylocmin = [];
        xlocmin = [];
        ylocmax = [];
        xlocmax = [];
        xindmin = [];
        xindmax = [];

        t1 = find(Time<=XLim(1),1,'last');
        t2 = find(Time>=XLim(2),1,'first');
        if isempty(t1)
            t1 = 1;
        end
        if isempty(t2)
            t2 = Ntp;
        end

        for ch = Nsch:-1:1
            if get(plothand(ch), 'Visible')
                % y values of the exact x point
                
                %YDATA = get(plothand(ch), 'YData');
                YDATA_fullres = Signal_postfilter(t1:BLIM:t2,selchan(ch));
                if ZscoreFilter.state
                    YDATA_fullres = nanzscore(YDATA_fullres)*ZscoreFilter.multiplier - chansep*ch;
                else
                    YDATA_fullres = YDATA_fullres - nanmedian(YDATA_fullres) - chansep*ch;
                end
                
                ypoint(ch) = interp1(Time, Signal_postfilter(:,selchan(ch)), selected_timepoint);
                
                % Locally search for min and max in each channel
                lookcenter = find(Time >= selected_timepoint, 1, 'first');
                loclookcenter = find(Time(t1:BLIM:t2) >= selected_timepoint, 1, 'first');
                
                if isempty(loclookcenter) || loclookcenter == 1
                    return
                end
                
                lookradius = ceil((t2 - t1) * 0.005);
                
                tmpunit = 'µV';
                if ZscoreFilter.state
                    tmpunit = 'µV (without zscore)';
                elseif EnvelopeFilter.state
                    tmpunit = 'µV²';
                end
                
                xr = diff(get(axehand, 'XLim'));
                
                if sel_type == 1
                    ypointdisplayed(ch) = YDATA_fullres(loclookcenter);
                    set(plotpip1hand(ch), 'XData', selected_timepoint, 'YData', ypointdisplayed(ch), 'Visible', 'on');
                    set(plottext1hand(ch), 'String', sprintf('%s Here\n%.6gs\n%.6g%s',ChanNames{selchan(ch)}, selected_timepoint,ypoint(ch),tmpunit), 'Position', [selected_timepoint + xr/100, ypointdisplayed(ch), 0], 'Visible', 'on');
                elseif sel_type == 2
                    if lookcenter - lookradius < 1 || lookcenter + lookradius > Ntp
                        return
                    end
                    tmp = Signal_postfilter(lookcenter-lookradius:lookcenter+lookradius,selchan(ch));
                    [ylocmin(ch), ind] = min(tmp);
                    xindmin(ch) = ind - 1 + loclookcenter - lookradius;
                    xlocmin(ch) = Time(t1-1+xindmin(ch));
                    if xindmin(ch) < 1 || xindmin(ch) > length(YDATA_fullres)
                        continue
                    end
                    set(plotpip2hand(ch), 'XData', xlocmin(ch), 'YData', YDATA_fullres(xindmin(ch)), 'Visible', 'on');
                    set(plottext2hand(ch), 'String', sprintf('%s LMin\n%.6gs\n%.6g%s',ChanNames{selchan(ch)}, xlocmin(ch), ylocmin(ch),tmpunit), 'Position', [xlocmin(ch) + xr/100, YDATA_fullres(xindmin(ch)), 0], 'Visible', 'on');
                elseif sel_type == 3
                    if lookcenter - lookradius < 1 || lookcenter + lookradius > Ntp
                        return
                    end
                    tmp = Signal_postfilter(lookcenter-lookradius:lookcenter+lookradius,selchan(ch));
                    [ylocmax(ch), ind] = max(tmp);
                    xindmax(ch) = ind - 1 + loclookcenter - lookradius;
                    xlocmax(ch) = Time(t1-1+xindmax(ch));
                    if xindmax(ch) < 1 || xindmax(ch) > length(YDATA_fullres)
                        continue
                    end
                    set(plotpip3hand(ch), 'XData', xlocmax(ch), 'YData', YDATA_fullres(xindmax(ch)), 'Visible', 'on');
                    set(plottext3hand(ch), 'String', sprintf('%s LMax\n%.6gs\n%.6g%s',ChanNames{selchan(ch)}, xlocmax(ch), ylocmax(ch),tmpunit), 'Position', [xlocmax(ch) + xr/100, YDATA_fullres(xindmax(ch)), 0], 'Visible', 'on');
                end
                
            end
        end
        
        
        
    end




    function update_psd()
        if ishandle(viewhand_psd)
            if ~isempty(selected_plothand) && ishandle(selected_plothand)
                channame = getappdata(selected_plothand, 'channame');
                chanind = getappdata(selected_plothand, 'chanind');
                tmp = Signal_postfilter(t1:t2,chanind);
                tmp = tmp(isfinite(tmp));
                if size(Signal_psd_source,1) == size(tmp,1) && size(Signal_psd_source,2) == size(tmp,2) && norm(Signal_psd_source - tmp) == 0
                    % These are the same signal.
                    return
                end
                tmp1 = fighand;
                set(0, 'CurrentFigure', viewhand_psd);
                
                set(viewhand_psd_axe, 'FontUnits', 'pixel', 'FontSize', 16);
                Signal_psd_source = tmp;
                [pxx, fxx] = pwelch(Signal_psd_source, [], [], pwelch_nfft, Fs);
                plot(fxx, 10*log10(pxx));
                xlabel('Frequency (Hz)');
                ylabel('PSD (dB/Hz)');
                fc1 = HighPassFilter.cutoff;
                fc2 = LowPassFilter.cutoff;
                if ~HighPassFilter.state
                    fc1 = 0;
                end
                if ~LowPassFilter.state
                    fc2 = Fs/2;
                end
                

                % Applicable filters are: ICA, Notch, HPF, LPF, Env
                tmp_af = {'ICA', 'Notch', 'HPF', 'LPF', 'Env'};
                tmp_ftr = false(1,5);
                if length(selica) ~= size(ica_A,2)
                    tmp_ftr(1) = 1;
                end
                if NotchFilter.state
                    tmp_ftr(2) = 1;
                end
                if HighPassFilter.state
                    tmp_ftr(3) = 1;
                end
                if LowPassFilter.state
                    tmp_ftr(4) = 1;
                end
                if EnvelopeFilter.state
                    tmp_ftr(5) = 1;
                end
                
                if any(tmp_ftr)
                    tmp_filttext = ['(after ' cell_to_string(tmp_af(tmp_ftr), ', ') ')'];
                else
                    tmp_filttext = ['(no filter applied)'];
                end
                
                try  %#ok<*TRYNC>
                    set(viewhand_psd_axe, 'XLim', [fc1 fc2]);
                end
                title(['Welch PSD in ' channame ', ' num2str(Time(t1))  '-' num2str(Time(t2)) ' s']);
                set(viewhand_psd, 'Name', ['PSD in ' channame ' from ' num2str(Time(t1)) ' to ' num2str(Time(t2)) ' s ' tmp_filttext]);
                set(0, 'CurrentFigure', tmp1);
            else
                set(viewhand_psd, 'Name', 'Select a channel first by clicking on its signal.');
            end
        end
    end


    function StitchSignalCell()
        if SignalIsStitched
            wh = waitbar(0, 'Stitching segments together...', 'Name', 'SignalViewer');
            tmp = Signal;
            averagesegmentduration = mean(diff(EventTimePoints,[],2)+1) / Fs;
            if stitch_mult > 0
                Fcut_minfreq = 2*stitch_mult/averagesegmentduration;
                if Fcut_minfreq < Fs / 10000
                    Fcut_minfreq = Fs / 10000;
                end
            end
            imax = size(EventTimePoints,1);
            for i = 1:imax
                % Do a detrend
                tmp(EventTimePoints(i,1):EventTimePoints(i,2),:) = detrend(Signal(EventTimePoints(i,1):EventTimePoints(i,2),:), 'linear');
                
                % Do a high-pass if stitch_mult > 0
                if stitch_mult > 0
                    reflect_len = averagesegmentduration * Fs / 2;
                    Fcut_stitch = Fcut_minfreq;
                    FO_stitch = 2;
                    tmp(EventTimePoints(i,1):EventTimePoints(i,2),:) = freqfilter(tmp(EventTimePoints(i,1):EventTimePoints(i,2),:), Fs, [Fcut_stitch, FO_stitch], 'high', 'butter', reflect_len);
                end
                
                if ishandle(wh)
                    waitbar(i/imax, wh);
                end
            end
            Signal = tmp;
            if ishandle(wh)
                delete(wh);
            end
            clear tmp
        end
    end


    function f_icasel_confirm(hObject, eventdata)
        if FilterBusy
            return;
        end
        
        if isempty(ica_W)
            % Use default ICA arguments
            FilterBusy = 1;
            set(h_icasel_confirm, 'Enable', 'off', 'String', 'Wait');
            set(h_bigtext, 'Visible', 'on', 'String', 'ICA: Calculating ICs...'); drawnow;
            try
                [~, ica_A, ica_W] = fastica(Signal.', 'stabilization', 'on', 'maxNumIterations', 200);
            catch exception
                ica_A = [];
                ica_W = [];
                set(h_bigtext, 'Visible', 'on', 'String', 'ICA: Failure. See command window for details.'); drawnow;
                disp('ICA failure:');
                disp(exception.message);
                
                pause(5.0);
                set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
                FilterBusy = 0;
                set(h_icasel_confirm, 'Enable', 'off', 'String', 'No ICA');
                set(h_icasel_reset, 'Enable', 'off');
                set(h_icasel_view_sources, 'Enable', 'off');
                set(h_icasel_view_mixmat, 'Enable', 'off');
                set(h_icasel_view_sepmat, 'Enable', 'off');
                set(h_icasel_view_export, 'Enable', 'off');
                return;
            end
            set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
            FilterBusy = 0;
        elseif size(Signal,2) ~= size(ica_W,2) || size(Signal,2) ~= size(ica_A,1)
            % ICA parameters are wrong. Disable it
            set(h_icasel_confirm, 'Enable', 'off', 'String', 'No ICA');
            set(h_icasel_reset, 'Enable', 'off');
            set(h_icasel_view_sources, 'Enable', 'off');
            set(h_icasel_view_mixmat, 'Enable', 'off');
            set(h_icasel_view_sepmat, 'Enable', 'off');
            set(h_icasel_view_export, 'Enable', 'off');
            return;
        end
        
        if ~ICA_Initialized
            FilterBusy = 1;
            set(h_icasel_confirm, 'Enable', 'off', 'String', 'Wait');
            set(h_bigtext, 'Visible', 'on', 'String', 'ICA: Separating sources...'); drawnow;
            ica_sig = ica_W * Signal.';
            mica = size(ica_A,2);
            selica = mica;
            icachans = string_to_cell(num2str(1:mica, 'ic%i '), ' ');
            set(h_icasel_list, 'String', icachans);
            set(h_icasel_list, 'Value', 1:mica);
            set(h_icasel_confirm, 'String', 'Mix', 'Enable', 'on');
            set(h_icasel_reset, 'Enable', 'on');
            set(h_icasel_view_sources, 'Enable', 'on');
            set(h_icasel_view_mixmat, 'Enable', 'on');
            set(h_icasel_view_sepmat, 'Enable', 'on');
            set(h_icasel_view_export, 'Enable', 'on');
            ICA_Initialized = 1;
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
            f_icasel_view_sources([], []);
        else
            tmp = get(h_icasel_list, 'Value');
            if length(tmp) == length(selica) && min(tmp == selica) == 1
                % No change
                return;
            else
                selica = tmp;
            end
            if length(selica) == size(ica_A,2)
                % All ICs selected
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Mixing new ICs...'); drawnow;
                Signal_postica = Signal;
                FilterBusy = 0;
                set(h_bigtext, 'Visible', 'off', 'String', '');
            else
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Mixing new ICs...'); drawnow;
                tmp = ica_A;
                tmp(:,setdiff(1:size(ica_A,2),selica)) = 0;
                Signal_postica = (tmp*ica_sig).';
                FilterBusy = 0;
                set(h_bigtext, 'Visible', 'off', 'String', '');
            end
            notch_update();
            filter_update();
        end
    end

    function f_icasel_reset(hObject, eventdata)
        if FilterBusy
            return;
        end
        set(h_icasel_list, 'Value', selica);
    end

    function f_icasel_view_sources(hObject, eventdata)
        if FilterBusy
            return;
        end
        if ~ishandle(viewhand_ica_sig)
            npad2 = floor(log10(size(ica_sig,1)))+1;
            ic_names = string_to_cell(num2str(1:size(ica_sig,1),['ic%0' num2str(npad2) 'i,']),',');
            viewhand_ica_sig = signalviewer(ica_sig.', Fs, ic_names, [], [], [], 'This figure window displays the ICA sources. Close this window to return to the original time signals. Exclude and remix ICs in the original time signals, not in here.');
        else
            figure(viewhand_ica_sig);
        end
    end

    function f_icasel_view_mixmat(hObject, eventdata)
        if FilterBusy
            return;
        end
        if ~ishandle(viewhand_ica_A)
            viewhand_ica_A = figure;
            set(viewhand_ica_A, 'Name', 'ICA Mixing Matrix (A)');
            tmp = imagesc(ica_A);
            viewhand_ica_A_ax = get(tmp, 'Parent');
            tmp = get(viewhand_ica_A_ax, 'CLim');
            set(viewhand_ica_A_ax, 'CLim', max(abs(tmp))*[-1 1]);
            set(viewhand_ica_A_ax, 'FontUnits', 'normalized', 'FontSize', 8);
            xlabel('Sources');
            ylabel('Channels');
            set(viewhand_ica_A_ax, 'YTick', 1:size(ica_A,1));
            set(viewhand_ica_A_ax, 'YTickLabel', ChanNames);
            set(viewhand_ica_A_ax, 'XTick', 1:size(ica_A,2));
        else
            figure(viewhand_ica_A);
        end
    end

    function f_icasel_view_sepmat(hObject, eventdata)
        if FilterBusy
            return;
        end
        if ~ishandle(viewhand_ica_W)
            viewhand_ica_W = figure;
            set(viewhand_ica_W, 'Name', 'ICA Separating Matrix (W)');
            tmp = imagesc(ica_W);
            viewhand_ica_W_ax = get(tmp, 'Parent');
            tmp = get(viewhand_ica_W_ax, 'CLim');
            set(viewhand_ica_W_ax, 'CLim', max(abs(tmp))*[-1 1]);
            set(viewhand_ica_W_ax, 'FontUnits', 'normalized', 'FontSize', 8);
            xlabel('Channels');
            ylabel('Sources');
            set(viewhand_ica_W_ax, 'YTick', 1:size(ica_W,1));
            set(viewhand_ica_W_ax, 'XTick', 1:size(ica_W,2));
            set(viewhand_ica_W_ax, 'XTickLabel', ChanNames);
        else
            figure(viewhand_ica_W);
        end
    end


    function f_icasel_view_export(hObject, eventdata)
        if FilterBusy
            return;
        end
        assignin('caller', 'signalviewer_ica_A', ica_A);
        assignin('caller', 'signalviewer_ica_W', ica_W);
        assignin('caller', 'signalviewer_ica_sig', ica_sig);
        assignin('caller', 'signalviewer_selica', selica);
        assignin('caller', 'signalviewer_ica_export_timestamp', now);
        set(h_hintbar, 'String', 'ICA variables exported to caller workspace.');
    end


    function f_signal_export(hObject, eventdata)
        assignin('caller', 'signalviewer_Signal', Signal);
        assignin('caller', 'signalviewer_Signal_postica', Signal_postica);
        assignin('caller', 'signalviewer_Signal_postnotch', Signal_postnotch);
        assignin('caller', 'signalviewer_Signal_postfilter', Signal_postfilter);
        assignin('caller', 'signalviewer_selchan', selchan);
        assignin('caller', 'signalviewer_signal_export_timestamp', now);
        set(h_hintbar, 'String', 'Signal variables exported to caller workspace.');
    end


    function f_axesfont_inc(hObject, eventdata)
        AxesFontSize = AxesFontSize + .001;
        set(axehand, 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
        set([plottext1hand plottext2hand plottext3hand], 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
    end

    function f_axesfont_dec(hObject, eventdata)
        if AxesFontSize > .002
            AxesFontSize = AxesFontSize - .001;
        end
        set(axehand, 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
        set([plottext1hand plottext2hand plottext3hand], 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
    end

    function f_windowhsize_inc(hObject, eventdata)
        %AxesPosition = [0.0500    0.0600    0.86    0.93];
        AxesPosition = AxesPosition + [-0.01 0 0.01 0];
        set(axehand, 'Position', AxesPosition);
    end

    function f_windowhsize_dec(hObject, eventdata)
        %AxesPosition = [0.0500    0.0600    0.86    0.93];
        if AxesPosition(3) > 0
            AxesPosition = AxesPosition - [-0.01 0 0.01 0];
        end
        set(axehand, 'Position', AxesPosition);
    end

    function f_plothand_buttondown(hObject, eventdata)
        selected_plothand = hObject;
        update_psd();
        
        acp = get(axehand, 'CurrentPoint');
        selected_timepoint = acp(1);
        if selected_timepoint == previously_selected_timepoint
            previously_selected_timepoint = selected_timepoint - 2;
            update_cursorline(0);
            return;
        else
            previously_selected_timepoint = selected_timepoint;
        end
        switch get(fighand,'SelectionType')
            case 'normal'
                update_cursorline(1);
            case 'alt'
                % Local min (hold Ctrl)
                update_cursorline(2);
            case 'extend'
                % Local max (hold Shift)
                update_cursorline(3);
        end
    end

    function f_main_close(hObject, eventdata, handles)
        if ishandle(viewhand_ica_sig)
            close(viewhand_ica_sig);
        end
        if ishandle(viewhand_ica_A)
            close(viewhand_ica_A);
        end
        if ishandle(viewhand_ica_W)
            close(viewhand_ica_W);
        end
        if ishandle(viewhand_psd)
            close(viewhand_psd);
        end
        delete(hObject);
    end

    function f_autofit(hObject, eventdata)
        autofit();
    end

    function autofit()
        Fs_int = round(Fs);
        tmpr = 1+Fs_int:size(Signal_postfilter,1)-Fs_int+1;
        if isempty(tmpr)
            tmpr = 1:size(Signal_postfilter,1);
        end
        if ZscoreFilter.state
            %tmp = max(max(abs(zscore(Signal_postfilter(tmpr,:))*ZscoreFilter.multiplier)));
            tmp = quantile(  max(nanzscore(Signal_postfilter(tmpr,:))*ZscoreFilter.multiplier) - min(nanzscore(Signal_postfilter(tmpr,:))*ZscoreFilter.multiplier), 0.5 );
        else
            %tmp = max(max(abs(Signal_postfilter(tmpr,:))));
            tmp = quantile(max(Signal_postfilter(tmpr,:)) - min(Signal_postfilter(tmpr,:)),0.5);
        end
        refit(tmp);
    end

    function refit(new_chansep)
        while chansep > new_chansep
            oldchansep = chansep;
            f_sepdown([], []);
            if oldchansep == chansep
                break
            end
        end
        while chansep < new_chansep
            oldchansep = chansep;
            f_sepup([], []);
            if oldchansep == chansep
                break
            end
        end
    end


end


function C = lightercolor(C)
C = 1-((1 - C)/4);
end

function B = rearrange_top_bottom(A)
n = size(A,1);
B = A;
for i = 2:n
    if mod(i,2) == 0
        B(i,:) = A(end-(i/2-1),:);
    else
        B(i,:) = A((i+1)/2,:);
    end
end
end

% function B = rearrange_top_mid(A)
% B = rearrange_top_bottom(A);
% B(2:2:end,:) = flipud(B(2:2:end,:));
% end

function [A2, T2] = downsamplecustom(A1, T1, MaxPoints)
n = length(T1);
m = ceil(n/MaxPoints);
q = floor(size(A1,1)/m);
A2 = A1(1:m:q*m,:);
for i = 2:m
    A2 = A2 + A1(i:m:q*m,:);
end
A2 = A2/m;
T2 = T1(1:m:q*m);
end


function c = string_to_cell(s,d)
% Convert a delimited string to a cell array
% E.g., input is    "blah 1" "blah 2", delimiter is ",
%           output:    {'blah 1', 'blah 2'}

% Copyright 2001-2004 The MathWorks, Inc.
% $Revision: 1.1.8.1 $ $Date: 2004/07/21 06:23:56 $

c = {};
while containsValidString(s),
    [s1, s] = strtok(s, d); %#ok<*STTOK>
    if containsValidString(s1)
        c = {c{:} s1}; %#ok<*CCAT>
    end
end

% ---------------------------------------
    function ok = containsValidString(s)
        % Decide whether there is still valid data in s.
        % I.e., if s only contains separators, quotes, spaces,
        % newlines, etc (in any combination), then it
        % is not valid.
        % This is to be decided in the context of
        % valid filenames, valid code symbols, etc.
        
        goodChars = [ ...
            'abcdefghijklmnopqrstuvwxyz' ...
            'ABCDEFGHIJKLMNOPQRSTUVWXYZ' ...
            '1234567890' ...
            '_~-.!#$%'];
        % !"#$%&'()*+,-./0123456789:;<=>?@
        % [\]^_`
        s2 = strtok(s, goodChars);
        % If s2 does not contain any of these characters,
        % s and s2 will be equal.
        ok = ~isequal(s2, s);
    end
end



function y = fix_nans_for_filtering(y)
% 20180514 Remove NaNs before processing
ny = isnan(y);

if nnz(ny) > 0
    y = bridge_nans(y, 'linear', '');
    
    % If any NaNs still exist at either ends, reflect each channel's data
    for ch = 1:size(y,2)
        
        if ~any(isfinite(y(:,ch)))
            continue
        end
        
        
        nyc = isnan(y(:,ch));
        if any(nyc)
            cgnan = get_contig_groups(find(~isfinite(y(:,ch))));
            cgfin = get_contig_groups(find(isfinite(y(:,ch))));
            cgfin = cgfin{1};
            tmp = y(cgfin(1):cgfin(end),ch);
            
            if nyc(1)
                % NaN at start
                cgnan1 = cgnan{1};
                nn = cgnan1(end)-cgnan1(1)+1;
                nm = ceil(nn / length(tmp));
                ytmp = [];
                for j = 1:nm
                    if mod(j,2) == 1
                        ytmp = [-flipud(tmp(2:end)); ytmp];
                    else
                        ytmp = [tmp(1:end-1); ytmp];
                    end
                end
                y(cgnan1(1):cgnan1(end),ch) = ytmp(end-nn+1:end);
            end
            
            if nyc(end)
                % NaN at end
                cgnan2 = cgnan{end};
                nn = cgnan2(end)-cgnan2(1)+1;
                nm = ceil(nn / length(tmp));
                ytmp = [];
                for j = 1:nm
                    if mod(j,2) == 1
                        ytmp = [ytmp; -flipud(tmp(1:end-1))];
                    else
                        ytmp = [ytmp; tmp(2:end)];
                    end
                end
                y(cgnan2(1):cgnan2(end),ch) = ytmp(1:nn);
            end
            
            
        end
    end
end

end



function chanidx = GetChannelsWithSimilarStdev(signal)
s = std(signal);
[h, c] = hist(s);
[~, i] = max(h);
l1 = s < c(i) + (c(2)-c(1))/2;
l2 = s > c(i) - (c(2)-c(1))/2;
chanidx = find(l1 & l2);
end

function utc_time = unixtimemillis ()
utc_time = java.lang.System.currentTimeMillis;
end

function z = nanzscore (x)
% https://www.mathworks.com/matlabcentral/answers/249566-zscore-a-matrix-with-nan
mu = mean(x,'omitnan');
sigma = std(x, 'omitnan');
z = bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma);
end


