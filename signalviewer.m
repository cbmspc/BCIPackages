% This program expects Signal to be in microvolts and SampleRate to be in
% samples per second
%
% Optional input arguments:
%   EventTimeStamps should be a Nx2 cell {1.2345, 'Abc'; 4.5, 'Def'} with first column in seconds, and 2nd column a string
%   ica_W = ICA separating matrix with orientation (Nsource x Nchan), i.e. ica_W * Signal.' = Source.'
%   ica_A = ICA mixing matrix with orientation (Nchan x Nsource), i.e. ica_A * Source.' = Signal.'
%   Example: [ica_sig, ica_A, ica_W] = fastica(Signal.', 'stabilization', 'on', 'maxNumIterations', 200);
%   If ica matrices are not specified, FastICA is used

function fighand = signalviewer(Signal, SampleRate, ChanNames, EventTimeStamps, ica_W, ica_A)
%t_program_start = tic;
if ~exist('ChanNames', 'var') || isempty(ChanNames)
    nchan = size(Signal,2);
    npad = floor(log10(nchan))+1;
    ChanNames = string_to_cell(num2str(1:nchan,['ch%0' num2str(npad) 'i,']),',');
end
ChanNames = ChanNames(:).';
if size(Signal,1) == length(ChanNames) && size(Signal,2) ~= length(ChanNames)
    Signal = Signal.';
end
Fs = SampleRate;

if ~exist('ica_W', 'var') || ~exist('ica_A', 'var')
    ica_W = [];
    ica_A = [];
end

viewhand_ica_sig = ceil(rand*1000000000);
viewhand_ica_A = ceil(rand*1000000000);
viewhand_ica_W = ceil(rand*1000000000);
viewhand_psd = ceil(rand*1000000000);
selected_plothand = -1;
FineSnapScale = 10;

% Remove channels that are completely flat
n = size(Signal,1);
if n > 491520
    % Sample only these many points in 4 areas
    chnc = true(1,size(Signal,2));
    b = floor(n/4)*(0:3)' * [1 1] + ones(4,1)*[1 122880];
    for ch = 1:size(Signal,2) %#ok<*FXUP>
        if chnc(ch)
            for i = 1:size(b,1)
                chnc(ch) = chnc(ch) & nanstd(Signal(b(i,1):b(i,2),ch),[],1) == 0;
                if ~chnc(ch)
                    break
                end
            end
        end
    end
else
    chnc = nanstd(Signal,[],1) == 0;
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
Signal_psd_source = Signal(:,1);

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
PermittedXZoomRanges = PermittedXZoomRanges(PermittedXZoomRanges > 16/Fs);
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
EventFontSize = 20;

AxesFontName = 'Consolas';
AxesFontSize = 10;

AxesPosition = [0.0500    0.0600    0.86    0.93];

chansep = 100;

screensize = get(groot,'Screensize');
fighand = figure;
clf
set(gcf, 'ToolBar', 'none', 'MenuBar', 'none');
set(fighand, 'CloseRequestFcn', @f_main_close);
set(fighand, 'Name', ['Signal Viewer: ' num2str(size(Signal,2)) ' channels, ' num2str(size(Signal,1)/Fs) ' seconds record duration, ' num2str(Fs) ' Hz sample rate']);
hold on
set(gcf,'Position',screensize);
YLim = [-chansep*Nsch-0.5*chansep, -chansep+0.5*chansep];
XLim = [0 10];
set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim, 'XLim', XLim, 'Position', AxesPosition, 'FontWeight', 'bold', 'FontName', AxesFontName, 'FontSize', AxesFontSize); %#ok<*NBRAK>

t1 = 1;
t2 = 2;
Time = (0:Ntp-1)/Fs;
Time_min = 0;
Time_max = (Ntp-1)/Fs;

for ch = Nsch:-1:1
    plothand(ch) = plot(Time(t1:t2), Signal(t1:t2,selchan(ch)) - nanmean(Signal(t1:t2,selchan(ch))) - chansep*ch);
    set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:));
    set(plothand(ch), 'ButtonDownFcn', @f_plothand_buttondown);
    setappdata(plothand(ch), 'chanind', selchan(ch));
    setappdata(plothand(ch), 'channame', ChanNames{selchan(ch)});
end
if ~isempty(plothand)
    selected_plothand = plothand(1);
end

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
    EventTimeStamps = sortrows(EventTimeStamps);
    EventTimes = cell2mat(EventTimeStamps(:,1));
    YPos = sort(mean(YLim)+diff(YLim)/20*([-8:2:8]), 'descend');
    NYPos = length(YPos);
    for i = size(EventTimeStamps,1):-1:1
        eventplothand(i) = plot( EventTimeStamps{i,1}*[1 1], [-10000000*(Nsch+1), 10000000], '-', 'Color', EventKolor ); %#ok<NASGU>
        eventtexthand(i) = text( EventTimeStamps{i,1}, YPos(mod(i-1,NYPos)+1), EventTimeStamps{i,2}, 'FontSize', EventFontSize);
    end
    set(eventtexthand(i), 'Visible', 'off');
end


HighPassFilter.state = 0;
HighPassFilter.cutoff = 1;
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
ZscoreFilter.on_chansep = chansep;

h_hugetext = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.35 0.5 0.25 0.08], 'String', 'BIG TEXT', 'FontSize', 48, 'Visible', 'off', 'BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1]);
h_bigtext = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.25 0.5 0.45 0.05], 'String', 'BIG TEXT', 'FontSize', 28, 'Visible', 'off');

h_xoomtitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.96 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Horizontal Zoom'); %#ok<NASGU>
h_xzoomout = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.94 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-');
h_xzoomlevel = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.94 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '10 s');
h_xzoomin = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.94 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+');

h_yoomtitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.91 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Vertical Zoom'); %#ok<NASGU>
h_yzoomout = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.89 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-');
h_yzoomlevel = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.89 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '64 ch');
h_yzoomin = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.89 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+');

h_sensititle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.86 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Sensitivity'); %#ok<NASGU>
h_sepup = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.84 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-');
h_sensitivity = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.84 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '100 uV');
h_sepdown = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.84 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+');

h_pantitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.81 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Pan'); %#ok<NASGU>
h_panleft = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.755 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '<');
h_panright = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.755 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '>');
h_panup = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.94 0.78 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '^');
h_pandown = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.94 0.755 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'v');

h_hold_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.72 0.045 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Pause Plotting');
h_hold_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.965 0.72 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');

h_psd_plot = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.70 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Show PSD');
h_autofit = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.70 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Auto Fit');


h_passfilttitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.68 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Filters'); %#ok<NASGU>
h_notch_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.66 0.017 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Notch');
h_notch_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.937 0.66 0.017 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_notch_order = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.954 0.66 0.010 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(NotchFilter.order));
h_notch_orderunit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.964 0.66 0.010 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ord'); %#ok<NASGU>
h_notch_qfactor = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.975 0.66 0.010 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(NotchFilter.qfactor));
h_notch_qfactorunit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.985 0.66 0.005 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Q'); %#ok<NASGU>


h_hpf_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.64 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'HP');
h_hpf_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.64 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_hpf_cutoff = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.64 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(HighPassFilter.cutoff));
h_hpf_unit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.64 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz'); %#ok<NASGU>

h_lpf_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.62 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'LP');
h_lpf_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.62 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_lpf_cutoff = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.62 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(LowPassFilter.cutoff));
h_lpf_unit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.62 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz'); %#ok<NASGU>

h_evf_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.60 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ENV');
h_evf_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.60 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_evf_cutoff = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.60 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(EnvelopeFilter.cutoff));
h_evf_unit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.60 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz'); %#ok<NASGU>

h_zscore_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.58 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Z-Score');
h_zscore_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.95 0.58 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');

h_fastdraw_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.56 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Fast Plot');
h_fastdraw_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.95 0.56 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ON');

h_chansel_title = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.525 0.030 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Plotted Channels'); %#ok<NASGU>
h_chansel_list = uicontrol(gcf, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.92 0.12 0.040, 0.400]);
h_chansel_confirm = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.10 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Plot');
h_chansel_reset = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.94 0.10 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'R');


h_icasel_title = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.96 0.525 0.020 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ICA Comps'); %#ok<NASGU>
h_icasel_list = uicontrol(gcf, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.96 0.12 0.030, 0.400]);
h_icasel_confirm = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.10 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Start');
h_icasel_reset = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.10 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'R');
h_icasel_view_sources = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'S');
h_icasel_view_mixmat = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A');
h_icasel_view_sepmat = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W');

h_axesfont_inc = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.05 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A+');
h_axesfont_dec = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.93 0.05 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A-');
h_windowhsize_inc = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.95 0.05 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W+');
h_windowhsize_dec = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.05 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W-');

h_xspan_text = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.025 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 't range [0,12345]');

set(h_icasel_reset, 'Enable', 'off');
set(h_icasel_view_sources, 'Enable', 'off');
set(h_icasel_view_mixmat, 'Enable', 'off');
set(h_icasel_view_sepmat, 'Enable', 'off');
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
set(h_chansel_confirm, 'Callback', @f_chansel_confirm);
set(h_chansel_reset, 'Callback', @f_chansel_reset);
set(h_psd_plot, 'Callback', @f_psd_plot);
set(h_autofit, 'Callback', @f_autofit);
set(h_icasel_confirm, 'Callback', @f_icasel_confirm);
set(h_icasel_reset, 'Callback', @f_icasel_reset);
set(h_icasel_view_sources, 'Callback', @f_icasel_view_sources);
set(h_icasel_view_mixmat, 'Callback', @f_icasel_view_mixmat);
set(h_icasel_view_sepmat, 'Callback', @f_icasel_view_sepmat);
set(h_axesfont_inc, 'Callback', @f_axesfont_inc);
set(h_axesfont_dec, 'Callback', @f_axesfont_dec);
set(h_windowhsize_inc, 'Callback', @f_windowhsize_inc);
set(h_windowhsize_dec, 'Callback', @f_windowhsize_dec);
set(gcf, 'KeyPressFcn', @f_fig_keypress);
set(gcf,'Position',screensize);


notch_update();
filter_update();
%set(h_bigtext, 'Visible', 'off', 'String', '');

set(h_xspan_text, 'String', ['t range [' num2str(round(min(Time))) ', ' num2str(round(max(Time))) '] s']);

autofit();


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
                    set(gca, 'XLim', XLim);
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
                        fprintf('a');
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
                
    end


    function f_hold_switch(hObject, eventdata) %#ok<*INUSD>
        if PlotHold
            set(h_hold_state, 'String', 'OFF');
            PlotHold = 0;
            set(h_hugetext, 'Visible', 'off');
            redraw();
        else
            set(h_hold_state, 'String', 'ON');
            set(h_hugetext, 'String', 'Plotting Paused', 'Visible', 'on');
            PlotHold = 1;
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
            set(h_hpf_state, 'String', 'ON');
        else
            set(h_hpf_state, 'String', 'OFF');
        end
        
        filter_update();
        enable_filter_switches();
    end

    function f_lpf_switch(hObject, eventdata)
        disable_filter_switches();
        if ~LowPassFilter.state
            LowPassFilter.state = 1;
        else
            LowPassFilter.state = 0;
        end
        
        if LowPassFilter.state
            set(h_lpf_state, 'String', 'ON');
        else
            set(h_lpf_state, 'String', 'OFF');
        end
        
        filter_update();
        enable_filter_switches();
    end

    function f_evf_switch(hObject, eventdata)
        disable_filter_switches();
        if ~EnvelopeFilter.state
            EnvelopeFilter.state = 1;
        else
            EnvelopeFilter.state = 0;
        end
        
        if EnvelopeFilter.state
            set(h_evf_state, 'String', 'ON');
        else
            set(h_evf_state, 'String', 'OFF');
        end
        
        filter_update();
        enable_filter_switches();
    end

    function f_hpf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if HighPassFilter.state
            filter_update();
        end
        enable_filter_switches();
    end

    function f_lpf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if LowPassFilter.state
            filter_update();
        end
        enable_filter_switches();
    end

    function f_evf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if EnvelopeFilter.state
            filter_update();
        end
        enable_filter_switches();
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
        set(h_notch_state, 'String', 'Wait'); drawnow;
        if NotchFilter.state
            NotchFilter.state = 0;
        else
            NotchFilter.state = 1;
        end
        notch_update();
        filter_update();
        enable_filter_switches();
    end

    function f_notch_order(hObject, eventdata)
        disable_filter_switches();
        ord = str2double(get(h_notch_order, 'String'));
        if isfinite(ord) && mod(ord,2) == 0 && ord >= 2 && ord <= 10
            set(h_notch_state, 'String', 'Wait'); drawnow;
            NotchFilter.order = ord;
            notch_update();
            filter_update();
        else
            set(h_notch_order, 'String', num2str(NotchFilter.order));
        end
        enable_filter_switches();
    end

    function f_notch_qfactor(hObject, eventdata)
        disable_filter_switches();
        q = str2double(get(h_notch_qfactor, 'String'));
        if isfinite(q) && q > 0
            set(h_notch_state, 'String', 'Wait'); drawnow;
            NotchFilter.qfactor = q;
            notch_update();
            filter_update();
        else
            set(h_notch_qfactor, 'String', num2str(NotchFilter.qfactor));
        end
        enable_filter_switches();
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
            set(h_notch_state, 'String', 'ON');
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', '');
        else
            Signal_postnotch = Signal_postica;
            set(h_notch_state, 'String', 'OFF');
        end
    end


    function f_zscore_switch(hObject, eventdata)
        set(h_zscore_switch, 'Enable', 'off');
        set(h_zscore_state, 'String', 'Wait'); drawnow;
        if ZscoreFilter.state
            ZscoreFilter.on_chansep = chansep;
            ZscoreFilter.state = 0;
            set(h_zscore_state, 'String', 'OFF');
            refit(ZscoreFilter.off_chansep);
        else
            ZscoreFilter.off_chansep = chansep;
            ZscoreFilter.state = 1;
            set(h_zscore_state, 'String', 'ON');
            refit(ZscoreFilter.on_chansep);
        end
        redraw();
        set(h_zscore_switch, 'Enable', 'on');
    end


    function f_fastdraw_switch(hObject, eventdata)
        set(h_fastdraw_switch, 'Enable', 'off');
        set(h_fastdraw_state, 'String', 'Wait'); drawnow;
        if ScreenLimitedDownsampling
            ScreenLimitedDownsampling = 0;
            set(h_fastdraw_state, 'String', 'OFF');
        else
            ScreenLimitedDownsampling = 1;
            set(h_fastdraw_state, 'String', 'ON');
        end
        redraw();
        set(h_fastdraw_switch, 'Enable', 'on');
    end


    function filter_update()
        hpf = str2double(get(h_hpf_cutoff, 'String'));
        lpf = str2double(get(h_lpf_cutoff, 'String'));
        evf = str2double(get(h_evf_cutoff, 'String'));
        
        if hpf <= 0.01
            HighPassFilter.state = 0;
            set(h_hpf_state, 'String', 'OFF');
        end
        if lpf >= Fs/2
            LowPassFilter.state = 0;
            set(h_lpf_state, 'String', 'OFF');
        end
        if hpf >= lpf
            HighPassFilter.state = 0;
            set(h_hpf_state, 'String', 'OFF');
            LowPassFilter.state = 0;
            set(h_lpf_state, 'String', 'OFF');
        end
        if evf >= Fs/2
            EnvelopeFilter.state = 0;
            set(h_evf_state, 'String', 'OFF');
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
            set([h_hpf_state h_lpf_state h_evf_state], 'String', 'Wait');
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
                            set(h_hpf_state, 'String', 'OFF');
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
                            set(h_lpf_state, 'String', 'OFF');
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
                            set(h_evf_state, 'String', 'OFF');
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
                set(h_hpf_state, 'String', 'ON');
            else
                set(h_hpf_state, 'String', 'OFF');
            end
            if LowPassFilter.state
                set(h_lpf_state, 'String', 'ON');
            else
                set(h_lpf_state, 'String', 'OFF');
            end
            if EnvelopeFilter.state
                set(h_evf_state, 'String', 'ON');
            else
                set(h_evf_state, 'String', 'OFF');
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
        set(gca, 'XLim', XLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
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
        set(gca, 'XLim', XLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
    end

    function f_yzoomin(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YRange = YRange / 2;
        YLim(1) = YLim(2) - YRange;
        set(gca, 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
    end

    function f_yzoomout(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        %YCenter = (YLim(1)+YLim(2))/2;
        YRange = YRange * 2;
        YLim(1) = YLim(2) - YRange;
        set(gca, 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
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
        set(gca, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
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
        set(gca, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
    end

    function f_panup(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YLim = YLim + YRange;
        set(gca, 'YLim', YLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
    end

    function f_pandown(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YLim = YLim - YRange;
        set(gca, 'YLim', YLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
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
        %set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames), 'YLim', YLim);
        YLim(2) = 0.5*chansep - FirstChViewable*chansep;
        YLim(1) = YLim(2) - chansep*Nchviewable - 0.5*chansep;
        set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        %redraw();
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
        %set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames), 'YLim', YLim);
        YLim(2) = 0.5*chansep - FirstChViewable*chansep;
        YLim(1) = YLim(2) - chansep*Nchviewable - 0.5*chansep;
        set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim);
        MovementBusy = 1;
        resnap_zoom();
        MovementBusy = 0;
        %redraw();
    end

    function resnap_pan()
        XRange = XLim(2)-XLim(1);
        YRange = YLim(2)-YLim(1);
        if XLim(1) < Time_min
            XLim(1) = Time_min;
            XLim(2) = XLim(1) + XRange;
        end
        if XLim(2) > Time_max
            XLim(2) = Time_max;
            XLim(1) = XLim(2) - XRange;
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
        set(gca, 'XLim', XLim, 'YLim', YLim);
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
        
        set(gca, 'XLim', XLim, 'YLim', YLim);
        filter_render();
        redraw();
    end

    function redraw()
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
                Signal4 = Signal_postfilter(t1:BLIM:t2,selchan(ch));
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
                    YDATA = zscore(YDATA)*ZscoreFilter.multiplier;
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
                set(gca, 'YColor', BusyYColor);
                drawnow
                tdrawupdate = tic;
            end
            
        end
        for ch = length(plothand):-1:Nsch+1
            set(plothand(ch), 'Visible', 'off');
        end
        
        if PlotHold
            set(gca, 'YColor', InactiveYColor);
        else
            set(gca, 'YColor', DefaultYColor);
        end
        
        [~,in] = min(abs((Time(t2)-Time(t1))./XTickSpacings - 20));
        tm1 = floor(Time(t1)/XTickSpacings(in))*XTickSpacings(in);
        tm2 = ceil(Time(t2)/XTickSpacings(in))*XTickSpacings(in);
        tms = tm1:XTickSpacings(in):tm2;
        tmslabel = cell(1,length(tms));
        for i = 1:length(tms)
            tmslabel{i} = sprintf('%g%s', tms(i)/XTickSpacingsAndUnits{in,2}, XTickSpacingsAndUnits{in,3});
        end
        set(gca, 'XTick', tms, 'XTickLabel', tmslabel);
        
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
        if XRange < 60
            set(h_xzoomlevel, 'String', sprintf('%g s', XRange));
        elseif XRange < 3600
            set(h_xzoomlevel, 'String', sprintf('%g min', XRange/60));
        else
            set(h_xzoomlevel, 'String', sprintf('%g hr', XRange/3600));
        end
        
        YRange = YLim(2)-YLim(1);
        set(h_yzoomlevel, 'String', sprintf('%g ch', YRange/chansep));
        
        tmp_sq = '';
        if EnvelopeFilter.state
            tmp_sq = '�';
        end
        if ZscoreFilter.state
            set(h_sensitivity, 'String', sprintf('%g sdev', chansep/ZscoreFilter.multiplier));
        elseif chansep < 1
            set(h_sensitivity, 'String', sprintf('%g nV%s', chansep*1000, tmp_sq));
        elseif chansep < 1000
            set(h_sensitivity, 'String', sprintf('%g �V%s', chansep, tmp_sq));
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
    end


    function f_chansel_confirm(hObject, eventdata)
        selchan = get(h_chansel_list, 'Value');
        Nsch = length(selchan);
        set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)));
        redraw();
    end

    function f_chansel_reset(hObject, eventdata)
        set(h_chansel_list, 'Value', selchan);
    end

    function f_psd_plot(hObject, eventdata)
        if ~ishandle(viewhand_psd)
            viewhand_psd = figure;
            update_psd();
        else
            figure(viewhand_psd);
        end
    end

    function update_psd()
        if ishandle(viewhand_psd)
            if ~isempty(selected_plothand) && ishandle(selected_plothand)
                channame = getappdata(selected_plothand, 'channame');
                chanind = getappdata(selected_plothand, 'chanind');
                tmp = Signal_postfilter(t1:t2,chanind);
                if size(Signal_psd_source,1) == size(tmp,1) && size(Signal_psd_source,2) == size(tmp,2) && norm(Signal_psd_source - tmp) == 0
                    % These are the same signal.
                    return
                end
                tmp1 = gcf;
                set(0, 'CurrentFigure', viewhand_psd);
                set(gca, 'FontSize', 16);
                Signal_psd_source = tmp;
                [pxx, fxx] = pwelch(Signal_psd_source, [], [], [], Fs);
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
                    set(gca, 'XLim', [fc1 fc2]);
                end
                title(['Welch PSD in ' channame ', ' num2str(Time(t1))  '-' num2str(Time(t2)) ' s']);
                set(viewhand_psd, 'Name', ['PSD in ' channame ' from ' num2str(Time(t1)) ' to ' num2str(Time(t2)) ' s ' tmp_filttext]);
                set(0, 'CurrentFigure', tmp1);
            else
                set(viewhand_psd, 'Name', 'Select a channel first by clicking on its signal.');
            end
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
            [~, ica_A, ica_W] = fastica(Signal.', 'stabilization', 'on', 'maxNumIterations', 200);
            set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
            FilterBusy = 0;
        elseif size(Signal,2) ~= size(ica_W,2) || size(Signal,2) ~= size(ica_A,1)
            % ICA parameters are wrong. Disable it
            set(h_icasel_confirm, 'Enable', 'off', 'String', 'No ICA');
            set(h_icasel_reset, 'Enable', 'off');
            set(h_icasel_view_sources, 'Enable', 'off');
            set(h_icasel_view_mixmat, 'Enable', 'off');
            set(h_icasel_view_sepmat, 'Enable', 'off');
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
            ICA_Initialized = 1;
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
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
            viewhand_ica_sig = signalviewer(ica_sig.', Fs, ic_names);
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
            set(viewhand_ica_A_ax, 'FontSize', 8);
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
            set(viewhand_ica_W_ax, 'FontSize', 8);
            xlabel('Channels');
            ylabel('Sources');
            set(viewhand_ica_W_ax, 'YTick', 1:size(ica_W,1));
            set(viewhand_ica_W_ax, 'XTick', 1:size(ica_W,2));
            set(viewhand_ica_W_ax, 'XTickLabel', ChanNames);
        else
            figure(viewhand_ica_W);
        end
    end


    function f_axesfont_inc(hObject, eventdata)
        AxesFontSize = AxesFontSize + 1;
        set(gca, 'FontSize', AxesFontSize);
    end

    function f_axesfont_dec(hObject, eventdata)
        if AxesFontSize > 1
            AxesFontSize = AxesFontSize - 1;
        end
        set(gca, 'FontSize', AxesFontSize);
    end

    function f_windowhsize_inc(hObject, eventdata)
        %AxesPosition = [0.0500    0.0600    0.86    0.93];
        AxesPosition = AxesPosition + [-0.01 0 0.01 0];
        set(gca, 'Position', AxesPosition);
    end

    function f_windowhsize_dec(hObject, eventdata)
        %AxesPosition = [0.0500    0.0600    0.86    0.93];
        if AxesPosition(3) > 0
            AxesPosition = AxesPosition - [-0.01 0 0.01 0];
        end
        set(gca, 'Position', AxesPosition);
    end

    function f_plothand_buttondown(hObject, eventdata)
        selected_plothand = hObject;
        update_psd();
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
        tmpr = 1+Fs:size(Signal_postfilter,1)-Fs+1;
        if isempty(tmpr)
            tmpr = 1:size(Signal_postfilter,1);
        end
        if ZscoreFilter.state
            tmp = max(max(abs(zscore(Signal_postfilter(tmpr,:))*ZscoreFilter.multiplier)));
        else
            tmp = max(max(abs(Signal_postfilter(tmpr,:))));
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

