% Time can also be the sampling frequency of Signal instead of an array of sample time
% EventTimeStamps should be a Nx2 cell {1.2345, 'Abc'; 4.5, 'Def'} with first column in seconds, and 2nd column a string
% 

function signalviewer(Signal, Time, Chan, EventTimeStamps)
if ~exist('Chan', 'var')
    nchan = size(Signal,2);
    npad = floor(log10(nchan))+1;
    if numel(Time) > 1 && size(Signal,2) == length(Time)
        nchan = size(Signal,1);
    end
    Chan = string_to_cell(num2str(1:nchan,['ch%0' num2str(npad) 'i,']),',');
end
Chan = Chan(:).';
Time = Time(:);
if size(Signal,1) == length(Chan) && size(Signal,2) ~= length(Chan)
    Signal = Signal.';
end
if size(Signal,1) ~= length(Time) && length(Time) == 1
    Fs = Time;
    Time = (0:size(Signal,1)-1).' / Fs;
else
    Fs = round(1/nanmean(diff(Time)));
end

% Remove channels that are completely flat
tmp = nanstd(Signal,[],1);
Signal = Signal(:,tmp>0);
Chan = Chan(tmp>0);

selchan = 1:length(Chan);
selica = [];


[Chan,ix] = sort(Chan);
Signal = Signal(:,ix);
Signal1 = Signal;
Signal2 = Signal;
Signal3 = Signal;
PowerLineFrequency = 60;
MaxFilterOrder = 4;
FilterOrder = MaxFilterOrder;
FilterReflectMult = 0;
FilterChunkSec = 100;
FilterBusy = 0;
ICA_Initialized = 0;
ica_sig = [];
ica_A = [];
icachans = {};

Nch = length(Chan(selchan));
Ntp = length(Time);

PermittedXZoomRanges   = [0.001 0.005 0.01 0.05 0.1 0.5 1 2 5 10 20 30 60 120 300 600 1200 1800 3600:3600:6*3600 8*3600 12*3600 24*3600 7*24*3600];
PermittedXZoomRanges = PermittedXZoomRanges(PermittedXZoomRanges > 16/Fs);
PermittedChanSepRanges = [1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000];
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
    
% Kolor = jet(32)*0.75;
% Kolor = Kolor(reshape(reshape(1:32,[],2).',1,[]),:);

Kolor = [
    0.75 0 0
    0 0.75 0
    0 0 0.75
];
Nkolor = size(Kolor,1);

EventKolor = [0.75 0.75 0];
EventFontSize = 20;

sig_sd = std(Signal,[],1);
sig_sd_max = max(sig_sd);

chansep = 100;

screensize = get(groot,'Screensize');
figure(40301);
clf
hold on
set(gcf,'Position',screensize);
YLim = [-chansep*Nch-0.5*chansep, -chansep+0.5*chansep];
XLim = [0 10];
set(gca, 'YTick', [-chansep*Nch:chansep:-chansep], 'YTickLabel', fliplr(Chan(selchan)), 'YLim', YLim, 'XLim', XLim, 'Position', [0.0500    0.0600    0.86    0.93], 'FontWeight', 'bold');

t1 = find(Time<=0,1);
if isempty(t1)
    % Signal does not start from 0 seconds
    % Pad
    Time = [0; Time];
    Ntp = Ntp + 1;
    Signal = [nan(1,Nch); Signal];
    t1 = 1;
end
t2 = find(Time>=10,1);
if isempty(t2)
    % Signal is less than 10 seconds long
    % Pad
    Time(end+1) = 10;
    Ntp = Ntp + 1;
    Signal(end+1,:) = NaN;
    t2 = find(Time>=10,1);
end

for ch = Nch:-1:1
    plothand(ch) = plot(Time(t1:t2), Signal(t1:t2,selchan(ch)) - nanmean(Signal(t1:t2,selchan(ch))) - chansep*ch);
    set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:));
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
                    EventTimeStamps(k,:) = {tmp.(tmpf{i})(j), tmpf{i}};
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
        eventplothand(i) = plot( EventTimeStamps{i,1}*[1 1], [-10000000*(Nch+1), 10000000], '-', 'Color', EventKolor );
        eventtexthand(i) = text( EventTimeStamps{i,1}, YPos(mod(i-1,NYPos)+1), EventTimeStamps{i,2}, 'FontSize', EventFontSize);
    end
    set(eventtexthand(i), 'Visible', 'off');
end


HighPassFilter.state = 0;
HighPassFilter.cutoff = 1;
LowPassFilter.state = 0;
LowPassFilter.cutoff = 30;
NotchFilter.state = 0;
ZscoreFilter.state = 0;
ZscoreFilter.multiplier = 100;

h_bigtext = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.35 0.5 0.25 0.05], 'String', 'BIG TEXT', 'FontSize', 28, 'Visible', 'off');

h_xoomtitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.96 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Horizontal Zoom');
h_xzoomout = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.94 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-');
h_xzoomlevel = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.94 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '10 s');
h_xzoomin = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.94 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+');

h_yoomtitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.90 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Vertical Zoom');
h_yzoomout = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.88 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-');
h_yzoomlevel = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.88 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '64 ch');
h_yzoomin = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.88 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+');

h_sensititle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.84 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Sensitivity');
h_sepup = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.82 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-');
h_sensitivity = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.82 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '100 uV');
h_sepdown = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.97 0.82 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+');

h_pantitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.78 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Pan');
h_panleft = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.725 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '<');
h_panright = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.725 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '>');
h_panup = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.94 0.75 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '^');
h_pandown = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.94 0.725 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'v');

h_passfilttitle = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.68 0.065 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Filters');
h_hpf_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.66 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'HP');
h_hpf_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.66 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_hpf_cutoff = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.66 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '1');
h_hpf_unit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.66 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz');
h_lpf_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.64 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'LP');
h_lpf_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.935 0.64 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_lpf_cutoff = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.955 0.64 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '30');
h_lpf_unit = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.64 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz');
h_notch_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.62 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Notch');
h_notch_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.95 0.62 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');
h_zscore_switch = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.60 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Z-Score');
h_zscore_state = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.95 0.60 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'OFF');


h_chansel_title = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.92 0.525 0.030 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Plotted Channels');
h_chansel_list = uicontrol(gcf, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.92 0.12 0.040, 0.400]);
h_chansel_confirm = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.10 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Confirm');
h_chansel_reset = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.92 0.08 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Reset');

h_icasel_title = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.96 0.525 0.020 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ICA Comps');
h_icasel_list = uicontrol(gcf, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.96 0.12 0.030, 0.400]);
h_icasel_confirm = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.10 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Run ICA');
h_icasel_reset = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.08 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Reset');

set(h_icasel_reset, 'Enable', 'off');
set(h_chansel_list, 'String', Chan);
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

set(h_hpf_switch, 'Callback', @f_hpf_switch);
set(h_lpf_switch, 'Callback', @f_lpf_switch);
set(h_hpf_cutoff, 'Callback', @f_hpf_cutoff);
set(h_lpf_cutoff, 'Callback', @f_lpf_cutoff);
set(h_notch_switch, 'Callback', @f_notch_switch);
set(h_zscore_switch, 'Callback', @f_zscore_switch);
set(h_chansel_confirm, 'Callback', @f_chansel_confirm);
set(h_chansel_reset, 'Callback', @f_chansel_reset);
set(h_icasel_confirm, 'Callback', @f_icasel_confirm);
set(h_icasel_reset, 'Callback', @f_icasel_reset);
set(gcf, 'KeyPressFcn', @f_fig_keypress);
set(gcf,'Position',screensize);

notch_update();
filter_update();
%set(h_bigtext, 'Visible', 'off', 'String', '');


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
                    end
                elseif Shift
                    f_panleft(hObject, []);
                    f_panleft(hObject, []);
                    f_panleft(hObject, []);
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
                    end
                elseif Shift
                    f_panright(hObject, []);
                    f_panright(hObject, []);
                    f_panright(hObject, []);
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

    function f_hpf_cutoff(hObject, eventdata)
        disable_filter_switches();
        filter_update();
        enable_filter_switches();
    end

    function f_lpf_cutoff(hObject, eventdata)
        disable_filter_switches();
        filter_update();
        enable_filter_switches();
    end

    function disable_filter_switches()
        set([h_hpf_switch, h_lpf_switch, h_hpf_cutoff, h_lpf_cutoff, h_notch_switch], 'Enable', 'off');
        drawnow
    end

    function enable_filter_switches()
        set([h_hpf_switch, h_lpf_switch, h_hpf_cutoff, h_lpf_cutoff, h_notch_switch], 'Enable', 'on');
    end

    function f_notch_switch(hObject, eventdata)
        set(h_notch_switch, 'Enable', 'off');
        set(h_notch_state, 'String', 'Wait'); drawnow;
        if NotchFilter.state
            NotchFilter.state = 0;
        else
            NotchFilter.state = 1;
        end
        notch_update();
        filter_update();
        set(h_notch_switch, 'Enable', 'on');
    end
    
    function notch_update()
        if NotchFilter.state
            set(h_notch_state, 'String', 'ON');
            FilterBusy = 1;
            set(h_bigtext, 'Visible', 'on', 'String', 'Applying notch filter...'); drawnow;
            ETime = (Time(1):1/Fs:Time(end)).';
            ESignal = interp1(Time, Signal1, ETime);
            tmp = freqfilter(ESignal, Fs, [PowerLineFrequency+[-2 2], 2], 'stop', 'butter', 1*Fs);
            Signal2 = interp1(ETime, tmp, Time);
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', '');
            clear ETime ESignal tmp
        else
            set(h_notch_state, 'String', 'OFF');
            Signal2 = Signal1;
        end
    end


    function f_zscore_switch(hObject, eventdata)
        set(h_zscore_switch, 'Enable', 'off');
        set(h_zscore_state, 'String', 'Wait'); drawnow;
        if ZscoreFilter.state
            ZscoreFilter.state = 0;
            set(h_zscore_state, 'String', 'OFF');
        else
            ZscoreFilter.state = 1;
            set(h_zscore_state, 'String', 'ON');
        end
        redraw();
        set(h_zscore_switch, 'Enable', 'on');
    end


    function filter_update()
        hpf = str2double(get(h_hpf_cutoff, 'String'));
        lpf = str2double(get(h_lpf_cutoff, 'String'));
        
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
        
        if HighPassFilter.cutoff ~= hpf || LowPassFilter.cutoff ~= lpf
            HighPassFilter.cutoff = hpf;
            LowPassFilter.cutoff = lpf;
        end
        if HighPassFilter.state || LowPassFilter.state
            SigChunkRendered(:) = 0;
            FilterOrder = MaxFilterOrder;
            filter_render();
            
%             set([h_hpf_state h_lpf_state], 'String', 'Wait'); drawnow
%             ETime = (Time(1):1/Fs:Time(end)).';
%             ESignal = interp1(Time, Signal2, ETime);
%             if HighPassFilter.state
%                 [Signal3a, FilterInfo] = freqfilter(ESignal, Fs, [hpf FilterOrder], 'high', 'butter', FilterReflectMult*FilterOrder/hpf*Fs);
%                 while any(FilterInfo.ButterUnstable)
%                     if FilterOrder <= 1
%                         Signal3a = ESignal;
%                         HighPassFilter.state = 0;
%                         set(h_hpf_state, 'String', 'OFF');
%                         break
%                     end
%                     FilterOrder = ceil(FilterOrder / 2);
%                     [Signal3a, FilterInfo] = freqfilter(ESignal, Fs, [hpf FilterOrder], 'high', 'butter', FilterReflectMult*FilterOrder/hpf*Fs);
%                 end
%             else
%                 Signal3a = ESignal;
%             end
%             
%             if LowPassFilter.state
%                 [Signal3, FilterInfo] = freqfilter(Signal3a, Fs, [lpf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/lpf*Fs);
%                 while any(FilterInfo.ButterUnstable)
%                     if FilterOrder <= 1
%                         Signal3 = Signal3a;
%                         LowPassFilter.state = 0;
%                         set(h_lpf_state, 'String', 'OFF');
%                         break
%                     end
%                     FilterOrder = ceil(FilterOrder / 2);
%                     [Signal3, FilterInfo] = freqfilter(Signal3a, Fs, [lpf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/lpf*Fs);
%                 end
%             else
%                 Signal3 = Signal3a;
%             end
%             Signal3 = interp1(ETime, Signal3, Time);
%             clear ETime ESignal Signal3a
%             
%             if HighPassFilter.state
%                 set(h_hpf_state, 'String', 'ON');
%             else
%                 set(h_hpf_state, 'String', 'OFF');
%             end
%             if LowPassFilter.state
%                 set(h_lpf_state, 'String', 'ON');
%             else
%                 set(h_lpf_state, 'String', 'OFF');
%             end
            
        else
            Signal3 = Signal2;
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
        chunkindexrange = xlim_to_sigchunk();
        if SigChunkRendered(chunkindexrange(1):chunkindexrange(2))
            %redraw();
            return;
        end
            
        if HighPassFilter.state || LowPassFilter.state
            set([h_hpf_state h_lpf_state], 'String', 'Wait');
            redraw();
            
            % Only need to render the chunks and one chunk to the left and
            % to the right
            
            
            e1 = max(Time(1), (chunkindexrange(1)-1-1)*FilterChunkSec);
            e2 = min(Time(end), (chunkindexrange(2)+1)*FilterChunkSec);
            
            ti1 = find(Time>=e1,1);
            ti2 = find(Time<=e2,1,'last');
            
            %fprintf('Rendering time range %g-%g s, Time(%i) to Time(%i)\n', e1, e2, ti1, ti2);
            
            %ETime = (Time(1):1/Fs:Time(end)).';
            ETime = (e1:1/Fs:e2).';
            
            ESignal = interp1(Time, Signal2, ETime);
            
            if HighPassFilter.state
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Applying high-pass filter...'); drawnow;
                [Signal3a, FilterInfo] = freqfilter(ESignal, Fs, [hpf FilterOrder], 'high', 'butter', FilterReflectMult*FilterOrder/hpf*Fs);
                while any(FilterInfo.ButterUnstable)
                    if FilterOrder <= 1
                        Signal3a = ESignal;
                        HighPassFilter.state = 0;
                        set(h_hpf_state, 'String', 'OFF');
                        SigChunkRendered(:) = 0;
                        break
                    end
                    FilterOrder = ceil(FilterOrder / 2);
                    [Signal3a, FilterInfo] = freqfilter(ESignal, Fs, [hpf FilterOrder], 'high', 'butter', FilterReflectMult*FilterOrder/hpf*Fs);
                end
            else
                Signal3a = ESignal;
            end
            
            if LowPassFilter.state
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Applying low-pass filter...'); drawnow;
                [Signal3b, FilterInfo] = freqfilter(Signal3a, Fs, [lpf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/lpf*Fs);
                while any(FilterInfo.ButterUnstable)
                    if FilterOrder <= 1
                        Signal3b = Signal3a;
                        LowPassFilter.state = 0;
                        set(h_lpf_state, 'String', 'OFF');
                        SigChunkRendered(:) = 0;
                        break
                    end
                    FilterOrder = ceil(FilterOrder / 2);
                    [Signal3b, FilterInfo] = freqfilter(Signal3a, Fs, [lpf FilterOrder], 'low', 'butter', FilterReflectMult*FilterOrder/lpf*Fs);
                end
            else
                Signal3b = Signal3a;
            end
            Signal3(ti1:ti2,:) = interp1(ETime, Signal3b, Time(ti1:ti2));
            clear ETime ESignal Signal3a Signal3b
            
            
            
            SigChunkRendered(chunkindexrange(1):chunkindexrange(2)) = 1;
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', '');
            
            
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
        end
        
        redraw();
    end











    function f_xzoomin(hObject, eventdata)
        if FilterBusy
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
        resnap_zoom();
    end

    function f_xzoomout(hObject, eventdata)
        if FilterBusy
            return;
        end
        XRange = XLim(2)-XLim(1);
        [~,in] = min(abs(XRange - PermittedXZoomRanges));
        if in < length(PermittedXZoomRanges)
            XRange = PermittedXZoomRanges(in+1);
        else
            XRange = PermittedXZoomRanges(in);
        end
        
        if XRange > Time(end)
            XRange = Time(end);
        end
        
        XLim(2) = XLim(1) + XRange;
        set(gca, 'XLim', XLim);
        resnap_zoom();
    end

    function f_yzoomin(hObject, eventdata)
        if FilterBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YRange = YRange / 2;
        YLim(1) = YLim(2) - YRange;
        set(gca, 'YLim', YLim);
        resnap_zoom();
    end

    function f_yzoomout(hObject, eventdata)
        if FilterBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YCenter = (YLim(1)+YLim(2))/2;
        YRange = YRange * 2;
        YLim(1) = YLim(2) - YRange;
        set(gca, 'YLim', YLim);
        resnap_zoom();
    end

    function f_panleft(hObject, eventdata)
        if FilterBusy
            return;
        end
        XRange = XLim(2)-XLim(1);
        XLim = XLim - XRange;
        set(gca, 'XLim', XLim);
        resnap_pan();
    end

    function f_panright(hObject, eventdata)
        if FilterBusy
            return;
        end
        XRange = XLim(2)-XLim(1);
        XLim = XLim + XRange;
        set(gca, 'XLim', XLim);
        resnap_pan();
    end

    function f_panup(hObject, eventdata)
        if FilterBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YLim = YLim + YRange;
        set(gca, 'YLim', YLim);
        resnap_pan();
    end

    function f_pandown(hObject, eventdata)
        if FilterBusy
            return;
        end
        YRange = YLim(2)-YLim(1);
        YLim = YLim - YRange;
        set(gca, 'YLim', YLim);
        resnap_pan();
    end

    function f_sepup(hObject, eventdata)
        if FilterBusy
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
        %YLim = [-chansep*Nch-0.5*chansep, -chansep+0.5*chansep];
        %set(gca, 'YTick', [-chansep*Nch:chansep:-chansep], 'YTickLabel', fliplr(Chan), 'YLim', YLim);
        YLim(2) = 0.5*chansep - FirstChViewable*chansep;
        YLim(1) = YLim(2) - chansep*Nchviewable - 0.5*chansep;
        set(gca, 'YTick', [-chansep*Nch:chansep:-chansep], 'YTickLabel', fliplr(Chan(selchan)), 'YLim', YLim);
        resnap_zoom();
        %redraw();
    end

    function f_sepdown(hObject, eventdata)
        if FilterBusy
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
        %YLim = [-chansep*Nch-0.5*chansep, -chansep+0.5*chansep];
        %set(gca, 'YTick', [-chansep*Nch:chansep:-chansep], 'YTickLabel', fliplr(Chan), 'YLim', YLim);
        YLim(2) = 0.5*chansep - FirstChViewable*chansep;
        YLim(1) = YLim(2) - chansep*Nchviewable - 0.5*chansep;
        set(gca, 'YTick', [-chansep*Nch:chansep:-chansep], 'YTickLabel', fliplr(Chan(selchan)), 'YLim', YLim);
        resnap_zoom();
        %redraw();
    end


    function resnap_pan()
        XRange = XLim(2)-XLim(1);
        YRange = YLim(2)-YLim(1);
        if XLim(1) < Time(1)
            XLim(1) = Time(1);
            XLim(2) = XLim(1) + XRange;
        end
        if XLim(2) > Time(end)
            XLim(2) = Time(end);
            XLim(1) = XLim(2) - XRange;
        end
        if YLim(1) < -chansep*Nch-chansep/2
            YLim(1) = -chansep*Nch-chansep/2;
            YLim(2) = YLim(1) + YRange;
        end
        if YLim(2) > -chansep/2
            YLim(2) = -chansep/2;
            YLim(1) = YLim(2) - YRange;
        end
        
        
        if XRange >= 1
            XLim(1) = floor(XLim(1));
        elseif XRange >= 0.1
            XLim(1) = floor(XLim(1)*10)/10;
        elseif XRange >= 0.01
            XLim(1) = floor(XLim(1)*100)/100;
        elseif XRange >= 0.001
            XLim(1) = floor(XLim(1)*1000)/1000;
        elseif XRange >= 0.0001
            XLim(1) = floor(XLim(1)*10000)/10000;
        elseif XRange >= 0.00001
            XLim(1) = floor(XLim(1)*100000)/100000;
        elseif XRange >= 0.000001
            XLim(1) = floor(XLim(1)*1000000)/1000000;
        end
        XLim(2) = XLim(1) + XRange;
        set(gca, 'XLim', XLim, 'YLim', YLim);
        filter_render();
        redraw();
    end

    function resnap_zoom()
        XRange = XLim(2)-XLim(1);
        YRange = YLim(2)-YLim(1);
        
        if YRange > chansep/2 - (-chansep*Nch-chansep/2)
            YRange = chansep/2 - (-chansep*Nch-chansep/2);
            YLim(1) = YLim(2) - YRange;
        end
        if YRange < chansep
            YRange = chansep;
            YLim(1) = YLim(2) - YRange;
        end
        
        YRange = ceil(YRange/chansep)*chansep;
        YLim(1) = YLim(2) - YRange;

        if XRange >= 1
            XLim(1) = floor(XLim(1));
        elseif XRange >= 0.1
            XLim(1) = floor(XLim(1)*10)/10;
        elseif XRange >= 0.01
            XLim(1) = floor(XLim(1)*100)/100;
        elseif XRange >= 0.001
            XLim(1) = floor(XLim(1)*1000)/1000;
        elseif XRange >= 0.0001
            XLim(1) = floor(XLim(1)*10000)/10000;
        elseif XRange >= 0.00001
            XLim(1) = floor(XLim(1)*100000)/100000;
        elseif XRange >= 0.000001
            XLim(1) = floor(XLim(1)*1000000)/1000000;
        end
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
        Nch = length(selchan);
        
        for ch = Nch:-1:1
            if -chansep/2-chansep*ch < YLim(1) || chansep/2-chansep*ch > YLim(2)
                % out of plotting range
                %fprintf('%i is out of range\n', ch);
                set(plothand(ch), 'Visible', 'off');
                continue;
            end
            YDATA = Signal3(t1:t2,selchan(ch));
            if ZscoreFilter.state
                YDATA = zscore(YDATA)*ZscoreFilter.multiplier;
            else
                YDATA = YDATA - nanmedian(YDATA);
            end
            set(plothand(ch), 'XData', Time(t1:t2), 'YData', YDATA - chansep*ch, 'Visible', 'on');
            set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:));
        end
        for ch = length(plothand):-1:Nch+1
            set(plothand(ch), 'Visible', 'off');
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
        
        if ZscoreFilter.state
            set(h_sensitivity, 'String', sprintf('%g sdev', chansep/ZscoreFilter.multiplier));
        elseif chansep < 1000
            set(h_sensitivity, 'String', sprintf('%g µV', chansep));
        elseif chansep < 1e6
            set(h_sensitivity, 'String', sprintf('%g mV', chansep/1000));
        else
            set(h_sensitivity, 'String', sprintf('%g V', chansep/1000000));
        end
        
        if EventEnable
            l = XLim(1) <= EventTimes & EventTimes <= XLim(2);
            set(eventtexthand(l), 'Visible', 'on');
            set(eventtexthand(~l), 'Visible', 'off');
        end
    end


    function f_chansel_confirm(hObject, eventdata)
        selchan = get(h_chansel_list, 'Value');
        Nch = length(selchan);
        set(gca, 'YTick', [-chansep*Nch:chansep:-chansep], 'YTickLabel', fliplr(Chan(selchan)));
        redraw();
    end

    function f_chansel_reset(hObject, eventdata)
        set(h_chansel_list, 'Value', selchan);
    end

    function f_icasel_confirm(hObject, eventdata)
        if FilterBusy
            return;
        end
        if ~ICA_Initialized
            FilterBusy = 1;
            set(h_icasel_confirm, 'Enable', 'off', 'String', 'Wait');
            set(h_bigtext, 'Visible', 'on', 'String', 'Calculating ICA...'); drawnow;
            [ica_sig, ica_A, ~] = fastica(Signal.', 'stabilization', 'on', 'maxNumIterations', 200);
            mica = size(ica_A,2);
            selica = mica;
            icachans = string_to_cell(num2str(1:mica, 'IC%i '), ' ');
            set(h_icasel_list, 'String', icachans);
            set(h_icasel_list, 'Value', 1:mica);
            set(h_icasel_confirm, 'String', 'Confirm', 'Enable', 'on');
            set(h_icasel_reset, 'Enable', 'on');
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
                Signal1 = Signal;
                FilterBusy = 0;
                set(h_bigtext, 'Visible', 'off', 'String', '');
            else
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Mixing new ICs...'); drawnow;
                tmp = ica_A;
                tmp(:,setdiff(1:size(ica_A,2),selica)) = 0;
                Signal1 = (tmp*ica_sig).';
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

end



function c = string_to_cell(s,d)
% Convert a delimited string to a cell array
% E.g., input is    "blah 1" "blah 2", delimiter is ",
%           output:    {'blah 1', 'blah 2'}

% Copyright 2001-2004 The MathWorks, Inc.
% $Revision: 1.1.8.1 $ $Date: 2004/07/21 06:23:56 $

c = {};
while containsValidString(s),
    [s1 s] = strtok(s, d);
    if containsValidString(s1)
        c = {c{:} s1};
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

