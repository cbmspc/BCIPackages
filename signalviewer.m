% This program expects Signal to be in microvolts and SampleRate to be in
% samples per second
%
% Optional input arguments:
% First, create a struct called opts. Every field is optional.
%
%   opts.EventTimeStamps can accept many formats.
%    PRIMARY SYNTAX: This should be a Nx2 cell {1.2345, 'Abc'; [4.5 5.0], 'Def'} 
%     with first column in seconds (either a scalar for start time or a 
%      1x2 numeric array for start and end times, can be a mix of both), 
%      and 2nd column a string
%     EXAMPLE: opts.EventTimeStamps = {25, 'Experiment start'; [45 65], 'Idle Epoch'; 103.2, 'Stim 1'; [105 115], 'Stim 1 Response'};
%     Hint: If the signal is stitched from multiple epochs, EventTimeStamps
%           must have the same length as the number of epochs, as only the
%           label portion of EventTimeStamps is used (the event times are
%           overwritten by the epoch start times)
%
%    ALTERNATE FORMAT SYNTAX (as examples below):
%       opts.EventTimeStamps = [30 40; 45 45; 50 60; 60 70];
%           Program automatically labels them Event 1, Event 2, ... 
%            Any event with 1st column equal to 2nd column (zero duration) is drawn as a vertical line.
%            Any event with 2nd column greater than 1st column is drawn as a patch.
%            Any event with 1st column greater than 2nd column will have its 2nd column set to the 1st column.
%       opts.EventTimeStamps = [10 20 0; 20 30 1; 30 40 0; 40 50 1; 50 60 0];
%           Program automatically labels them using the 3rd column (in this
%            example, as Event 0, Event 1, Event 0, Event 1, Event 0).
%           Hint: If you only want to label durationless events, you still need three columns. 
%                 Just duplicate the 1st column to the 2nd column, and use the 3rd column as labels.
%
%   opts.EventColors = RGB triplets with the same number of rows as opts.EventTimeStamps
%     Instead of using automatic coloring for event lines/patches and labels, you can specify 
%       the color of line (or patch) and text label for each event. Note that the colors will 
%       still be whitened according to the transparency settings.
%
%   opts.EventColorRules = Regular Expression rules (Nx2 cell array: first column is the rule, second column is the RGB triplet) to automatically color based on the event names
%     Example: opts.EventColorRules = {'^Idle' [1.0 0.5 0.0]; '^Move' [0 0 1]}
%
%   opts.EventLineRules = Regular Expression rules (Nx3 cell array: first column is the rule, second column is the LineWidth (a number), third column is the LineStyle (a character)) to automatically set based on the event names
%     Example: opts.EventLineRules = {'^Idle' 3 ':'; '^Move' 2 '--'}
%
%   opts.EventFontSize = Event font size in points.
%     Example: opts.EventFontSize = 16
%
%   opts.ica_W = ICA separating matrix with orientation (Nsource x Nchan), i.e. ica_W * Signal.' = Source.'
%   opts.ica_A = ICA mixing matrix with orientation (Nchan x Nsource), i.e. ica_A * Source.' = Signal.'
%    Example: [ica_sig, ica_A, ica_W] = fastica(Signal.', ...
%                   'stabilization', fastica_stabilization, ...
%                   'maxNumIterations', fastica_maxNumIterations, ...
%                   'approach', fastica_approach, ...
%                   'g', fastica_g, ...
%                   'interactivePCA', fastica_interactivePCA); 
%              These optional paramters can be set in opts.
%    If ica matrices are not specified, FastICA is used
%
%   opts.ltmat = One or more linear transformation matrix with orientation
%   (NchanOUT x NchanIN), such that transformedSignal.' = ltmat * originalSignal.'
%     Use the 3rd dimension to specify more than one transformation matrix:
%       Example: ltmat(:,:,1) = randn(15,15); ltmat(:,:,2) = randn(15,15)*10;
%   opts.ltmat_do_not_detrend = 0 or 1
%     By default, signal is detrended first to remove the mean (DC shift)
%     for LTMAT. Set this to 1 to disable automatic detrending.
%     
%     For example, if signal has 15 channels, ltmat transforms the input
%     one at a time: out(15x1) = ltmat(15x15) * in(15x1)
%   If this was ICA, ltmat = ica_A * ica_W, because 
%     1. Source.' = ica_W * Signal.'
%     2. Signal.' = ica_A * Source.'
%     3. transformedSignal.' = ica_A * ica_W * originalSignal.'
%     4. Hence, ltmat replaces ica_A * ica_W
%    Because the ltmat option replaces ICA, specifying a valid ltmat is
%    treated as if ica_W and ica_A are both specified in the options, where
%    ica_A = ltmat and ica_W = identity matrix.
%
%   opts.FooterMessage = 'footer messaage'
%
%   opts.pwelch_nfft = 4096
%   NFFT for pwelch
%
%   opts.bridgenans_method = 'linear'
%   Interpolation method (see documentation for interp1) used to fill in
%   the NaN data points
%
%   opts.stitch_mult = 2
%   Stitching algorithm does a high pass on each epoch before joining them
%
%   opts.nanaround_stitch_samples = 4
%   This many samples at and before/after each epoch's edge are NaN'd.
%   This changes the signal.
%
%   opts.blankaround_stitch_samples = 9
%   This many samples to the left and to the right of the signal are NaN'd to
%   hide the stitching artifact for display only; does not affect filtering. 
%
%   opts.nofiltchannames = a cell array of channel names
%   If specified, will not allow CAR, notch, bandpass, or envelope on these
%   channel names (case sensitive channel name matching)
%
%   opts.selectchannames = a cell array of channel names
%   If specified, select the channel names to plot after viewer is opened
%   By default, the viewer opens with all channels plotted.
%   
%   opts.car = 1
%   If specified, turns on Common Average Reference after the viewer is
%   opened. If you want to specify which channels to CAR, select the
%   channels first using opts.selectchannames
%
%   opts.notch = 60
%   If specified, turns on the notch filter at the specified frequency and
%   the 2nd and 3rd harmonics. Order and Q factor can be adjusted after the
%   viewer is loaded. Frequency can be changed in the viewer by clicking
%   the Notch button.
%
%   opts.notch = [60 120 180]
%   If a list is specified, turns on the notch filter at these frequencies.
%   If a frequency is above Nyquist, it is automatically aliased.
%
%   opts.bandpass = [1 35]
%   If specified, turns on Butterworth band-pass filter at the specified
%   corner frequencies. If one of the corner frequencies is not valid, it
%   degrades into either a high-pass filter, a low-pass filter, or the
%   filter is turned off.
%
%   opts.envelope = 2
%   If specified, turns on envelope filter at the specified frequency. This
%   is implemented by squaring and then low-pass filtering.
%
%   opts.zscore = 1
%   If specified, turns on zscore display after the viewer is opened
%
%   opts.sensitivity = 100
%   If specified, sets vertical sensitivity (in microvolts) after viewer is opened
%   By default, the viewer autoscales after applying filters. This can be disabled by setting opts.sensitivity = -1
%
%   opts.xlim = [0 100]
%   If specified, sets the horizontal time range (in seconds) after viewer is opened
%   By default, the viewer opens with the first 10 seconds displayed.
%   If [-inf inf] is specified, sets the range to the entire signal.
%
%   opts.psd_ylim_auto = 1
%   If specified, PSD's ylim will automatically expand and shrink when
%   navigating across different channels. Otherwise, the default behavior
%   is to only expand the vertical (spectral power dB) and never shrink
%   until the "PSD" button is clicked again
%
%   opts.use_jet_colors = 1
%   If specified, the jet color map is used for line colors instead of the "lines" map.
%
%   opts.fastdraw = 1
%   If specified, turns on "Fast Draw" that subsamples the signal before
%   plotting to speed up. Some details are lost in Fast Draw.
%
%   opts.saveas = full path to a file ending in .png
%   If specified, at the end of opening the signal viewer, automatically
%   saves the current view into the png file. The file must not already exist.
%
%   opts.screensize = [1 1 2560 1440]
%   If specified, override the default screensize detection
%   get(groot,'ScreenSize'). This is useful if you want to screenshots to
%   have the same size between different computers
%
%   Signal Hash = The hash (currently using the MD5 hashing algorithm) of
%   the input signal after correcting for orientation, stitching, etc. 
%   Note that SampleRate, ChanNames, Events, and other options are not
%   parts of the hash. Two signals with the same hash are extremely likely
%   to be identical. Signal Hash is displayed on the title bar.
%


function argout1 = signalviewer(Signal, SampleRate, ChanNames, opts, argin5, argin6, argin7)
%t_program_start = tic;
try
    mfnp = mfilename('fullpath');
    if isempty(regexp(mfnp, '\\.m$', 'match'))
        mfnp = [mfnp '.m'];
    end
    tmp = dir(mfnp);
    lmoddate = datestr(tmp.datenum, 'YYmmdd.HHMM'); %#ok<*DATST>
    clear tmp
catch
    lmoddate = '';
end

export_workspace_dest = 'base';

if ~iscell(Signal) && isnumeric(Signal) && size(Signal,3) > 1
    % Assume ch x time x trial
    tmp = cell(1,size(Signal,3));
    for i = 1:size(Signal,3) %#ok<*FXUP>
        tmp{i} = Signal(:,:,i).';
    end
    Signal = tmp;
end

pwelch_nfft = [];

panfrac_noncursor_mode = 1.00; % Default is panning left/right by 100% of the screen. This can be changed on the GUI
% default to pan left/right this many screens

panfrac_cursor_mode = 0.80;
% default to pan left/right this many screens when Cursor mode is ON

hpanamount_mode = 0; % 0 = fraction of screen, 1 = seconds

bridgenans_method = 'pchip';
% See methods in interp1 documentation

stitch_mult = 2;
% Stitching algorithm does a high pass on each epoch before joining them
% To turn the stitching algorithm off, pass in opts.stitch_mult = 0

nanaround_stitch_samples = ceil(0.0125 * SampleRate);
% This many samples to the left and to the right of the signal are NaN'd,
% including the edge itself (meaning if you set this to 0, nothing is
% NaN'd), and the resulting signal is bridged_nan'd again using the method
% specified above; this changes the signal. 

blankaround_stitch_samples = ceil(0.02 * SampleRate);
% This many samples to the left and to the right of the signal are NaN'd to
% hide the stitching artifact for display only; does not affect filtering. 
% Example: blankaround_stitch_samples = 15
% To turn the blanking off, pass in opts.blankaround_stitch_samples = 0

set_chansep_to = 0;
% After the viewer is opened, scale to this vertical sensitivity

set_xlim_to = [0 10];
% After the viewer is opened, set to this xlim

set_selectchannames_to = {};

fastica_stabilization = 'on';
fastica_maxNumIterations = 200;
fastica_approach = 'defl';
fastica_g = 'pow3';
fastica_interactivePCA = 'off';
disable_ica = false;

ZoomedInMarker = 'x';
ZoomedOutMarker = 'none';
MarkerZoomThreshold = 100;

EventTextMinDistApart = 0.00; %Po240625 was 0.05, but text labels are opaque now so it doesn't matter anymore

DecimalSecondsResolution = ceil(log10(SampleRate));

VerticalOverlapAllow = 100;
% 1 = Each channel stays in its own lane. 
% 2 = Each channel can go 50% above or below.
% 100 = Pretty much uncapped

% Default channel names to be excluded from filtering
nofiltchannames = {'TRIGGER', 'TRIG', 'D1', 'D2', 'DIGITAL', 'DIAG'};

isfirsttime_commandentry = 1;

infolabels_deviations = [];
infolabels_text_highdeviations = 'Channel signal deviations (sorted high to low, double click to copy): ';

lookradius_fraction = 0.015;
psd_ylim_auto = 0;
psd_xscale = 'linear';
use_jet_colors = 0;
fastdraw = 0;
psd_dblim_allchans = [inf -inf];
psd_ylim_restrict = [-260 220];
use_ltmat = 0;
ltmat_do_not_detrend = 0;

Ctrl = 0;
Alt = 0;
Shift = 0;
centered_on_eventnum = 0;

SavedPointsTable = {};
% column 1: time in seconds
% column 2: channel index
% column 3: channel name
% column 4: numerical value (after envelope filter)
% column 5: whether this is min(2) or max(3)





if nargin >= 4 && exist('opts', 'var') && isstruct(opts)

    if isfield(opts, 'fastica_stabilization') && ~isempty(opts.fastica_stabilization) && ischar(opts.fastica_stabilization)
        fastica_stabilization = opts.fastica_stabilization;
    end

    if isfield(opts, 'fastica_maxNumIterations') && ~isempty(opts.fastica_maxNumIterations) && isnumeric(opts.fastica_maxNumIterations) && isfinite(opts.fastica_maxNumIterations) && opts.fastica_maxNumIterations > 0
        fastica_maxNumIterations = opts.fastica_maxNumIterations;
    end

    if isfield(opts, 'fastica_approach') && ~isempty(opts.fastica_approach) && ischar(opts.fastica_approach)
        fastica_approach = opts.fastica_approach;
    end

    if isfield(opts, 'fastica_g') && ~isempty(opts.fastica_g) && ischar(opts.fastica_g)
        fastica_g = opts.fastica_g;
    end

    if isfield(opts, 'fastica_interactivePCA') && ~isempty(opts.fastica_interactivePCA) && ischar(opts.fastica_interactivePCA)
        fastica_interactivePCA = opts.fastica_interactivePCA;
    end
    
    if isfield(opts, 'disable_ica') && isscalar(opts.disable_ica) && opts.disable_ica
        disable_ica = true;
    end

    if isfield(opts, 'SavedPointsTable') && iscell(SavedPointsTable) && size(SavedPointsTable,2) == 5 && size(SavedPointsTable,1) >= 1
        SavedPointsTable = opts.SavedPointsTable;
    end

    if isfield(opts, 'EventTimeStamps')
        EventTimeStamps = opts.EventTimeStamps;
    elseif isfield(opts, 'eventtimestamps')
        EventTimeStamps = opts.eventtimestamps;
    else
        clear EventTimeStamps;
    end
    if isfield(opts, 'ltmat') && size(opts.ltmat,1) == size(opts.ltmat,2) && size(opts.ltmat,1) > 0 && isnumeric(opts.ltmat)
        opts.ica_A = opts.ltmat(:,:,1);
        opts.ica_W = eye(size(opts.ica_A,1));
        ltmat = opts.ltmat;
        use_ltmat = 1;
        ltmat_do_not_detrend = 0;
        if isfield(opts, 'ltmat_do_not_detrend') && isscalar(opts.ltmat_do_not_detrend)
            if opts.ltmat_do_not_detrend
                ltmat_do_not_detrend = 1;
            end
        end
    else
        ltmat = [];
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
    if isfield(opts, 'stitch_mult') && isnumeric(opts.stitch_mult) && isscalar(opts.stitch_mult) && isfinite(opts.stitch_mult) && opts.stitch_mult >= 0
        stitch_mult = opts.stitch_mult;
    end
    if isfield(opts, 'bridgenans_method') && ~isempty(opts.bridgenans_method) && ischar(opts.bridgenans_method)
        bridgenans_method = opts.bridgenans_method;
    end
    if isfield(opts, 'nanaround_stitch_samples') && isnumeric(opts.nanaround_stitch_samples) && isscalar(opts.nanaround_stitch_samples) && isfinite(opts.nanaround_stitch_samples) && opts.nanaround_stitch_samples >= 0
        nanaround_stitch_samples = opts.nanaround_stitch_samples;
    end
    if isfield(opts, 'blankaround_stitch_samples') && isnumeric(opts.blankaround_stitch_samples) && isscalar(opts.blankaround_stitch_samples) && isfinite(opts.blankaround_stitch_samples) && opts.blankaround_stitch_samples >= 0
        blankaround_stitch_samples = opts.blankaround_stitch_samples;
    end
    if isfield(opts, 'sensitivity') && isnumeric(opts.sensitivity) && isscalar(opts.sensitivity) && isfinite(opts.sensitivity) && opts.sensitivity >= 0
        set_chansep_to = opts.sensitivity;
    end
    if isfield(opts, 'xlim')
        if isnumeric(opts.xlim) && numel(opts.xlim) == 2
            if isfinite(opts.xlim(1)) && isfinite(opts.xlim(2)) && opts.xlim(1) < opts.xlim(2)
                set_xlim_to = opts.xlim;
            elseif opts.xlim(1) == -inf && opts.xlim(2) == inf
                set_xlim_to = [-inf inf];
            end
        elseif ischar(opts.xlim) && strcmpi(opts.xlim,'all')
            set_xlim_to = [-inf inf];
        end
    end
    if isfield(opts, 'psd_ylim_auto') && isscalar(opts.psd_ylim_auto)
        if opts.psd_ylim_auto
            psd_ylim_auto = 1;
        end
    end
    if isfield(opts, 'use_jet_colors') && isscalar(opts.use_jet_colors)
        if opts.use_jet_colors
            use_jet_colors = 1;
        end
    end
    if isfield(opts, 'fastdraw') && isscalar(opts.fastdraw)
        if opts.fastdraw
            fastdraw = 1;
        end
    end
    if isfield(opts, 'nofiltchannames') && iscell(opts.nofiltchannames)
        nofiltchannames = opts.nofiltchannames;
    end
    if isfield(opts, 'selectchannames') && iscell(opts.selectchannames)
        set_selectchannames_to = opts.selectchannames;
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


if iscell(Signal) && numel(Signal) ~= length(Signal)
    %2024-10-15 If Signal is a cell, force it to be 1D
    Signal = Signal(:).';
end


% This whole if statement happens if Signal is meant to be stitched, i.e.
% it came in as a cell (each cell element is an epoch), or if it came in as
% a 3D array (3rd dimension being the epochs).
% EventTimeStamps must have the same length as Signal
SignalIsStitched = 0;
if iscell(Signal)
    NumEpochs = length(Signal);
    if exist('EventTimeStamps', 'var') && length(EventTimeStamps) == NumEpochs
        % Signal comprises epochs intended to be stitched together.
        % EventTimeStamps must be an array with the same length as the
        % number of epochs, or else it is disregarded completely. 

        if size(EventTimeStamps,2) == NumEpochs && size(EventTimeStamps,1) == 2
            % User has transposed EventTimeStamps. Fix it.
            EventTimeStamps = EventTimeStamps.';
        end

        if iscell(EventTimeStamps) && numel(EventTimeStamps) == length(EventTimeStamps) && numel(EventTimeStamps) == NumEpochs && NumEpochs > 0 && min(cellfun(@ischar,EventTimeStamps))
            % Special case when the user-supplied EventTimeStamps is a 1D
            % cell array of strings. In this case, it is assumed that the
            % user intends to label each epoch with the corresponding text.
            % Does not matter if the user-supplied EventTimeStamps is
            % horizontal or vertical.
            StitchedEpochNames = EventTimeStamps;
            EventTimeStamps = cell(NumEpochs,2);
            EventTimeStamps(:,2) = StitchedEpochNames;
        elseif iscell(EventTimeStamps) && numel(EventTimeStamps) == length(EventTimeStamps) && numel(EventTimeStamps) == NumEpochs && NumEpochs > 0 && min(cellfun(@isnumeric,EventTimeStamps))
            % Special case when the user-supplied EventTimeStamps is a 1D
            % cell array of numbers. In this case, it is assumed that the
            % user intends to label each epoch with the corresponding
            % numbers. Does not matter if the user-supplied EventTimeStamps
            % is horizontal or vertical. This is similar to the case above,
            % but with the extra step of converting to strings
            StitchedEpochNames = cellfun(@num2str,EventTimeStamps);
            EventTimeStamps = cell(NumEpochs,2);
            EventTimeStamps(:,2) = StitchedEpochNames;
        elseif isnumeric(EventTimeStamps) && numel(EventTimeStamps) == length(EventTimeStamps) && numel(EventTimeStamps) == NumEpochs
            % Special case when the user-supplied EventTimeStamps is a 1D
            % numeric array. In this case, it is assumed that the user
            % intends to label each epoch with the corresponding numbers.
            % Does not matter if the user-supplied EventTimeStamps is
            % horizontal or vertical. This is similar to the two cases
            % above, but with the extra step of converting to strings
            StitchedEpochNames = cellfun(@num2str,num2cell(EventTimeStamps),'UniformOutput',false);
            EventTimeStamps = cell(NumEpochs,2);
            EventTimeStamps(:,2) = StitchedEpochNames;
        else
            % Not a completely valid entry. At least fix the format.
            if iscell(EventTimeStamps)
                % Salvage what we can
                if size(EventTimeStamps,2) > 2
                    EventTimeStamps = EventTimeStamps(:,1:2);
                elseif size(EventTimeStamps,2) == 1
                    EventTimeStamps(:,2) = EventTimeStamps(:,1);
                end
                if size(EventTimeStamps,1) > NumEpochs
                    EventTimeStamps = EventTimeStamps(1:NumEpochs,:);
                elseif size(EventTimeStamps,1) < NumEpochs
                    EventTimeStamps(size(EventTimeStamps,1)+1:NumEpochs,:) = {};
                end
            else
                % Completely rebuild
                warning('EventTimeStamps does not have the same length as the number of epochs and is ignored.')
                EventTimeStamps = cell(NumEpochs,2);
                EventTimeStamps(:,2) = string_to_cell(sprintf('Epoch %i|', 1:NumEpochs), '|');
            end
        end

        for i = 1:NumEpochs
            if ~ischar(EventTimeStamps{i,2})
                if isempty(EventTimeStamps{i,2})
                    EventTimeStamps{i,2} = sprintf('Epoch %i', i);
                elseif isnumeric(EventTimeStamps{i,2})
                    EventTimeStamps{i,2} = num2str(EventTimeStamps{i,2});
                else
                    EventTimeStamps{i,2} = sprintf('Epoch %i', i);
                end
            end
        end
        
    else
        % User supplied EventTimeStamps does not meet specification.
        if exist('EventTimeStamps', 'var') && ~isempty(EventTimeStamps)
            warning('EventTimeStamps does not have the same length as the number of epochs and is ignored.')
        end
        EventTimeStamps = cell(NumEpochs,2);
        EventTimeStamps(:,2) = string_to_cell(sprintf('Epoch %i|', 1:NumEpochs), '|');
    end

    EventTimeStamps_timeportion = nan(NumEpochs,2); %zeros(length(Signal),2);
    tmp = 0;
    for i = 1:NumEpochs
        if size(Signal{i},1) == 0
            % Skip labeling this empty event, but still increment count
            continue
        end
        EventTimeStamps_timeportion(i,1:2) = [tmp, tmp + size(Signal{i},1) / SampleRate];
        tmp = tmp + size(Signal{i},1) / SampleRate;
    end
    
    % 2022-05-24: Keep the time points in each cell before merging
    SignalIsStitched = 1;
    EventTimePoints = EventTimeStamps_timeportion(:,1:2) * SampleRate;
    EventTimePoints(:,1) = EventTimePoints(:,1) + 1;
    
    EventTimeStamps(:,1) = num2cell(EventTimeStamps_timeportion(:,1)); % Only care about the start times
    
    Signal = cat(1,Signal{:});
else
    EventTimePoints = [];
end

if isempty(Signal)
    error('Signal is empty');
end

Signal = double(Signal);

if exist('ChanNames', 'var') && isstruct(ChanNames)
    error('ChanNames must be specified (even if it is empty as {}). Correct syntax: signalviewer(Signal, SampleRate, ChanNames, Opts).');
end

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
    ChanNames = string_to_cell(num2str(1:nchan,['c%0' num2str(npad) 'i,']),',');
end
ChanNames = ChanNames(:).';
if size(Signal,1) == length(ChanNames) && size(Signal,2) ~= length(ChanNames)
    Signal = permute(Signal,[2 1 3 4 5 6 7 8 9]);
end
Fs = SampleRate;

% 2022-05-24
Fcut_minfreq = Fs/10000;

averageepochduration = [];
splitnotch_factors = [];

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

% 2024-03-28: Basic warning if orientation is wrong
if ~iscell(Signal) && ~isstruct(Signal) && size(Signal,2) > size(Signal,1) && size(Signal,2) > Fs/4
    warning('The input signal %s has more columns (%i) than rows (%i). It is being interpreted as %.3g seconds long and %i channels. If this is not correct, transpose the matrix first.', inputname(1), size(Signal,2), size(Signal,1), size(Signal,1)/Fs, size(Signal,2));
end
if size(Signal,2) > 256
    warning('Plotting more than 256 channels (you have %i channels) can take a very long time. Recommend splitting the data into multiple plots.', size(Signal,2));
end



viewhand_ica_sig = ceil(rand*1000000000);
viewhand_ica_A = ceil(rand*1000000000);
viewhand_ica_W = ceil(rand*1000000000);
viewhand_psd = ceil(rand*1000000000);
viewhand_psd_axe = -1;
psd_held_hand = [];
psd_now_hand = [];
%viewhand_psd_peak = -1;
psd_held_lpxx = [];
psd_held_fxx = [];
psd_held_Time1 = NaN;
psd_held_Time2 = NaN;
psd_held_chancolor = [1 1 1];
psd_held_channame = '';
psd_held_filttext = '';
psd_held_linewidth = 2;
psd_held_colorfadefactor = 0.8;
psd_now_pxx = [];
psd_now_lpxx = [];
psd_now_fxx = [];
psd_now_Time1 = NaN;
psd_now_Time2 = NaN;
psd_now_chancolor = [1 1 1];
psd_now_channame = '';
psd_now_filttext = '';
psd_now_linewidth = 1;
selected_plothand = -1;
selected_timepoint = -1;
previously_selected_timepoint = -2;
FineSnapScale = 100;
selected_cursortype = 0;
CursorEnable = 0;
% LinemarkerEnable = 0;
InfoLabelEnable = 0;

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
chnc = nanmax(Signal) - nanmin(Signal) == 0; %#ok<*NANMIN,*NANMAX>

% Identify not-connected channels by zero-length channel names
if length(chnc) == length(ChanNames)
    chnc = chnc | cellfun(@isempty,ChanNames);
end

if exist('ica_W', 'var') && exist('ica_A', 'var') && ~isempty(ica_W) && ~isempty(ica_A)
    if size(ica_W,2) == size(Signal,2) && size(ica_A,1) == size(Signal,2)
        ica_W = ica_W(:,~chnc);
        ica_A = ica_A(~chnc,:);
    end
    if size(ltmat,2) == size(Signal,2) && size(ltmat,1) == size(Signal,2)
        ltmat = ltmat(~chnc,~chnc,:);
    end
else
    ica_W = [];
    ica_A = [];
end

if size(Signal,2) ~= numel(ChanNames)
    error('The numbers of channels in Signal and in ChanNames disagree.');
end

Signal = Signal(:,~chnc);
ChanNames = ChanNames(~chnc);

% Fix NaNs if any
Signal_nonfinite = ~isfinite(Signal);
Signal = fix_nans_for_filtering(Signal, bridgenans_method);



selchan = 1:length(ChanNames);
selica = [];


if isempty(ChanNames)
    error('There is nothing to plot. Either the input is an empty matrix or all channels are completely flat.');
end

fprintf('SignalViewer: Loading data... ');
try
    hashOpt.Method = 'MD5';
    hashOpt.Format = 'hex';
    signalHashStr = datahash(Signal,hashOpt);
catch
    signalHashStr = '';
end
fprintf('done.\n');
fprintf('              Allocating filter memories... ');

%[ChanNames,ix] = sort(ChanNames);
%Signal = Signal(:,ix);
%ica_W = ica_W(:,ix);
%ica_A = ica_A(ix,:);
Signal_postica = Signal;
Signal_postreref = Signal;
Signal_postnotch = Signal;
Signal_postbutter = Signal;
Signal_postenvelope = Signal;
Signal_postdecimation = Signal(:,1); % This is for each channel's rendering only
Signal_psd_source = Signal(:,1) * 0;

PerChannelFilterStates = false(4,length(ChanNames)); % Reref, Notch, Butter, Envelope
fprintf('done.\n');

fprintf('              Setting up plot elements... ');

% If a signal is known to be bandlimited due to low-pass and/or envelope
% filter, reduce the number of points to plot to speed up
SigBandwidth = Fs/2;
%BandLimitedInterpolation = 1;
BLIM = 1;

% If there are more signal points than horizontal screen resolution, reduce
% the number of points to plot to speed up
ScreenLimitedDownsampling = 0;
SLD_H = 32768;

DefaultYColor = [0.15 0.15 0.15];
BusyYColor = [1 0 0];
InactiveYColor = [0.65 0.65 0.65];


PlotHold = 0;


NotchFrequencies = 60;
MaxFilterOrder = 4;
FilterOrder = MaxFilterOrder;
%FilterChunkSec = 600;
FilterBusy = 0;
MovementBusy = 0;
ICA_Initialized = 0;
ica_sig = [];
icachans = {};

Nsch = length(ChanNames(selchan));
Ntp = size(Signal,1);

if ~SignalIsStitched
    splitnotch_factors = 1 - 1./(1+exp((1:Ntp).'-Ntp/2));
end

 
PermittedXZoomRanges   = [0.001 0.005 0.01 0.05 0.1 0.5 1 2 5 10 20 30 60 120 300 600 1200 1800 3600:3600:6*3600 8*3600 12*3600 24*3600 7*24*3600];
PermittedXZoomRanges = PermittedXZoomRanges(PermittedXZoomRanges > 1/Fs);
%PermittedChanSepRanges = [1e-6 2e-6 5e-6 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000];
PermittedChanSepRanges = [(1:9)*1e-6, (1:9)*1e-5, (1:9)*1e-4, (1:9)*1e-3, (1:9)*1e-2, (1:9)*1e-1, (1:9)*1e0, (1:9)*1e1, (1:9)*1e2, (1:9)*1e3, (1:9)*1e4, (1:9)*1e5, (1:9)*1e6, 1e7];
XTickSpacingsAndUnits = {
    1e-9    1e-9    'ns'
    2e-9    1e-9    'ns'
    10e-9   1e-9    'ns'
    20e-9   1e-9    'ns'
    100e-9  1e-9    'ns'
    200e-9  1e-9    'ns'
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
    60      1   's'
    300     1   's'
    600     1   's'
    1200    1   's'
    1800    1   's'
    3600    1   's'
    7200    1   's'
    3600*4  1   's'
    3600*8  1   's'
    3600*24 1   's'
    3600*24*7 1 's'
};
XTickSpacings = cat(1,XTickSpacingsAndUnits{:,1});

tmp = ceil(Nsch/3)*3;
Kolor = jet(tmp)*0.50;
tmp = reshape(reshape(rearrange_top_bottom((1:tmp).'),tmp/3,[]).',1,[]);
Kolor = Kolor(tmp,:);
if ~use_jet_colors
    try
        % 2024-10-08 Try to use the "lines" color map
        Kolor = lines(size(Signal,2));
    end
end
Nkolor = size(Kolor,1);

%EventKolor = [0.75 0.75 0];
%EventKolors = [0.75 0.75 0; 0 0.75 0.75; 0.75 0 0.75];
%Will generate EventKolors when the number of events is confirmed
EventPatchAlpha = 0.05;
EventFontName = 'Calibri';
EventFontSize = 16;
EventFontWeight = 'bold';
EventLineWidth = 0.1;
EventLineStyle = ':';

AxesFontName = 'Consolas';
AxesFontSize = 12;
CursorTextFontSize = 10;

PlotLineWidth = 0.5;

%InfoLabelFontName = 'Calibri';
InfoLabelFontName = 'Consolas';
%InfoLabelFontSizeScale = 0.8;
InfoLabelFontSizeScale = 2.0;
InfoLabelMargin = 0.01;


AxesPosition = [0.0300    0.0600    0.89    0.93];
ICA_AxesPosition = [0.0300    0.0600    0.96    0.93];

chansep = 100;

screensize = get(groot,'ScreenSize');
if isfield(opts, 'screensize') && length(opts.screensize) == 4
    screensize = opts.screensize;
end
tmp = mod(ceil(unixtimemillis()),604800000)+100000001;
while ishandle(tmp)
    tmp = tmp + 1;
end
if isempty(signalHashStr)
    signalHashStr = num2str(tmp);
end

fighand = figure(tmp);
if nargout > 0
    argout1 = fighand;
end
clf
set(fighand, 'ToolBar', 'none', 'MenuBar', 'none');
set(fighand, 'CloseRequestFcn', @f_main_close);


if SignalIsStitched
    tmp = sprintf('Signal is stitched from %g epochs with %.2g s average length.', size(EventTimePoints,1), averageepochduration);
else
    tmp = '';
end
PSigName = '';
if ~isempty(inputname(1))
    PSigName = [' «' inputname(1) '»'];
end
if ~isempty(FooterMessage)
    PSigName = [PSigName ' ' FooterMessage];
end
set(fighand, 'Name', [PSigName ' ' tmp ' Hash ' signalHashStr '.']);
hold on
set(fighand,'Position',screensize);
YLim = [-chansep*Nsch-0.5*chansep, -chansep+0.5*chansep];
XLim = [0 10];

%changed_CND = 0;
set(gca, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)), 'YLim', YLim, 'XLim', XLim, 'Position', AxesPosition, 'FontWeight', 'bold', 'FontName', AxesFontName, 'FontUnits', 'pixel', 'FontSize', AxesFontSize); %#ok<*NBRAK>
axehand = get(fighand,'CurrentAxes');

t1 = 1;
t2 = 2;
Time = (0:Ntp-1)/Fs;
Time_min = 0;
Time_max = (Ntp-1)/Fs;

PipMarkerSize = 8;

cursorlinehand = plot([0 0], [0 0], '-', 'LineWidth', 3, 'Color', [0.8 0.8 0.8]);
for ch = Nsch:-1:1
    plotpip1hand(ch) = plot(0, 0, '+', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', PipMarkerSize, 'LineWidth', 2, 'Visible', 'off');
    plotpip2hand(ch) = plot(0, 0, 'v', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', PipMarkerSize, 'LineWidth', 2, 'Visible', 'off');
    plotpip3hand(ch) = plot(0, 0, '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', PipMarkerSize, 'LineWidth', 2, 'Visible', 'off');
    plothand(ch) = plot(Time(t1:t2), Signal(t1:t2,selchan(ch)) - nanmean(Signal(t1:t2,selchan(ch))) - chansep*ch); %#ok<*NANMEAN>
    set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:), 'LineWidth', PlotLineWidth);
    set(plothand(ch), 'ButtonDownFcn', @f_plothand_buttondown);
    setappdata(plothand(ch), 'chanind', selchan(ch));
    setappdata(plothand(ch), 'channame', ChanNames{selchan(ch)});
    if ismember(ChanNames{selchan(ch)}, nofiltchannames)
        set(plothand(selchan(ch)), 'Color', [0 0 0]);
    end
end
for ch = Nsch:-1:1
    plottext_cursor_hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', CursorTextFontSize);
    plottext_lmin_hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', CursorTextFontSize);
    plottext_lmax_hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', CursorTextFontSize);
end
for ch = Nsch:-1:1
    plotinfolabel1hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Color', (Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', AxesFontSize*InfoLabelFontSizeScale, 'FontName', InfoLabelFontName, 'FontWeight', 'bold', 'Margin', InfoLabelMargin);
    plotinfolabel2hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Color', (Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', AxesFontSize*InfoLabelFontSizeScale, 'FontName', InfoLabelFontName, 'Margin', InfoLabelMargin);
    plotinfolabel3hand(ch) = text(0, 0, '', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', lightercolor(Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Color', (Kolor(mod(selchan(ch)-1,Nkolor)+1,:)), 'Visible', 'off', 'FontUnits', 'pixel', 'FontSize', AxesFontSize*InfoLabelFontSizeScale, 'FontName', InfoLabelFontName, 'Margin', InfoLabelMargin);
    set(plotinfolabel1hand(ch), 'FontUnits', 'normalized');
    set(plotinfolabel2hand(ch), 'FontUnits', 'normalized');
    set(plotinfolabel3hand(ch), 'FontUnits', 'normalized');
    if ismember(ChanNames{selchan(ch)}, nofiltchannames)
        set(plotinfolabel1hand(selchan(ch)), 'BackgroundColor', [0.9 0.9 0.9], 'Color', [0 0 0]);
    end
end
if ~isempty(plothand)
    selected_plothand = plothand(1);
end

set(axehand, 'FontUnits', 'normalized');
AxesFontSize = get(axehand, 'FontSize');

fprintf('done.\n');
fprintf('              Labeling events... ');

%tmax = max(Time);
%SigChunkRendered = false(ceil(tmax / FilterChunkSec),length(ChanNames));
%ChunkIndexMax = size(SigChunkRendered,1);

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
        for i = 1:size(tmp,1)
            EventTimeStamps(i,:) = {[tmp(i,1) tmp(i,2)], sprintf('Event %g', tmp(i,3))};
        end
    elseif isnumeric(EventTimeStamps) && size(EventTimeStamps,2) == 2
        % boundary format without augmented label
        tmp = EventTimeStamps;
        EventTimeStamps = {};
        for i = 1:size(tmp,1)
            EventTimeStamps(i,:) = {[tmp(i,1) tmp(i,2)], sprintf('Event %i', i)};
        end
    elseif isnumeric(EventTimeStamps) && size(EventTimeStamps,2) == 1
        % simple time stamps with no labels
        tmp = EventTimeStamps;
        EventTimeStamps = {};
        %k = 0;
        for i = 1:size(tmp,1)
            %k = k + 1;
            EventTimeStamps(i,:) = {tmp(i,1), sprintf('Event %i', i)};
        end
    end
    
    
end

if exist('EventTimeStamps','var') && ~isempty(EventTimeStamps) && iscell(EventTimeStamps) && size(EventTimeStamps,2) == 2
    EventEnable = 1;
else
    EventEnable = 0;
end


if EventEnable
    Nevents = size(EventTimeStamps,1);
    EventTimes = nan(Nevents,2);
    tmp2 = ceil(Nevents/2)*2;
    EventKolors = jet(tmp2);
    EventLineWidths = repmat(EventLineWidth, [Nevents 1]);
    EventLineStyles = repmat({EventLineStyle}, [Nevents 1]);
    tmp = 1:tmp2;
    tmp = reshape([tmp(tmp2/2+1:end); tmp(1:tmp2/2)], 1, []);
    EventKolors = EventKolors(tmp,:);

    %Po240712: Make them darker so easier to see
    EventKolors = EventKolors * 2/3;

    %Po250505
    if isfield(opts, 'EventFontSize') && isscalar(opts.EventFontSize) && isnumeric(opts.EventFontSize) && opts.EventFontSize > 0 && isfinite(opts.EventFontSize)
        EventFontSize = opts.EventFontSize;
    end

    if isfield(opts, 'EventColors') && isnumeric(opts.EventColors) && isequal(size(EventKolors),size(opts.EventColors)) && isreal(opts.EventColors) && min(min(opts.EventColors)) >= 0
        if max(max(opts.EventColors)) <= 1
            EventKolors = opts.EventColors;
        else
            EventKolors = opts.EventColors ./ max(max(opts.EventColors));
        end
    elseif isfield(opts, 'EventColorRules') && iscell(opts.EventColorRules) && size(opts.EventColorRules,2) == 2
        for i = 1:size(opts.EventColorRules,1)
            foundind = find(~cellfun(@isempty,regexp(EventTimeStamps(:,2), opts.EventColorRules{i,1}, 'match', 'once')));
            for j = 1:length(foundind)
                EventKolors(foundind(j),:) = opts.EventColorRules{i,2};
                if max(max(EventKolors(foundind(j),:))) > 1
                    EventKolors(foundind(j),:) = EventKolors(foundind(j),:) ./ max(max(EventKolors(foundind(j),:)));
                end
            end
        end
    end

    if isfield(opts, 'EventLineRules') && iscell(opts.EventLineRules) && size(opts.EventLineRules,2) == 3
        for i = 1:size(opts.EventLineRules,1)
            foundind = find(~cellfun(@isempty,regexp(EventTimeStamps(:,2), opts.EventLineRules{i,1}, 'match', 'once')));
            for j = 1:length(foundind)
                EventLineWidths(foundind(j),1) = opts.EventLineRules{i,2};
                EventLineStyles(foundind(j),1) = opts.EventLineRules(i,3);
            end
        end
    end

    for i = 1:Nevents
        tmp = EventTimeStamps{i,1};
        EventTimes(i,:) = [tmp(1) tmp(end)];
        if EventTimes(i,1) > EventTimes(i,2)
            EventTimes(i,2) = EventTimes(i,1);
        end
    end

    [EventTimes,I] = sortrows(EventTimes);
    EventTimeStamps = EventTimeStamps(I,:);
    EventKolors = EventKolors(I,:);
    EventLineWidths = EventLineWidths(I,:);
    EventLineStyles = EventLineStyles(I,:);

    % try
    %     EventTimeStamps = sortrows(EventTimeStamps);
    % catch
    %     [B,I] = sortrows(EventTimeStamps(:,1));
    %     EventTimeStamps = [B EventTimeStamps(I,2)];
    %     clear B I
    % end
    % EventTimes = cell2mat(EventTimeStamps(:,1));
    % Initial y positions for event text labels (redraw will overwrite the positions later)
    YPos = sort(mean(YLim)+diff(YLim)/20*([-8:2:8]), 'descend');
    NYPos = length(YPos);
    imax = size(EventTimeStamps,1);
    ndigits = floor(log10(imax))+1;
    fprintf(['%' num2str(ndigits) 'i/%' num2str(ndigits) 'i '], 0, imax);
    tlastreported = now;
    for i = 1:imax
        %eventplothand(i) = plot( EventTimeStamps{i,1}*[1 1], [-10000000*(Nsch+1), 10000000], '-', 'Color', EventKolor ); %#ok<NASGU>
        if EventTimes(i,1) == EventTimes(i,2)
            EventKolor = EventKolors(mod(i-1,size(EventKolors,1))+1,:);
            %EventKolor = EventKolor/2+0.50;
            eventplothand(i) = plot( EventTimes(i,1)*[1 1], [-10000000*(Nsch+1), 10000000], 'LineStyle', EventLineStyles{i}, 'Color', EventKolor, 'LineWidth', EventLineWidths(i) );
        else
            EventKolor = EventKolors(mod(i-1,size(EventKolors,1))+1,:);
            eventplothand(i) = patch( [EventTimes(i,1) EventTimes(i,2) EventTimes(i,2) EventTimes(i,1)], [-10000000*(Nsch+1), -10000000*(Nsch+1), 10000000, 10000000], '-', 'FaceColor', EventKolor, 'EdgeColor', 'none', 'FaceAlpha', EventPatchAlpha); 
        end
        larrow = '«';
        rarrow = '';
        horali = 'left';
        if EventTimes(i,1) - XLim(1) > (XLim(2) - XLim(1))*0.5
            larrow = '';
            rarrow = '»';
            horali = 'right';
        end

        %Po240712: Changed event text color and background color
        %EventTextBackgroundKolor = EventKolor;
        EventTextBackgroundKolor = 'none';
        EventTextKolor = EventKolor;
        %if mean(rgb2gray(EventKolor)) < 0.5
        %    EventTextKolor = [1 1 1];
        %else
        %    EventTextKolor = [0 0 0];
        %end
        eventtexthand(i) = text( EventTimes(i,1), YPos(mod(i-1,NYPos)+1), [larrow EventTimeStamps{i,2} rarrow], 'FontName', EventFontName, 'FontUnits', 'pixel', 'FontSize', EventFontSize, 'FontWeight', EventFontWeight, 'HorizontalAlignment', horali, 'BackgroundColor', EventTextBackgroundKolor, 'Color', EventTextKolor);
        if now - tlastreported > 1/86400
            fprintf([repmat('\b',[1 ndigits*2+2]) '%' num2str(ndigits) 'i/%' num2str(ndigits) 'i '], i, imax);
            tlastreported = now;
        end
    end
    set(eventtexthand(i), 'Visible', 'off');
    fprintf([repmat('\b',[1 ndigits*2+2]) '%' num2str(ndigits) 'i/%' num2str(ndigits) 'i '], i, imax);
end


%2024-10-03: Move patches to background so that user can still click on traces for PSD.
o_patches = findobj(get(axehand,'Children'),'Type','Patch');
set(axehand,'Children',[setdiff(get(axehand,'Children'), o_patches); o_patches]);

fprintf('done.\n');
fprintf('              Setting up buttons... ');

BandPassFilter.state = 0;
BandPassFilter.cutoff = [0 Fs/2];
%HighPassFilter.state = 0;
%HighPassFilter.cutoff = ceil(Fcut_minfreq*101)/100;
%LowPassFilter.state = 0;
%LowPassFilter.cutoff = floor(Fs/2-1);
EnvelopeFilter.state = 0;
EnvelopeFilter.cutoff = 2;
NotchFilter.state = 0;
RerefFilter.state = 0;
RerefFilter.type = 'car';
RerefFilter.chanidx = [];
NotchFilter.order = 0; % Must be even (0 disables)
NotchFilter.qfactor = 16;
NotchFilter.pendingqfactor = 16;
notchordermax = 10;
notchqfactormax = 999;
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
fontcolor_on3 = [1 0.25 0];
fontcolor_off = [0 0 0];
backgroundcolor_buttondefault = [0.7 0.7 0.7];
backgroundcolor_buttonalert = [1 1 0];


h_hugetext = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.35 0.5 0.25 0.16], 'String', 'BIG TEXT', 'FontUnits', 'normalized', 'FontSize', 0.4, 'Visible', 'off', 'BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1]);
h_bigtext = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.25 0.5 0.45 0.05], 'String', 'BIG TEXT', 'FontUnits', 'normalized', 'FontSize', 0.4, 'Visible', 'off');

NormalizedControlFontSize = 0.8;

h_xoomtitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.96 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'Horizontal Zoom', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_xzoomout = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.94 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_xzoomlevel = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.94 0.94 0.04 0.02], 'BackgroundColor', [0.7 0.7 0.9], 'String', '10 s', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_xzoomin = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.94 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_yzoomtitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.91 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', '# chans on screen', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_yzoomout = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.89 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_yzoomlevel = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.94 0.89 0.04 0.02], 'BackgroundColor', [0.7 0.7 0.9], 'String', '64 ch', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_yzoomin = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.89 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_sensititle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.86 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'Y Sensitivity', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_sepup = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.84 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '-', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_sensitivity = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.94 0.84 0.04 0.02], 'BackgroundColor', [0.7 0.7 0.9], 'String', '100 uV', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_sepdown = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98 0.84 0.015 0.02], 'BackgroundColor', [0.7 0.7 0.7], 'String', '+', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_pantitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.81 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'Pan', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>
h_panleft = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.955 0.775 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '<', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_panright = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.975 0.775 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '>', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_panup = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.80 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', '^', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_pandown = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.775 0.020 0.025], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'v', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_hold_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.75 0.070 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Pause Plotting', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8);
%h_hold_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.97 0.74 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_verticaloverlapdisallow_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.7315 0.070 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Reduce Overlap', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8);

h_psd_plot = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.70 0.03 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'PSD', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8);
h_autofit = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.965 0.70 0.03 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Auto Fit', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.5);


h_passfilttitle = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.68 0.07 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Filters', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_ica_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.66 0.025 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ICA', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Tooltip', 'Use the bottom right panel to control Independent Component Analysis.');
h_car_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.95 0.66 0.025 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'CAR', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Tooltip', 'The currently selected channels are re-referenced and will not change until you click this button again.');
h_car_chancount = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.975 0.66 0.020 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', '0ch', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8);

h_notch_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.64 0.025 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Notch', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Tooltip', 'Notch filter turns on automatically after entering a valid set of order and Q factor. Try entering order 4 and Q factor 16.');
%h_notch_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.942 0.64 0.017 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_notch_order = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.95 0.64 0.010 0.017], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(NotchFilter.order), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_notch_orderunit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.96 0.64 0.019 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'order', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_notch_qfactor = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.98 0.64 0.010 0.017], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(NotchFilter.qfactor), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_notch_qfactorunit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.99 0.64 0.005 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Q', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>

h_bpf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.62 0.025 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Butter', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Tooltip', 'Bandpass/Lowpass/Highpass Butterworth filter turns on automatically after entering a pair of valid cutoff frequencies.');
h_bpf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.95 0.62 0.030 0.017], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2strcompact(BandPassFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_bpf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.98 0.62 0.015 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

%h_hpf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.64 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'HighPass', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_hpf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.96 0.64 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(HighPassFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_hpf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.98 0.64 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

%h_lpf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.62 0.035 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'LowPass', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_lpf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.96 0.62 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(LowPassFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_lpf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.98 0.62 0.015 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_evf_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.60 0.035 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Envelope', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_evf_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.94 0.60 0.02 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_evf_cutoff = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.96 0.60 0.020 0.017], 'BackgroundColor', [0.99 0.99 0.99], 'String', num2str(EnvelopeFilter.cutoff), 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_evf_unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.98 0.60 0.015 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hz', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize); %#ok<NASGU>

h_zscore_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.58 0.030 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Zscore', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
%h_zscore_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.945 0.58 0.015 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_off, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_fastdraw_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.965 0.58 0.030 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Fast', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Value', 0, 'FontWeight', fontweight_off, 'ForegroundColor', fontcolor_off);
%h_fastdraw_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.985 0.58 0.015 0.017], 'BackgroundColor', [0.7 0.7 0.7], 'String', text_on, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_chansel_warnconfirm = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.56 0.07 0.12], 'ForegroundColor', [1 1 0], 'BackgroundColor', [0.2, 0.2, 0.2], 'String', 'Apply or revert channel selection first.', 'FontUnits', 'normalized', 'FontSize', 0.15, 'Visible', 'off');


h_chansel_title = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.525 0.030 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Plotted Channels', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.5);
h_chansel_list = uicontrol(fighand, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.925 0.14 0.040, 0.380], 'FontUnits', 'pixel');
h_chansel_commandentry = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.124 0.039 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Adv Select', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Tooltip', 'Enter regular expression patterns to step 1: select channels and then step 2: deselect channels. Leave blank to skip a step.');
h_chansel_confirm = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.108 0.02 0.015], 'BackgroundColor', backgroundcolor_buttondefault, 'String', 'Apply', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
%h_chansel_auto = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.945 0.108 0.02 0.015], 'BackgroundColor', backgroundcolor_buttondefault, 'String', 'Auto', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_chansel_reset = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.945 0.108 0.02 0.015], 'BackgroundColor', backgroundcolor_buttondefault, 'String', 'Revert', 'FontUnits', 'normalized', 'Visible', 'on');

h_infolabel_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.925 0.091 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'I-box', 'FontUnits', 'normalized', 'FontName', 'Times New Roman');
h_cursor_switch = uicontrol(fighand, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [0.945 0.091 0.02 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Cursor', 'FontUnits', 'normalized', 'FontName', 'Times New Roman');
%h_cursor_state = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.955 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'off', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

h_eventlabels_showhide = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.075 0.039 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Hide Events', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);

%h_linemarker_switch = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.06 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'X', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);

h_signal_export = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.935 0.06 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Export Sig', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);


h_icasel_title = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.965 0.525 0.030 0.030], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'ICA Comps', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.5);
h_icasel_list = uicontrol(fighand, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0.965 0.14 0.030, 0.380], 'FontUnits', 'pixel');
h_icasel_confirm = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.966 0.10 0.019 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Start', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_reset = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.985 0.10 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'R', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_sources = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.965 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'S', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_mixmat = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.975 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_sepmat = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.985 0.08 0.01 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_icasel_view_export = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.965 0.06 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Export ICA', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);

h_axesfont_inc = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.925 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'A', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_axesfont_dec = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.93375 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'a', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_linewidth_inc = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.9425 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'L', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_linewidth_dec = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.95125 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'l', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_eventfont_inc = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'E', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_eventfont_dec = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.96785 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'e', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_windowhsize_inc = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.9775 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'W', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);
h_windowhsize_dec = uicontrol(fighand, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.98625 0.04 0.00875 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'w', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9);

h_xspan_text = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.925 0.000 0.07 0.035], 'BackgroundColor', [0.7 0.7 0.9], 'String', 'full t range [0,12345]', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.33);

h_xspan_edittext1intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.020 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'XLim(1) =', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_xspan_edit1 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.04 0.020 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_xspan_samplenumber1intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.010 0.03 0.010], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'sample#', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_xspan_samplenumber1 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.04 0.010 0.03 0.010], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Enable', 'off');
h_xspan_edittext1unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.07 0.020 0.004 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 's', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>

h_hpanamount_intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.075 0.020 0.048 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Arrow key step size: ', 'HorizontalAlignment', 'left', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8); %#ok<NASGU>
h_hpanamount_edit = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.125 0.020 0.020 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '100%', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.8); 

h_lmoddate = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.0 0.000 0.095 0.009], 'BackgroundColor', [0.7 0.7 0.7], 'String', ['SignalViewer version: ' lmoddate], 'FontUnits', 'normalized', 'FontSize', 0.8, 'FontName', 'Calibri', 'HorizontalAlignment', 'left'); %#ok<NASGU>
h_hintbar = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.15 0.020 0.705 0.015], 'BackgroundColor', [0.7 0.9 0.7], 'String', '', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'FontName', 'Verdana', 'HorizontalAlignment', 'left');

h_xspan_edittext2intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.856 0.020 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'XLim(2) =', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_xspan_edit2 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.886 0.020 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize);
h_xspan_samplenumber2intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.856 0.010 0.03 0.01], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'sample#', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>
h_xspan_samplenumber2 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.886 0.010 0.03 0.01], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'Enable', 'off');
h_xspan_edittext2unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.916 0.020 0.004 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 's', 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize*0.9); %#ok<NASGU>

% h_yspan_edittext1intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.41 0.025 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Selected channel YMin =', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edit1 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.44 0.025 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edittext1unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.47 0.025 0.025 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'µV', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% 
% h_yspan_edittext2intro = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.50 0.025 0.03 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Selected channel YMax =', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edit2 = uicontrol(fighand, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.53 0.025 0.03 0.015], 'BackgroundColor', [0.99 0.99 0.99], 'String', '00000', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);
% h_yspan_edittext2unit = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.56 0.025 0.025 0.015], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'µV', 'FontUnits', 'normalized', 'FontSize', ControlFontSize);


h_footer_message = uicontrol(fighand, 'Style', 'text', 'Units', 'normalized', 'Position', [0.10 0.000 0.75 0.018], 'BackgroundColor', [0.8 0.8 1.0], 'String', FooterMessage, 'FontUnits', 'normalized', 'FontSize', NormalizedControlFontSize, 'FontName', 'Verdana');
if isempty(FooterMessage)
    set(h_footer_message, 'Visible', 'off');
end

if ~EventEnable
    set([h_eventfont_inc h_eventfont_dec h_eventlabels_showhide], 'Enable', 'off', 'Visible', 'off');
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

set(h_hintbar, 'Enable', 'Inactive', 'ButtonDownFcn', @f_hintbar);
set(h_xspan_edit1, 'Callback', @f_xspan_edit1);
set(h_xspan_edit2, 'Callback', @f_xspan_edit2);
set(h_hpanamount_edit, 'Callback', @f_hpanamount_edit);

set(h_hold_switch, 'Callback', @f_hold_switch);
set(h_verticaloverlapdisallow_switch, 'Callback', @f_verticaloverlapdisallow_switch);

set(h_bpf_switch, 'Callback', @f_bpf_switch);
set(h_bpf_cutoff, 'Callback', @f_bpf_cutoff);
set(h_ica_switch, 'Callback', @f_ica_switch);
set(h_car_switch, 'Callback', @f_car_switch);
%set(h_hpf_switch, 'Callback', @f_hpf_switch);
%set(h_lpf_switch, 'Callback', @f_lpf_switch);
%set(h_hpf_cutoff, 'Callback', @f_hpf_cutoff);
%set(h_lpf_cutoff, 'Callback', @f_lpf_cutoff);
set(h_evf_switch, 'Callback', @f_evf_switch);
set(h_evf_cutoff, 'Callback', @f_evf_cutoff);
set(h_notch_switch, 'Callback', @f_notch_switch);
set(h_notch_order, 'Callback', @f_notch_order);
set(h_notch_qfactor, 'Callback', @f_notch_qfactor);
set(h_zscore_switch, 'Callback', @f_zscore_switch);
set(h_fastdraw_switch, 'Callback', @f_fastdraw_switch);
set(h_chansel_list, 'Callback', @f_chansel_list);
set(h_chansel_commandentry, 'Callback', @f_chansel_commandentry);
set(h_eventlabels_showhide, 'Callback', @f_eventlabels_showhide);
set(h_chansel_confirm, 'Callback', @f_chansel_confirm);
set(h_chansel_reset, 'Callback', @f_chansel_reset);
set(h_psd_plot, 'Callback', @f_psd_plot);
set(h_autofit, 'Callback', @f_autofit);
set(h_icasel_list, 'Callback', @f_icasel_list);
set(h_icasel_confirm, 'Callback', @f_icasel_confirm);
set(h_icasel_reset, 'Callback', @f_icasel_reset);
set(h_icasel_view_sources, 'Callback', @f_icasel_view_sources);
set(h_icasel_view_mixmat, 'Callback', @f_icasel_view_mixmat);
set(h_icasel_view_sepmat, 'Callback', @f_icasel_view_sepmat);
set(h_icasel_view_export, 'Callback', @f_icasel_view_export);
%set(h_chansel_auto, 'Callback', @f_chansel_auto);
set(h_infolabel_switch, 'Callback', @f_infolabel_enable);
set(h_cursor_switch, 'Callback', @f_cursor_enable);
%set(h_linemarker_switch, 'Callback', @f_linemarker_switch);
set(h_signal_export, 'Callback', @f_signal_export);
set(h_axesfont_inc, 'Callback', @f_axesfont_inc);
set(h_axesfont_dec, 'Callback', @f_axesfont_dec);
set(h_linewidth_inc, 'Callback', @f_linewidth_inc);
set(h_linewidth_dec, 'Callback', @f_linewidth_dec);
set(h_eventfont_inc, 'Callback', @f_eventfont_inc);
set(h_eventfont_dec, 'Callback', @f_eventfont_dec);
set(h_windowhsize_inc, 'Callback', @f_windowhsize_inc);
set(h_windowhsize_dec, 'Callback', @f_windowhsize_dec);


for ip01 = 1:length(plottext_lmin_hand)
    set(plottext_lmin_hand(ip01), 'ButtonDownFcn', @f_localminlabel_savepoint);
end
for ip01 = 1:length(plottext_lmax_hand)
    set(plottext_lmax_hand(ip01), 'ButtonDownFcn', @f_localmaxlabel_savepoint);
end

set(axehand, 'ButtonDownFcn', @f_axe_buttondown)
set(fighand, 'KeyPressFcn', @f_fig_keypress);
set(fighand, 'WindowScrollWheelFcn', @f_fig_scrollwheel);
set(fighand,'Position',screensize);



reref_update();
notch_update();
bandpass_update();
envelope_update();
render_update();
%set(h_bigtext, 'Visible', 'off', 'String', '');

%set(h_xspan_text, 'String', ['full t range [' num2str(round(min(Time))) ', ' num2str(round(max(Time))) '] s']);
set(h_xspan_text, 'String', sprintf('full t range [%i, %i] s\nSample rate %g Hz\n%i chans, %i epochs', round(min(Time)), round(max(Time)), Fs, length(ChanNames), size(EventTimePoints,1)));

if isfield(opts, 'screensize') && length(opts.screensize) == 4
else
    try
        pause(0.00001);
        set(fighand, 'WindowState', 'maximized');
    catch
        try
            pause(0.00001);
            oldWarningState = warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');
            frame_h = get(handle(fighand),'JavaFrame'); %#ok<JAVFM>
            set(frame_h,'Maximized',1);
            warning(oldWarningState);
        end
    end
end

%Po240614: Hold until startup is done
f_hold_switch(-10000, []);

if fastdraw
    ScreenLimitedDownsampling = 1;
    set(h_fastdraw_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
end

if ~isempty(set_selectchannames_to)
    selchan = chan2idx(ChanNames, set_selectchannames_to, 1);
    set(h_chansel_list, 'Value', selchan);
    f_chansel_confirm([], []);
end

if isfield(opts, 'car')
    if isscalar(opts.car) && opts.car
        f_car_switch([], []);
    end
end

if isfield(opts, 'notch')
    if ~isempty(opts.notch) && isnumeric(opts.notch)
        NotchFrequencies = opts.notch(:).';
        NotchFrequencies(~(isfinite(NotchFrequencies) & NotchFrequencies > 0 & mod(NotchFrequencies,Fs/2) > 0)) = NaN;
        NotchFrequencies = NotchFrequencies(~isnan(NotchFrequencies));
        %f_notch_switch(fighand, []);
        set(h_notch_order, 'String', '4');
        f_notch_order(fighand, []);
    end
end
hasvalidbpf = 0;
if ~isfield(opts, 'bandpass') || numel(opts.bandpass) ~= 2 || ~isnumeric(opts.bandpass)
    opts.bandpass = [0 Fs/2];
end

if isfield(opts, 'highpass')
    if isscalar(opts.highpass) && opts.highpass(1) > 0 && opts.highpass(1) < Fs/2
        opts.bandpass(1) = opts.highpass;
        hasvalidbpf = 1;
    end
end
if isfield(opts, 'lowpass')
    if isscalar(opts.lowpass) && opts.lowpass(1) > 0 && opts.lowpass(1) < Fs/2
        opts.bandpass(2) = opts.lowpass;
        hasvalidbpf = 1;
    end
end
if isfield(opts, 'bandpass')
    if numel(opts.bandpass) == 2 && (opts.bandpass(1) > 0 || opts.bandpass(2) < Fs/2)
        hasvalidbpf = 1;
    end
end
if hasvalidbpf
    set(h_bpf_cutoff, 'String', num2strcompact(opts.bandpass));
    f_bpf_cutoff([], []);
end

if isfield(opts, 'envelope')
    if isscalar(opts.envelope) && isnumeric(opts.envelope) && opts.envelope < Fs/2
        set(h_evf_cutoff, 'String', num2str(opts.envelope));
        f_evf_cutoff([], []);
        f_evf_switch([], []);
    end
end

if isfield(opts, 'zscore')
    if isscalar(opts.zscore) && opts.zscore
        f_zscore_switch(-10000, []);
    end
end

%drawnow

set(fighand, 'HandleVisibility', 'callback');

if set_xlim_to(1) ~= XLim(1) || set_xlim_to(2) ~= XLim(2)
    if isfinite(set_xlim_to(1)) && isfinite(set_xlim_to(2))
        XLim(1) = set_xlim_to(1);
        XLim(2) = set_xlim_to(2);
    else
        XLim(1) = 0;
        XLim(2) = size(Signal,1)/SampleRate;
    end
    set(axehand, 'XLim', XLim);
    MovementBusy = 1;
    resnap_pan();
    MovementBusy = 0;
    set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
end

if set_chansep_to >= 0
    autofit();
end
if set_chansep_to > 0
    refit(set_chansep_to);
end


% Hide disabled buttons
if disable_ica
    set([h_icasel_title h_icasel_list h_icasel_confirm h_icasel_reset h_icasel_view_sources h_icasel_view_mixmat h_icasel_view_sepmat h_icasel_view_export], 'Visible', 'off');
elseif ~isempty(ica_W)
    %2024-10-10 If ICA matrices are in the options, auto-populate ICA.
    f_icasel_confirm([], 'startup');
end

% Done with startup

fprintf('done.\n');
fprintf('              Done with startup: Figure %i.\n', get(fighand, 'Number'));


f_hold_switch(-100000, []);
figure(fighand);

if isfield(opts,'saveas') && ischar(opts.saveas) && ~isfile(opts.saveas) && length(opts.saveas) > 4 && strcmpi(opts.saveas(end-3:end),'.png')
    try
        saveas(fighand, opts.saveas);
    catch
        warning('Unable to save into %s', opts.saveas);
    end
end


    function f_fig_scrollwheel(hObject, eventdata)
        if eventdata.VerticalScrollCount < 0
            % Mouse whell scroll up

        elseif eventdata.VerticalScrollCount > 0
            % Mouse whell scroll down
        end

    end


    function f_localminlabel_savepoint(hObject, eventdata)
        if isequal(eventdata.Button, 1)
            if size(SavedPointsTable,2) >= 2 && size(SavedPointsTable,1) >= 1
                tcpair = cell2mat(SavedPointsTable(:,1:2));
                L = getappdata(hObject, 'tvalue')==tcpair(:,1) & getappdata(hObject, 'chanind')==tcpair(:,2);
            else
                L = false;
            end
            if any(L)
                % Already exists. Clicking removes from saved list
                SavedPointsTable = SavedPointsTable(~L,:);
                % autosave the table
                assignin('base',['signalviewer_SavedPointsTable_' signalHashStr],SavedPointsTable);
                % update the label
                set(hObject, 'String', getappdata(hObject,'originalString'));
            else
                SavedPointsTable = [SavedPointsTable ; [{getappdata(hObject, 'tvalue')} ...
                    {getappdata(hObject, 'chanind')} ...
                    {getappdata(hObject, 'ChanName')} ...
                    {getappdata(hObject, 'yvalue')} ...
                    {2}]
                    ];
                % autosave the table
                assignin('base',['signalviewer_SavedPointsTable_' signalHashStr],SavedPointsTable);
                % update the label
                set(hObject, 'String', sprintf('%s\n%s',getappdata(hObject,'originalString'),'(SAVED)'));
            end
        end
    end

    function f_localmaxlabel_savepoint(hObject, eventdata)
        if isequal(eventdata.Button, 1)
            if size(SavedPointsTable,2) >= 2 && size(SavedPointsTable,1) >= 1
                tcpair = cell2mat(SavedPointsTable(:,1:2));
                L = getappdata(hObject, 'tvalue')==tcpair(:,1) & getappdata(hObject, 'chanind')==tcpair(:,2);
            else
                L = false;
            end
            if any(L)
                % Already exists. Clicking removes from saved list
                SavedPointsTable = SavedPointsTable(~L,:);
                % autosave the table
                assignin('base',['signalviewer_SavedPointsTable_' signalHashStr],SavedPointsTable);
                % update the label
                set(hObject, 'String', getappdata(hObject,'originalString'));
            else
                SavedPointsTable = [SavedPointsTable ;
                    [{getappdata(hObject, 'tvalue')} ...
                    {getappdata(hObject, 'chanind')} ...
                    {getappdata(hObject, 'ChanName')} ...
                    {getappdata(hObject, 'yvalue')} ...
                    {3} ]
                    ];
                % autosave the table
                assignin('base',['signalviewer_SavedPointsTable_' signalHashStr],SavedPointsTable);
                % update the label
                set(hObject, 'String', sprintf('%s\n%s',getappdata(hObject,'originalString'),'(SAVED)'));
            end
        end
    end


    function f_axe_buttondown(hObject, eventdata)
        if isequal(eventdata.Button, 2)
            % Middle mouse click on the plot area
            autofit();
        end
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
                set(h_hintbar, 'String', '');
                %uin = inputdlg('Center on specific time in seconds, specify a time range, or enter "all":');
                uin = inputdlg('Center on specific time in seconds, specify a time range, or enter "all":', '', 1, {sprintf('%g   %g', XLim(1), XLim(2))});
                gotosec = str2double(uin);
                if isfinite(gotosec)
                    if FilterBusy
                        return;
                    end
                    XRange = XLim(2)-XLim(1);
                    XLim = gotosec + [-1 1]*XRange/2;
                    set(axehand, 'XLim', XLim);
                    resnap_pan();
                    set(h_hintbar, 'String', ['Zoomed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
                else
                    % Maybe the user entered two numbers?
                    if ~isempty(uin)
                        if strcmpi(uin{1}, 'all')
                            uin{1} = sprintf('%g %g', floor(min(Time)), ceil(max(Time)));
                        end
                        gotosec = sscanf(uin{1}, '%g %g');
                    end
                    if numel(gotosec) == 2 && isfinite(gotosec(1)) && isfinite(gotosec(2))
                        XLim(1:2) = gotosec(1:2);
                        set(axehand, 'XLim', XLim);
                        resnap_pan();
                        set(h_hintbar, 'String', ['Zoomed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
                    end
                end
            case 'c'
                %2024-11-06: Center on channel for PSD
                set(h_hintbar, 'String', '');
                tmp_answer = inputdlg('Choose a channel for PSD:', 'PSD chan', [1 40]);
                
                if isempty(tmp_answer) || isempty(tmp_answer{1})
                    return
                end
                choose_psd_channel(tmp_answer{1});

            case 'l'
                %2024-11-19: New feature to switch to log10 frequency scale
                %(only works in the PSD window)
                if isequal(hObject, viewhand_psd)
                    if strcmpi(psd_xscale, 'linear')
                        psd_xscale = 'log';
                    else
                        psd_xscale = 'linear';
                    end
                    update_psd();
                end

            case 'h'
                %2024-10-09: New feature to hold the PSD (only works in the PSD window)
                if isequal(hObject, viewhand_psd)
                    if ~isempty(psd_now_lpxx) && length(psd_now_fxx) == length(psd_now_lpxx)
                        if isequal(psd_held_lpxx,psd_now_lpxx)
                            % Identical. Clear hold.
                            psd_held_lpxx = [];
                            psd_held_fxx = [];
                            psd_held_Time1 = NaN;
                            psd_held_Time2 = NaN;
                            psd_held_chancolor = [1 1 1];
                            psd_held_channame = '';
                            psd_held_filttext = '';
                            set(h_hintbar, 'String', 'PSD plot hold is released.');
                        else
                        psd_held_lpxx = psd_now_lpxx;
                        psd_held_fxx = psd_now_fxx;
                        psd_held_Time1 = psd_now_Time1;
                        psd_held_Time2 = psd_now_Time2;
                        psd_held_chancolor = psd_now_chancolor + ([1 1 1] - psd_now_chancolor) * psd_held_colorfadefactor;
                        psd_held_channame = psd_now_channame;
                        psd_held_filttext = psd_now_filttext;
                        set(h_hintbar, 'String', 'Current PSD plot is held. Use arrow keys to move to a different time or channel to compare!');
                        end
                    end
                else
                    set(h_hintbar, 'String', 'PSD hold (the h button) only works in the PSD window.');
                end
            case 'p'
                %2024-11-06: Quickly toggle PSD plot
                if ~ishandle(viewhand_psd)
                    f_psd_plot(hObject, []);
                else
                    set(h_hintbar, 'String', '');
                    close(viewhand_psd);
                end
            case 'leftarrow'
                if Ctrl
                    f_xzoomout(hObject, []);
                elseif Alt && EventEnable
                    if Shift
                        u = find(EventTimes(:,1) < XLim(1), 1, 'last');
                    else
                        u = find(EventTimes(:,1) < (XLim(1) + XLim(2))/2, 1, 'last');
                    end
                    if ~isempty(u)
                        if u == centered_on_eventnum
                            % We already centered on this. Move one over
                            u = u - 1;
                            if u <= 0
                                u = 1;
                            end
                        end
                        centered_on_eventnum = u;
                        XRange = XLim(2) - XLim(1);
                        if Shift
                            % 20240917: If the Shift key is also held
                        XLim(1) = EventTimes(u,1);
                        XLim(2) = EventTimes(u,1) + XRange;
                        else
                            XLim(1) = EventTimes(u,1) - XRange/2;
                            XLim(2) = EventTimes(u,1) + XRange/2;
                        end
                        resnap_pan();
                        if EventTimes(u,1) == EventTimes(u,2)
                            et = sprintf(['%.' num2str(DecimalSecondsResolution) 'f s'], EventTimes(u,1));
                        else
                            et = sprintf(['%.' num2str(DecimalSecondsResolution) 'f -- %.' num2str(DecimalSecondsResolution) 'f s'], EventTimes(u,1), EventTimes(u,2));
                        end
                        if Shift
                            set(h_hintbar, 'String', ['Left-aligned on event #' num2str(u) '/' num2str(Nevents) ' [' et ']'  ': ' EventTimeStamps{u,2} '.']);
                        else
                            set(h_hintbar, 'String', ['Centered on event #' num2str(u) '/' num2str(Nevents) ' [' et ']'  ': ' EventTimeStamps{u,2} '.   Hint: Hold both Alt and Shift to align to the left instead of center.']);
                        end
                    elseif XLim(1) < min(EventTimes(:,1)) && XLim(1) > 0
                        set(h_hintbar, 'String', 'There are no more events to the left!');
                    end
                elseif Alt && ~EventEnable
                    if XLim(1) > 0
                        f_panleft(hObject, 'bigger');
                    end
                elseif Shift
                    f_panleft(hObject, 'smaller');
                else
                    f_panleft(hObject, []);
                end
            case 'rightarrow'
                if Ctrl
                    f_xzoomin(hObject, []);
                elseif Alt && EventEnable
                    if Shift
                        u = find(EventTimes(:,1) > XLim(1), 1, 'first');
                    else
                        u = find(EventTimes(:,1) > (XLim(1)+XLim(2))/2, 1, 'first');
                    end
                    if ~isempty(u)
                        if u == centered_on_eventnum
                            % We already centered on this. Move one over
                            u = u + 1;
                            if u > Nevents
                                u = Nevents;
                            end
                        end
                        centered_on_eventnum = u;
                        XRange = XLim(2) - XLim(1);
                        if Shift
                            % 20240917: If the Shift key is also held
                            XLim(1) = EventTimes(u,1);
                            XLim(2) = EventTimes(u,1) + XRange;
                        else
                            XLim(1) = EventTimes(u,1) - XRange/2;
                            XLim(2) = EventTimes(u,1) + XRange/2;
                        end
                        resnap_pan();
                        if EventTimes(u,1) == EventTimes(u,2)
                            et = sprintf(['%.' num2str(DecimalSecondsResolution) 'f s'], EventTimes(u,1));
                        else
                            et = sprintf(['%.' num2str(DecimalSecondsResolution) 'f -- %.' num2str(DecimalSecondsResolution) 'f s'], EventTimes(u,1), EventTimes(u,2));
                        end
                        if Shift
                            set(h_hintbar, 'String', ['Left-aligned on event #' num2str(u) '/' num2str(Nevents) ' [' et ']'  ': ' EventTimeStamps{u,2} '.']);
                        else
                            set(h_hintbar, 'String', ['Centered on event #' num2str(u) '/' num2str(Nevents) ' [' et ']'  ': ' EventTimeStamps{u,2} '.   Hint: Hold both Alt and Shift to align to the left instead of center.']);
                        end
                    elseif XLim(2) > max(EventTimes(:,1)) && XLim(2) < Time_max
                        set(h_hintbar, 'String', 'There are no more events to the right!');
                    end
                elseif Alt && ~EventEnable
                    if XLim(2) < Time_max
                        f_panright(hObject, 'bigger');
                    end
                elseif Shift
                    f_panright(hObject, 'smaller');
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
        %set(h_hintbar, 'String', 'Ctrl left/right: change time scale. Ctrl up/down: change sensitivity. Shift left/right: scroll slowly. Alt left/right: go to events.');
    end


    function choose_psd_channel(chan)
        found = false;
        if ~found
            for s = randperm(length(plothand))
                if strcmp(getappdata(plothand(s), 'channame'), chan)
                    selected_plothand = plothand(s);
                    found = true;
                    break
                end
            end
        end
        if ~found
            for s = randperm(length(plothand))
                if strcmpi(getappdata(plothand(s), 'channame'), chan)
                    selected_plothand = plothand(s);
                    found = true;
                    break
                end
            end
        end
        if ~found && ~isempty(chan2idx(ChanNames, chan, 1))
            if ~PlotHold
                selchan = get(h_chansel_list, 'Value');
                i = chan2idx(ChanNames,chan,1);
                if isscalar(i) && ~any(selchan==i)
                    selchan = union(selchan,i);
                    set(h_chansel_list, 'Value', selchan);
                    f_chansel_confirm([], []);
                    for s = randperm(length(plothand))
                        if strcmp(getappdata(plothand(s), 'channame'), chan)
                            selected_plothand = plothand(s);
                            found = true;
                            break
                        end
                    end
                    if ~found
                        for s = randperm(length(plothand))
                            if strcmpi(getappdata(plothand(s), 'channame'), chan)
                                selected_plothand = plothand(s);
                                found = true;
                                break
                            end
                        end
                    end
                end
            else
                set(h_hintbar, 'String', 'Cannot center on an unselected channel when plotting is paused. Unpause first.');
            end
        end


        if found
            loc = 0 - chansep*getappdata(selected_plothand, 'chanind');
            while loc > YLim(2)
                % Located above the current pan area
                YLim = YLim + chansep;
            end
            while loc < YLim(1)
                % Located below the current pan area
                YLim = YLim - chansep;
            end
            if ishandle(viewhand_psd)
                update_psd();
            else
                f_psd_plot([], []);
            end
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

    end


    function f_hold_switch(hObject, eventdata) %#ok<*INUSD>
        if ~isempty(hObject) && isscalar(hObject) && hObject <= -100000
            % This always disables plot hold
            PlotHold = 1;
        elseif ~isempty(hObject) && isscalar(hObject) && hObject <= -10000
            % This always enables plot hold
            PlotHold = 0;
        end

        if PlotHold
            set(h_hold_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            PlotHold = 0;
            set(h_hugetext, 'Visible', 'off');
            redraw();
            set(h_hintbar, 'String', '');
        else
            set(h_hold_switch, 'Value', 1, 'ForegroundColor', fontcolor_on2, 'FontWeight', fontweight_on2);
            if ~isempty(hObject) && hObject == -10000
                set(h_hugetext, 'String', 'Starting up...', 'Visible', 'on');
                drawnow
            elseif ~isempty(hObject) && hObject == -11000
                set(h_hugetext, 'String', 'Rescaling...', 'Visible', 'on');
                drawnow
            else
                set(h_hugetext, 'String', 'Plotting Paused', 'Visible', 'on');
                drawnow
            end
            PlotHold = 1;
            set(h_hintbar, 'String', 'Plotting paused');
        end
    end


    function f_verticaloverlapdisallow_switch(hObject, eventdata)
        if VerticalOverlapAllow <= 2 % means disallowed
            set(h_verticaloverlapdisallow_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            VerticalOverlapAllow = 100;
            redraw();
        else
            set(h_verticaloverlapdisallow_switch, 'Value', 1, 'ForegroundColor', fontcolor_on2, 'FontWeight', fontweight_on2);
            VerticalOverlapAllow = 1.2; % can go 10% over
            redraw();
        end
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
        
        envelope_update();
        render_update();
        enable_filter_switches();
        if EnvelopeFilter.state
            set(h_hintbar, 'String', 'Turned on envelope filter.');
        else
            set(h_hintbar, 'String', 'Turned off envelope filter.');
        end
    end


    function f_bpf_switch(hObject, eventdata)
        if isequal(get(h_bpf_switch, 'ForegroundColor'), fontcolor_on1)
            set(h_bpf_cutoff, 'String', '');
            f_bpf_cutoff([], []);
        else
            set(h_hintbar, 'String', get(h_bpf_switch, 'Tooltip'));
            set(h_bpf_switch, 'Value', 0);
        end
    end


    function f_bpf_cutoff(hObject, eventdata)
        disable_filter_switches();

        %Po240516 Automatically render filter if a valid range is entered.
        validate_bpf_cutoff();
        if is_nontrivial_bpf_cutoff()
            BandPassFilter.state = 1;
            set(h_bpf_switch, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1, 'Value', 1);
            bandpass_update();
            envelope_update();
            render_update();
        elseif BandPassFilter.state
            bandpass_update();
            envelope_update();
            render_update();
            BandPassFilter.state = 0;
            set(h_bpf_switch, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off, 'Value', 0);
        end

        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_ica_switch(hObject, eventdata)
        if isequal(get(h_ica_switch, 'ForegroundColor'), fontcolor_on1)
            if use_ltmat
                set(h_icasel_list, 'Value', 1);
            else
                set(h_icasel_list, 'Value', 1:length(get(h_icasel_list, 'String')));
            end
            f_icasel_confirm([], []);
        else
            set(h_hintbar, 'String', get(h_ica_switch, 'Tooltip'));
            set(h_ica_switch, 'Value', 0);
        end
    end

    function f_car_switch(hObject, eventdata)
        disable_filter_switches();
        if RerefFilter.state && isequal(RerefFilter.chanidx, setdiff(selchan, chan2idx(ChanNames, nofiltchannames, 1)))
            RerefFilter.state = 0;
            RerefFilter.chanidx = [];
        else
            RerefFilter.state = 1;
            RerefFilter.chanidx = selchan;
            % Po240610: Do not CAR/notch/filter/envelope matching channels:
            RerefFilter.chanidx = setdiff(RerefFilter.chanidx, chan2idx(ChanNames, nofiltchannames, 1));
        end
        reref_update();
        notch_update();
        bandpass_update();
        envelope_update();
        render_update();
        enable_filter_switches();
        if RerefFilter.state
            set(h_hintbar, 'String', 'Common average reference (CAR) enabled on the currently selected channels (denoted with ©). (Remember to re-CAR if changing channels)');
        else
            set(h_hintbar, 'String', 'Turned off CAR.');
        end
    end

    function f_evf_cutoff(hObject, eventdata)
        disable_filter_switches();
        if EnvelopeFilter.state
            envelope_update();
            render_update();
        end
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function disable_filter_switches()
        set([h_bpf_switch, h_ica_switch, h_car_switch, h_evf_switch, h_bpf_cutoff, h_evf_cutoff, h_notch_switch, h_notch_order, h_notch_qfactor], 'Enable', 'off');
        drawnow
    end

    function enable_filter_switches()
        set([h_bpf_switch, h_ica_switch, h_car_switch, h_evf_switch, h_bpf_cutoff, h_evf_cutoff, h_notch_switch, h_notch_order, h_notch_qfactor], 'Enable', 'on');
    end

    function disable_movement_switches()
        set([h_panleft h_panright h_panup h_pandown h_xzoomout h_xzoomin h_yzoomout h_yzoomin h_sepup h_sepdown], 'Enable' , 'off');
    end

    function enable_movement_switches()
        set([h_panleft h_panright h_panup h_pandown h_xzoomout h_xzoomin h_yzoomout h_yzoomin h_sepup h_sepdown], 'Enable' , 'on');
    end

    function f_notch_switch(hObject, eventdata)
        if isequal(get(h_notch_switch, 'ForegroundColor'), fontcolor_on1)
            set(h_notch_order, 'String', '0');
            f_notch_order([], []);
        else

            %2024-12-17 If PSD data tips have been used, maybe the user
            %want to notch these?
            %psd_datatips = [];
            psd_datatip_freqs = [];
            if ishandle(viewhand_psd_axe)
                try
                    psd_datatips = findobj(viewhand_psd_axe, 'Type', 'DataTip');
                    psd_datatip_freqs = unique(psd_now_fxx([psd_datatips(:).DataIndex]));
                catch
                    %psd_datatips = [];
                    psd_datatip_freqs = [];
                end
            end


            %2024-11-07 Allow the user to specify notch frequencies
            tmp_query1 = {
                'Enter a list of frequencies in Hz to notch (separated by space). If you enter only one frequency, the first three harmonics are also included automatically. If PSD window is open, frequencies indicated by the datatips have been added to this list.'
                ['Enter the Notch order. The order must be a positive even integer no more than ' num2str(notchordermax) '. Enter 0 to save the parameters without turning on the notch filter.']
                ['Enter the Q factor. The Q factor must be a positive number no more than ' num2str(notchqfactormax) '. Q factors are scaled automatically.']
                };
            tmp_defaultnotchfreqs = NotchFrequencies;
            if ~isempty(psd_datatip_freqs)
                tmp_defaultnotchfreqs = unique(round([NotchFrequencies(:).' psd_datatip_freqs(:).']));
            end
            
            tmp_default1 = {
                regexprep(num2str(tmp_defaultnotchfreqs), '\s+', ' ')
                num2str(NotchFilter.order)
                num2str(NotchFilter.qfactor)
                };
            tmp_answer = inputdlg(tmp_query1, 'Notch Filter Setup', [1 65], tmp_default1);
            if ~isempty(tmp_answer) && ~isempty(tmp_answer{1})
                b = regexp(tmp_answer{1}, '(\d*\.?\d*)', 'tokens');
                if ~isempty(b)
                    tmp_array = zeros(1,length(b));
                    for i = 1:length(b)
                        tmp_array(i) = str2double(b{i}{1});
                    end
                    tmp_array(~(isfinite(tmp_array) & tmp_array > 0 & mod(tmp_array,Fs/2) > 0)) = NaN;
                end
                tmp_array = unique(tmp_array(~isnan(tmp_array)));
                NotchFrequencies = tmp_array;
                tmp_additional = '';
                if any(NotchFrequencies >= Fs/2)
                    tmp_additional = '. Frequencies above Nyquist are adjusted for aliasing automatically.';
                end
                set(h_hintbar, 'String', ['Notch frequencies set to: ' num2str(NotchFrequencies) ' Hz' tmp_additional]);
            else
                set(h_hintbar, 'String', get(h_notch_switch, 'Tooltip'));
            end
            set(h_notch_switch, 'Value', 0);
            if ~isempty(tmp_answer) && ~isempty(tmp_answer{2}) && ~isempty(tmp_answer{3})
                set(h_notch_qfactor, 'String', tmp_answer{3});
                set(h_notch_order, 'String', tmp_answer{2});
                f_notch_order(hObject, eventdata);
            end
        end
    end

    function f_notch_order(hObject, eventdata)
        disable_filter_switches();
        %Po241106 Automatically render filter if parameters are valid.
        notch_update();
        bandpass_update();
        envelope_update();
        render_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    function f_notch_qfactor(hObject, eventdata)
        disable_filter_switches();
        %Po241106 Automatically render filter if parameters are valid.
        notch_update();
        bandpass_update();
        envelope_update();
        render_update();
        enable_filter_switches();
        set(h_hintbar, 'String', '');
    end

    
    function newchannames = channames_with_car_symbol()
        newchannames = ChanNames;
        inds = find(PerChannelFilterStates(1,:));
        for i = inds
            newchannames{i} = [newchannames{i} ' ©'];
        end
    end


    function reref_update()
        Signal_postreref = Signal_postica;
        PerChannelFilterStates(1:end,:) = false;

        if length(RerefFilter.chanidx) < 2
            %Po240611: Cannot CAR less than 2 channels
            RerefFilter.chanidx = [];
            RerefFilter.state = 0;
        end

        if RerefFilter.state
            FilterBusy = 1;
            set(h_bigtext, 'Visible', 'on', 'String', ['Re-referencing...']); drawnow;
            Signal_postreref(:,RerefFilter.chanidx) = car(Signal_postreref(:,RerefFilter.chanidx));
            set(h_bigtext, 'Visible', 'on', 'String', ['Finishing re-reference...']); drawnow;
            set(h_car_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            PerChannelFilterStates(1,RerefFilter.chanidx) = true;
            PerChannelFilterStates(2:end,RerefFilter.chanidx) = false; % Once this channel is re-referenced, all subsequent filters on this channel are invalidated
            set(h_chansel_list, 'String', channames_with_car_symbol());
            set(h_car_chancount, 'String', [num2str(length(RerefFilter.chanidx)) 'ch']);
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', '');
        else
            set(h_car_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            set(h_chansel_list, 'String', ChanNames);
            set(h_car_chancount, 'String', 'off');
        end
    end


    function notch_update()
        ord = str2double(get(h_notch_order, 'String'));
        qfactor = str2double(get(h_notch_qfactor, 'String'));
        
        if isfinite(ord) && mod(ord,2) == 1 && ord > 0
            ord = ord + 1;
        end

        if isfinite(qfactor) && qfactor > 0 && qfactor <= notchqfactormax
            NotchFilter.pendingqfactor = qfactor;
        end
        
        if isfinite(ord) && mod(ord,2) == 0 && ord >= 2 && ord <= notchordermax && isfinite(qfactor) && qfactor > 0 && qfactor <= notchqfactormax
            if ~isequal(NotchFilter.order,ord) || ~isequal(NotchFilter.qfactor,NotchFilter.pendingqfactor)
                PerChannelFilterStates(2:end,:) = false; %Filter parameters changed. Invalidating.
                NotchFilter.order = ord;
                NotchFilter.qfactor = NotchFilter.pendingqfactor;
            end
            NotchFilter.state = 1;
        else
            % Entering an invalid notch parameter will turn it off
            NotchFilter.order = 0;
            NotchFilter.state = 0;
            set(h_notch_order, 'String', num2str(NotchFilter.order));
        end

        if NotchFilter.state
            notch_failed = false;
            FilterBusy = 1;
            set(h_bigtext, 'Visible', 'on', 'String', ['Preparing multi-harmonic notch filters...']); drawnow;
            clear d Hd
            %2024-11-07: Notch filter to allow any list of frequencies and
            %can automatically calculate from aliasing
            if isscalar(NotchFrequencies)
                % If unspecified, filter the first 3 harmonics
                Freqs = NotchFrequencies(1)*(1:3);
            else
                Freqs = unique(NotchFrequencies(:)).';
            end
            % Get rid of frequencies exactly at multiples of Nyquist
            Freqs(mod(Freqs, Fs/2) == 0) = NaN;
            Freqs = Freqs(~isnan(Freqs));

            % Aliasing formula
            PerceivedFreqs = abs(Freqs - Fs .* round(Freqs ./ Fs));


            for h = 1:length(Freqs)
                % Create one for each frequency
                set(h_bigtext, 'Visible', 'on', 'String', ['Preparing ' num2str(Freqs(h)) ' Hz notch filter...']); drawnow;
                d = fdesign.notch('N,F0,Q',NotchFilter.order,PerceivedFreqs(h)/(Fs/2),NotchFilter.qfactor/Freqs(1)*Freqs(h));
                Hd{h} = design(d);
            end
            warning('off', 'signal:filtfilt:ParseSOS');
            warning('off', 'signal:filtfilt:ParseB');

            for ch = size(Signal_postreref,2):-1:1
                set(h_bigtext, 'Visible', 'on', 'String', ['Applying notch filters (' num2str(ch) ' chans to go)']); drawnow;
                if PerChannelFilterStates(2,ch) % Po240528: Only if this channel hasn't been done
                    continue
                end
                if ~ismember(ch, selchan) % Po240528: Only if this channel is plotted
                    continue
                end
                if ismember(ChanNames{ch}, nofiltchannames)
                    % Po240610: Do not CAR/notch/filter/envelope
                    continue
                end
                Signal_postnotch(:,ch) = Signal_postreref(:,ch);
                try
                    %2024-11-07 Notch in both directions and take the
                    %better half respectively

                    % Flip
                    tmp_flipped = flipud(Signal_postnotch(:,ch));

                    % Forward
                    for h = 1:length(Hd)
                        Signal_postnotch(:,ch) = filter(Hd{h},Signal_postnotch(:,ch));
                    end

                    % Reverse
                    for h = 1:length(Hd)
                        tmp_flipped = filter(Hd{h},tmp_flipped);
                    end
                    tmp_flipped = flipud(tmp_flipped);

                    % Combine
                    Signal_postnotch(:,ch) = Signal_postnotch(:,ch) .* splitnotch_factors + tmp_flipped .* (1-splitnotch_factors);

                    % Original code:
                    %for h = 1:length(Hd)
                    %    Signal_postnotch(:,ch) = filtfilt(Hd{h}.sosMatrix,Hd{h}.ScaleValues,Signal_postnotch(:,ch));
                    %end
                    PerChannelFilterStates(2,ch) = true;
                    PerChannelFilterStates(3:end,ch) = false; % Once this channel is notched, all subsequent filters on this channel are invalidated
                catch ex
                    warning(ex.identifier,'Notch filter crashed with this message: %s Notch filter will be turned off', ex.message);
                    notch_failed = true;
                    break
                end
            end

            if notch_failed
                Signal_postnotch = Signal_postreref;
                PerChannelFilterStates(2:end,:) = false;
                set(h_bigtext, 'Visible', 'on', 'String', ['WARNING: Notch filter crashed.']); drawnow;
                set(h_notch_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
                NotchFilter.order = 0;
                NotchFilter.state = 0;
                set(h_notch_order, 'String', num2str(NotchFilter.order));
                set(h_notch_qfactor, 'String', num2str(NotchFilter.pendingqfactor));
                set(h_hintbar, 'String', 'Notch filter crashed and is turned off.');
            else
                set(h_bigtext, 'Visible', 'on', 'String', ['Finishing notch filter...']); drawnow;
                set(h_notch_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
                set(h_notch_order, 'String', num2str(NotchFilter.order));
                set(h_notch_qfactor, 'String', num2str(NotchFilter.pendingqfactor));
            end
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', '');
        else
            if ~isequal(Signal_postnotch, Signal_postreref)
                Signal_postnotch = Signal_postreref;
                PerChannelFilterStates(2:end,:) = false;
            end
            set(h_notch_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            set(h_notch_order, 'String', num2str(NotchFilter.order));
            set(h_notch_qfactor, 'String', num2str(NotchFilter.pendingqfactor));
        end
    end


    function str = num2strcompact(num)
        str = regexprep(num2str(num),' {3,}', '  ');
    end


    function bpf = parse_bpf_cutoff_string()
        try
            b = regexp(get(h_bpf_cutoff,'String'), '(-?[0-9./*+-]+)([ ,]+)(-?[0-9./*+-]+)', 'tokens', 'once');
        catch
            b = {};
        end
        if ~isempty(b)
            bpf = [str2num(b{1}) str2num(b{3})]; %#ok<ST2NM>
        else
            bpf = [];
        end
    end


    function bandpass_update()
        % Input validation
        bpf = parse_bpf_cutoff_string();
        %bpf = str2num(get(h_bpf_cutoff, 'String'), 'Evaluation', 'restricted'); %#ok<ST2NM>
        if length(bpf) ~= 2 || numel(bpf) ~= 2
            bpf = [0 Fs/2];
        end
        if ~isfinite(bpf(1)) || bpf(1) < 0 || imag(bpf(1)) ~= 0
            bpf(1) = 0;
        end
        if ~isfinite(bpf(2)) || bpf(2) > Fs/2 || imag(bpf(2)) ~= 0
            bpf(2) = Fs/2;
        end

        if bpf(2) <= bpf(1)
            bpf(1) = 0;
            bpf(2) = Fs/2;
            BandPassFilter.state = 0;
        end

        if BandPassFilter.cutoff(1) ~= bpf(1) || BandPassFilter.cutoff(2) ~= bpf(2)
            PerChannelFilterStates(3:end,:) = false; %Filter parameters changed. Invalidating.
            BandPassFilter.cutoff = bpf(1:2);
            set(h_bpf_cutoff, 'String', num2strcompact(BandPassFilter.cutoff));
        end
        
        % If after validation still enabled, do filtering
        if BandPassFilter.state
            FilterBusy = 1;
            set(h_bigtext, 'Visible', 'on', 'String', ['Preparing Butterworth filter...']); drawnow;

            for ch = size(Signal_postnotch,2):-1:1
                if PerChannelFilterStates(3,ch) % Po240528: Only if this channel hasn't been done
                    continue
                end
                if ~ismember(ch, selchan) % Po240528: Only if this channel is plotted
                    continue
                end
                if ismember(ChanNames{ch}, nofiltchannames)
                    % Po240610: Do not CAR/notch/filter/envelope
                    continue
                end
                if bpf(1) > 0 && bpf(2) < Fs/2
                    set(h_bigtext, 'Visible', 'on', 'String', ['Applying band-pass filter (' num2str(ch) ' chans to go)']); drawnow;
                    [Signal_postbutter(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [bpf FilterOrder], 'pass', 'butter');
                elseif bpf(1) > 0
                    set(h_bigtext, 'Visible', 'on', 'String', ['Applying high-pass filter (' num2str(ch) ' chans to go)']); drawnow;
                    [Signal_postbutter(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [bpf(1) FilterOrder], 'high', 'butter');
                elseif bpf(2) < Fs/2
                    set(h_bigtext, 'Visible', 'on', 'String', ['Applying low-pass filter (' num2str(ch) ' chans to go)']); drawnow;
                    [Signal_postbutter(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [bpf(2) FilterOrder], 'low', 'butter');
                else
                    Signal_postbutter = Signal_postnotch;
                    break;
                end
                while any(FilterInfo.ButterUnstable)
                    if FilterOrder <= 1
                        Signal_postbutter(:,ch) = Signal_postnotch(:,ch);
                        BandPassFilter.state = 0;
                        break
                    end
                    FilterOrder = ceil(FilterOrder / 2);
                    if bpf(1) > 0 && bpf(2) < Fs/2
                        [Signal_postbutter(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [bpf FilterOrder], 'pass', 'butter');
                    elseif bpf(1) > 0
                        [Signal_postbutter(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [bpf(1) FilterOrder], 'high', 'butter');
                    elseif bpf(2) < Fs/2
                        [Signal_postbutter(:,ch), FilterInfo] = freqfilter(Signal_postnotch(:,ch), Fs, [bpf(2) FilterOrder], 'low', 'butter');
                    end
                end
                if ~BandPassFilter.state
                    %Signal_postbutter = Signal_postnotch;
                    break
                end
                PerChannelFilterStates(3,ch) = true;
                PerChannelFilterStates(4:end,ch) = false; % Once this channel is buttered, all subsequent filters on this channel are invalidated
            end

        end
        if ~BandPassFilter.state
            if ~isequal(Signal_postbutter, Signal_postnotch)
                Signal_postbutter = Signal_postnotch;
                PerChannelFilterStates(3:end,:) = false;
            end
            set(h_bpf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        FilterBusy = 0;
    end


%Po240524: Butterworth filter and Envelope filter are now separate.
    function envelope_update()
        % Input validation
        evf = str2double(get(h_evf_cutoff, 'String'));

        if ~isfinite(evf)
            evf = Fs;
        end

        if evf >= Fs/2 || evf <= 0
            EnvelopeFilter.state = 0;
            set(h_evf_cutoff, 'String', num2str(floor(Fs/2*99)/100));
            set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end

        if EnvelopeFilter.cutoff ~= evf
            PerChannelFilterStates(4:end,:) = false; %Filter parameters changed. Invalidating.
            EnvelopeFilter.cutoff = evf;
        end

        % If after validation still enabled, do filtering
        if EnvelopeFilter.state
            FilterBusy = 1;
            set(h_bigtext, 'Visible', 'on', 'String', ['Preparing envelope filter...']); drawnow;
            for ch = size(Signal_postbutter,2):-1:1
                if PerChannelFilterStates(4,ch) % Po240528: Only if this channel hasn't been done
                    continue
                end
                if ~ismember(ch, selchan) % Po240528: Only if this channel is plotted
                    continue
                end
                if ismember(ChanNames{ch}, nofiltchannames)
                    % Po240610: Do not CAR/notch/filter/envelope
                    continue
                end
                set(h_bigtext, 'Visible', 'on', 'String', ['Applying envelope filter (' num2str(ch) ' chans to go)']); drawnow;
                [Signal_postenvelope(:,ch), FilterInfo] = freqfilter(Signal_postbutter(:,ch).^2, Fs, [evf FilterOrder], 'low', 'butter');
                while any(FilterInfo.ButterUnstable)
                    if FilterOrder <= 1
                        Signal_postenvelope(:,ch) = Signal_postbutter(:,ch);
                        EnvelopeFilter.state = 0;
                        set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
                        break
                    end
                    FilterOrder = ceil(FilterOrder / 2);
                    [Signal_postenvelope(:,ch), FilterInfo] = freqfilter(Signal_postbutter(:,ch).^2, Fs, [evf FilterOrder], 'low', 'butter');
                end
                if ~EnvelopeFilter.state
                    %Signal_postbutter = Signal_postbutter;
                    break
                end
                PerChannelFilterStates(4,ch) = true;
            end
            if EnvelopeFilter.state
                SigBandwidth = min(SigBandwidth,evf);
            end

        end
        if ~EnvelopeFilter.state
            if ~isequal(Signal_postenvelope, Signal_postbutter)
                Signal_postenvelope = Signal_postbutter;
                PerChannelFilterStates(4:end,:) = false;
            end
            set(h_evf_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
        end
        FilterBusy = 0;
    end


    function render_update() 
        set(h_bigtext, 'Visible', 'on', 'String', 'Rendering...'); drawnow;
        redraw();
        set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
    end        


    function f_zscore_switch(hObject, eventdata)
        set(h_zscore_switch, 'Enable', 'off');
        %set(h_zscore_state, 'String', 'Wait'); 
        drawnow;
        plot_was_held = 0;
        if ~PlotHold
            f_hold_switch(-11000, []);
            plot_was_held = 1;
        end
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
        %redraw();
        if plot_was_held
            if isempty(hObject) || numel(hObject) ~= 1 || hObject ~= -10000
                % Turn off plot hold (except during startup sequence)
                f_hold_switch(-100000, []);
            end
        else
            redraw();
        end
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


    function L = validate_bpf_cutoff() %#ok<STOUT>
        bpf = parse_bpf_cutoff_string();
        %bpf = str2num(get(h_bpf_cutoff, 'String'), 'Evaluation', 'restricted'); %#ok<ST2NM>
        changesmade = 0;
        if length(bpf) ~= 2 || numel(bpf) ~= 2
            bpf = [0 Fs/2];
            changesmade = 1;
        end
        if ~isfinite(bpf(1)) || bpf(1) < 0 || imag(bpf(1)) ~= 0
            bpf(1) = 0;
            changesmade = 1;
        end
        if ~isfinite(bpf(2)) || bpf(2) > Fs/2 || imag(bpf(2)) ~= 0
            bpf(2) = Fs/2;
            changesmade = 1;
        end
        if changesmade
        set(h_bpf_cutoff, 'String', num2strcompact(bpf));
        end
    end

    function L = is_nontrivial_bpf_cutoff()
        L = false;
        bpf = parse_bpf_cutoff_string();
        %bpf = str2num(get(h_bpf_cutoff, 'String'), 'Evaluation', 'restricted'); %#ok<ST2NM>
        if (bpf(1) > 0 || bpf(2) < Fs/2) && bpf(2) > bpf(1)
            L = true;
        end
        return
    end


    function f_xzoomin(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        set(h_hintbar, 'String', 'Zooming in. Please wait.'); drawnow;

        acp = get(axehand, 'CurrentPoint');
        XCursor = acp(1); % this is the time point in seconds where the cursor is pointing at when scrolled
        XRange = XLim(2)-XLim(1);
        % Ratio of cursor position:
        XCursorRatio = (XCursor - XLim(1)) / XRange;
        %if XCursorRatio < 0
        %    XCursorRatio = 0;
        %elseif XCursorRatio > 1
        %    XCursorRatio = 1;
        %end
        
        [~,in] = min(abs(XRange - PermittedXZoomRanges));
        if in > 1
            XRange = PermittedXZoomRanges(in-1);
        else
            XRange = PermittedXZoomRanges(in);
        end
        XLim(1) = XCursor - XCursorRatio*XRange;
        XLim(2) = XCursor + (1-XCursorRatio)*XRange;
        %XLim(2) = XLim(1) + XRange;
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
        set(h_hintbar, 'String', 'Zooming out. Please wait.'); drawnow;

        acp = get(axehand, 'CurrentPoint');
        XCursor = acp(1); % this is the time point in seconds where the cursor is pointing at when scrolled
        XRange = XLim(2)-XLim(1);
        % Ratio of cursor position:
        XCursorRatio = (XCursor - XLim(1)) / XRange;
        %if XCursorRatio < 0
        %    XCursorRatio = 0;
        %elseif XCursorRatio > 1
        %    XCursorRatio = 1;
        %end

        [~,in] = min(abs(XRange - PermittedXZoomRanges));
        if in < length(PermittedXZoomRanges)
            XRange = PermittedXZoomRanges(in+1);
        else
            XRange = PermittedXZoomRanges(in);
        end
        
        if XRange > Time_max
            XRange = Time_max;
        end
        XLim(1) = XCursor - XCursorRatio*XRange;
        XLim(2) = XCursor + (1-XCursorRatio)*XRange;
        if XLim(1) < 0 && XLim(2) > Time_max
            XLim = [0 Time_max];
        elseif XLim(1) < 0
            XLim = XLim - XLim(1);
        elseif XLim(2) > Time_max
            XLim(2) = Time_max;
        end
        %XLim(2) = XLim(1) + XRange;
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
        set(h_hintbar, 'String', 'Zooming in. Please wait.'); drawnow;
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
        set(h_hintbar, 'String', 'Zooming out. Please wait.'); drawnow;
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
        set(h_hintbar, 'String', 'Panning. Please wait.'); drawnow;
        panfrac = panfrac_noncursor_mode;
        panmod = 1;
        if ~isempty(eventdata)
            switch eventdata
                case 'bigger'
                    panmod = 10;
                case 'smaller'
                    panmod = 0.1;
            end
        elseif CursorEnable
            panfrac = panfrac_cursor_mode;
        end
        XRange = XLim(2)-XLim(1);
        if hpanamount_mode == 1
            XLim = XLim - panfrac*panmod;
        else
            XLim = XLim - XRange*panfrac*panmod;
        end
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
        set(h_hintbar, 'String', 'Panning. Please wait.'); drawnow;
        panfrac = panfrac_noncursor_mode;
        panmod = 1;
        if ~isempty(eventdata)
            switch eventdata
                case 'bigger'
                    panmod = 10;
                case 'smaller'
                    panmod = 0.1;
            end
        elseif CursorEnable
            panfrac = panfrac_cursor_mode;
        end
        XRange = XLim(2)-XLim(1);
        if hpanamount_mode == 1
            XLim = XLim + panfrac*panmod;
        else
            XLim = XLim + XRange*panfrac*panmod;
        end
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end

    function select_next_plothand()
        try
            f = find(selchan == getappdata(selected_plothand, 'chanind'));
            if ~isempty(f)
                s = mod(f+1 -1,length(selchan))+1;
                for i = randperm(length(plothand))
                    if getappdata(plothand(i),'chanind') == selchan(s)
                        selected_plothand = plothand(i);
                        break
                    end
                end
                while 0 - chansep*getappdata(selected_plothand, 'chanind') < YLim(1) && YLim(1) > -chansep*size(Signal,2)*3
                    YLim = YLim - chansep;
                end
                while 0 - chansep*getappdata(selected_plothand, 'chanind') > YLim(2) && YLim(2) < 0+chansep*size(Signal,2)*3
                    YLim = YLim + chansep;
                end
            end
        end
    end
    function select_prev_plothand()
        try
            f = find(selchan == getappdata(selected_plothand, 'chanind'));
            if ~isempty(f)
                s = mod(f-1 -1,length(selchan))+1;
                for i = randperm(length(plothand))
                    if getappdata(plothand(i),'chanind') == selchan(s)
                        selected_plothand = plothand(i);
                        break
                    end
                end
                while 0 - chansep*getappdata(selected_plothand, 'chanind') > YLim(2) && YLim(2) < 0+chansep*size(Signal,2)*3
                    YLim = YLim + chansep;
                end
                while 0 - chansep*getappdata(selected_plothand, 'chanind') < YLim(1) && YLim(1) > -chansep*size(Signal,2)*3
                    YLim = YLim - chansep;
                end
            end
        end
    end


    function f_panup(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        set(h_hintbar, 'String', 'Panning. Please wait.'); drawnow;

        %2024-10-08 Pan up/down behaves differently in the main window vs. in the PSD window
        if isequal(hObject, viewhand_psd)
            select_prev_plothand();
            update_psd();
        else
            YRange = YLim(2)-YLim(1);
            YLim = YLim + YRange;
            set(axehand, 'YLim', YLim);
        end
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
        set(h_hintbar, 'String', 'Panning. Please wait.'); drawnow;

        %2024-10-08 Pan up/down behaves differently in the main window vs. in the PSD window
        if isequal(hObject, viewhand_psd)
            select_next_plothand();
            update_psd();
        else
            YRange = YLim(2)-YLim(1);
            YLim = YLim - YRange;
            set(axehand, 'YLim', YLim);
        end
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


    function f_hintbar(hObject, eventdata)
        persistent t_lastclicked
        if isempty(t_lastclicked)
            t_lastclicked = 0;
        end
        if now - t_lastclicked < 0.25 / 86400 %#ok<*TNOW1>
            clipboard('copy', strrep(get(h_hintbar,'String'), infolabels_text_highdeviations, ''));
            t_lastclicked = 0;
        else
            t_lastclicked = now;
        end
        
    end


    function f_xspan_edit1(hObject, eventdata)
        if FilterBusy || MovementBusy
            set(h_xspan_edit1, 'String', num2str(XLim(1)));
            set(h_xspan_samplenumber1, 'String', num2str(max(1,XLim(1)*SampleRate+1)));
            return;
        end
        set(h_hintbar, 'String', 'Changing horizontal range. Please wait.'); drawnow;
        v = str2double(get(h_xspan_edit1, 'String'));
        if isfinite(v) && imag(v) == 0
            if  v < XLim(2)
                XLim(1) = v;
            else
                XLim(2) = v + (XLim(2)-XLim(1));
                XLim(1) = v;
            end
        end
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end



    function f_xspan_edit2(hObject, eventdata)
        if FilterBusy || MovementBusy
            set(h_xspan_edit2, 'String', num2str(XLim(2)));
            set(h_xspan_samplenumber2, 'String', num2str(min(Ntp,XLim(2)*SampleRate+1)));
            return;
        end
        set(h_hintbar, 'String', 'Changing horizontal range. Please wait.'); drawnow;
        v = str2double(get(h_xspan_edit2, 'String'));
        if isfinite(v) && imag(v) == 0
            if v > XLim(1)
                XLim(2) = v;
            else
                XLim(1) = v - (XLim(2)-XLim(1));
                XLim(2) = v;
            end
        end
        set(axehand, 'XLim', XLim);
        MovementBusy = 1;
        resnap_pan();
        MovementBusy = 0;
        set(h_hintbar, 'String', ['XLim changed to ' num2str(XLim(1)) ' -- ' num2str(XLim(2)) ' seconds']);
    end


    function f_hpanamount_edit(hObject, eventdata)
        try
            uin = get(h_hpanamount_edit,'String');
            uin = regexprep(uin, '\s+$', '');
            uin = regexprep(uin, 'weeks', 'w');
            uin = regexprep(uin, 'week', 'w');
            uin = regexprep(uin, 'wks', 'w');
            uin = regexprep(uin, 'wk', 'w');
            uin = regexprep(uin, 'days', 'd');
            uin = regexprep(uin, 'day', 'd');
            uin = regexprep(uin, 'hours', 'h');
            uin = regexprep(uin, 'hour', 'h');
            uin = regexprep(uin, 'hrs', 'h');
            uin = regexprep(uin, 'hr', 'h');
            uin = regexprep(uin, 'minutes', 'm');
            uin = regexprep(uin, 'minute', 'm');
            uin = regexprep(uin, 'mins', 'm');
            uin = regexprep(uin, 'min', 'm');
            uin = regexprep(uin, 'seconds', 's');
            uin = regexprep(uin, 'second', 's');
            uin = regexprep(uin, 'secs', 's');
            uin = regexprep(uin, 'sec', 's');
            uin = regexprep(uin, 'milliseconds', 'ms');
            uin = regexprep(uin, 'millisecond', 'ms');
            uin = regexprep(uin, 'millisec', 'ms');
            uin = regexprep(uin, 'microseconds', 'us');
            uin = regexprep(uin, 'microsecond', 'us');
            uin = regexprep(uin, 'microsec', 'us');
            uin = regexprep(uin, 'nanoseconds', 'ns');
            uin = regexprep(uin, 'nanosecond', 'ns');
            uin = regexprep(uin, 'nanosec', 'ns');
            uin = regexprep(uin, 'picoseconds', 'ps');
            uin = regexprep(uin, 'picosecond', 'ps');
            uin = regexprep(uin, 'picosec', 'ps');
            uin = cell_to_string(regexp(uin, '[0123456789.smunphdwe-% ]', 'match'), '');
            % Get the unit part
            uinb = regexp(uin, '^([0123456789e-]+)\s?([pnmshdw%]+)$', 'tokens', 'once');
            if ~isempty(uinb)
                uinunit = uinb{2};
                uinval = uinb{1};
            else
                uinunit = 's';
                uinval = uin;
            end

            % Try parse the value part
            uinval = str2double(uinval);
            if isnan(uinval) || uinval < 0
                error('cannot parse');
            end
            if isequal(uinunit, '%')
                panfrac_noncursor_mode = uinval/100;
                hpanamount_mode = 0;
                set(h_hpanamount_edit, 'String', sprintf('%g%s',uinval,uinunit));
                set(h_hintbar, 'String', sprintf('Left and Right arrow keys will now pan the chart by %g%% of screen. (You can still hold Shift to pan 1/10th of the normal amount.)', uinval));
            else
                switch uinunit
                    case 'ps'
                        uinval = uinval * 1e-12;
                    case 'ns'
                        uinval = uinval * 1e-9;
                    case 'us'
                        uinval = uinval * 1e-6;
                    case 'ms'
                        uinval = uinval * 1e-3;
                    case 'm'
                        uinval = uinval * 60;
                    case 'h'
                        uinval = uinval * 3600;
                    case 'd'
                        uinval = uinval * 86400;
                    case 'w'
                        uinval = uinval * 86400*7;
                end
                panfrac_noncursor_mode = uinval;
                hpanamount_mode = 1;
                set(h_hpanamount_edit, 'String', sprintf('%gs',uinval));
                set(h_hintbar, 'String', sprintf('Left and Right arrow keys will now pan the chart by %g seconds. (You can still hold Shift to pan 1/10th of the normal amount.)', uinval));
            end
        catch
            uin = '100%';
            panfrac_noncursor_mode = 1.00;
            hpanamount_mode = 0;
            set(h_hpanamount_edit,'String',uin);
            set(h_hintbar, 'String', sprintf('Left and Right arrow keys will now pan the chart by %g%% of screen. (You can still hold Shift to pan 1/10th of the normal amount.)', panfrac_noncursor_mode*100));
        end
    end


    function f_sepup(hObject, eventdata)
        if FilterBusy || MovementBusy
            return;
        end
        set(h_hintbar, 'String', 'Decreasing vertical sensitivity. Please wait.'); drawnow;
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
        set(h_hintbar, 'String', 'Increasing vertical sensitivity. Please wait.'); drawnow;
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
        XRange = round(XRange*1e9)/1e9; %Po240712: Fixed rounding error
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
        set(h_xspan_samplenumber1, 'String', num2str(max(1,XLim(1)*SampleRate+1)));
        set(h_xspan_edit2, 'String', num2str(XLim(2)));
        set(h_xspan_samplenumber2, 'String', num2str(min(Ntp,XLim(2)*SampleRate+1)));
        redraw();
        
    end


    function XLim = round_xlim(XLim, XRange)
        xl1 = round(XLim(1)*1e9)/1e9; %Po240712: Fixed rounding error
        if XRange >= 1
            XLim(1) = floor(xl1*FineSnapScale)/FineSnapScale;
        elseif XRange >= 0.1
            XLim(1) = floor(xl1*10*FineSnapScale)/10/FineSnapScale;
        elseif XRange >= 0.01
            XLim(1) = floor(xl1*100*FineSnapScale)/100/FineSnapScale;
        elseif XRange >= 0.001
            XLim(1) = floor(xl1*1000*FineSnapScale)/1000/FineSnapScale;
        elseif XRange >= 0.0001
            XLim(1) = floor(xl1*10000*FineSnapScale)/10000/FineSnapScale;
        elseif XRange >= 0.00001
            XLim(1) = floor(xl1*100000*FineSnapScale)/100000/FineSnapScale;
        elseif XRange >= 0.000001
            XLim(1) = floor(xl1*1000000*FineSnapScale)/1000000/FineSnapScale;
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
        set(h_xspan_samplenumber1, 'String', num2str(max(1,XLim(1)*SampleRate+1)));
        set(h_xspan_edit2, 'String', num2str(XLim(2)));
        set(h_xspan_samplenumber2, 'String', num2str(min(Ntp,XLim(2)*SampleRate+1)));
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
        sd = zeros(1,Nsch);
        tdrawupdate = tic;
        for ch = Nsch:-1:1
            
            if -chansep/2-chansep*ch < YLim(1) || chansep/2-chansep*ch > YLim(2)
                % out of plotting range
                %fprintf('%i is out of range\n', ch);
                set(plothand(ch), 'Visible', 'off');
                continue;
            end
            
            if ~PlotHold

                if InfoLabelEnable %Turn off the I-box because it lags
                    f_infolabel_enable([],[]);
                end

                if ~isempty(EventTimePoints)
                    % NaN around blankaround_stitch_samples
                    l1l = EventTimePoints(:,1) >= t1 & EventTimePoints(:,1) <= t2;
                    l2l = EventTimePoints(:,2) >= t1 & EventTimePoints(:,2) <= t2;
                    lcomlist = find(l1l | l2l);
                else
                    lcomlist = [];
                end

                if ~isempty(lcomlist)
                    tmp = Signal_postenvelope(:,selchan(ch));
                    for i2 = 1:length(lcomlist)
                        i = lcomlist(i2);
                        if blankaround_stitch_samples > 0 && blankaround_stitch_samples <= EventTimePoints(i,2)-EventTimePoints(i,1)+1
                            tmp(round(EventTimePoints(i,1) + [0:blankaround_stitch_samples-1]),:) = NaN;
                            tmp(round(EventTimePoints(i,2) - [blankaround_stitch_samples-1:-1:0]),:) = NaN;
                        end
                    end
                    Signal_postdecimation = tmp(t1:BLIM:t2,:);
                else
                    Signal_postdecimation = Signal_postenvelope(t1:BLIM:t2,selchan(ch));
                end

                %Po250324: NaN the regions where there used to be nan
                Signal_postdecimation(Signal_nonfinite(t1:BLIM:t2,selchan(ch)),:) = NaN;
                
                Time4 = Time(t1:BLIM:t2);
                if ScreenLimitedDownsampling && length(Time4) > 2*SLD_H
                    tdiff = Time(t2) - Time(t1);
                    Fs_pref = 2^nextpow2(SLD_H / tdiff);
                    if Fs_pref < Fs/BLIM
                        [Signal_postdecimation, Time4] = downsamplecustom(Signal_postdecimation, Time4, SLD_H);
                    end
                end
                
                YDATA = Signal_postdecimation;
                sd(ch) = std(YDATA);

                
                
                if ZscoreFilter.state
                    YDATA = nanzscore(YDATA)*ZscoreFilter.multiplier;
                else
                    YDATA = YDATA - nanmedian(YDATA); %#ok<*NANMEDIAN>
                end

                %Po241217: Cap the traces above/below chansep
                YDATA(YDATA>chansep/2*VerticalOverlapAllow) = chansep/2*VerticalOverlapAllow;
                YDATA(YDATA<-chansep/2*VerticalOverlapAllow) = -chansep/2*VerticalOverlapAllow;

                set(plothand(ch), 'XData', Time4, 'YData', YDATA - chansep*ch, 'Visible', 'on');
                if size(YDATA,1) <= MarkerZoomThreshold
                    set(plothand(ch), 'Marker', ZoomedInMarker);
                else
                    set(plothand(ch), 'Marker', ZoomedOutMarker);
                end
                set(plothand(ch), 'Color', Kolor(mod(selchan(ch)-1,Nkolor)+1,:));
                setappdata(plothand(ch), 'chanind', selchan(ch));
                setappdata(plothand(ch), 'channame', ChanNames{selchan(ch)});
                if ismember(ChanNames{selchan(ch)}, nofiltchannames)
                    set(plothand(ch), 'Color', [0 0 0]);
                end
                update_psd();
                update_infolabels();

                set(h_chansel_title, 'String', sprintf('%i chans plotted', length(selchan)));

            else
                set(plothand(ch), 'Visible', 'off');
            end

            
            if toc(tdrawupdate) > 0.2
                set(axehand, 'YColor', dynamicBusyYColor(ch,Nsch));
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
            %2023-05-20: Try to set labels to avoid the Yticks
            %yt = get(axehand, 'YTick');
            YPos = YLim(1)+(YLim(2)-YLim(1))/100*[1 2 3];
            YPos = YPos(YPos>YLim(1) & YPos<YLim(2));
            %YPos = stagger_odds_evens(YPos);


            linevisible = get(eventplothand, 'Visible');
            if any([linevisible{:}] == 0)
                set(eventplothand, 'Visible', 'on');
            end

            if isempty(YPos)
                YPos = sort(mean(YLim)+diff(YLim)/20*([-8:2:8]), 'descend');
            end
            NYPos = length(YPos);
            for i = size(EventTimeStamps,1):-1:1
                if EventTimes(i,1) < XLim(1) || EventTimes(i,1) > XLim(2)
                    continue
                end
                % Po240516: Repositioned event labels to top of chart
                tmp = get(eventtexthand(i), 'Position');
                tmp(2) = YLim(2); %YPos(mod(i-1,NYPos)+1);
                %larrow = '«';
                larrow = '';
                rarrow = '';
                horali = 'left';
                verali = 'top';
                %if EventTimeStamps{i,1} - XLim(1) > (XLim(2) - XLim(1))*0.5
                %    larrow = '';
                %    rarrow = '»';
                %    horali = 'right';
                %    tmp(2) = 0; %tmp(2)+(YLim(2)-YLim(1))/100*96;
                %end
                set(eventtexthand(i), 'Position', tmp, 'String', [larrow EventTimeStamps{i,2} rarrow], 'HorizontalAlignment', horali, 'VerticalAlignment', verali, 'Margin', 0.01);
                %set(eventtexthand(i), 'Position', tmp, 'String', [larrow EventTimeStamps{i,2} rarrow], 'HorizontalAlignment', horali, 'Rotation', 90, 'LineStyle', ':', 'EdgeColor', EventKolor);
            end

            % Po240516: A second round to un-overlap labels at or too close
            % to the same timestamp
            if iscell(get(eventtexthand, 'Position'))
                eventtexthandpositions = cell2mat(get(eventtexthand, 'Position'));
            else
                eventtexthandpositions = get(eventtexthand, 'Position');
            end
            eventtexthandpositions(:,3) = 1:size(eventtexthandpositions,1);
            %eventtexthandpositions = eventtexthandpositions(XLim(1) <= cell2mat(EventTimeStamps(:,1)) & XLim(2) >= cell2mat(EventTimeStamps(:,1)),:);
            eventtexthandpositions = eventtexthandpositions(XLim(1) <= EventTimes(:,1) & XLim(2) >= EventTimes(:,1),:);

            if size(eventtexthandpositions,1) > 1
                tmp2 = eventtexthandpositions(:,1:2);
                tmp2(:,1) = tmp2(:,1) ./ (XLim(2) - XLim(1));
                tmp2(:,2) = tmp2(:,2) ./ (YLim(2) - YLim(1));
                tmp = squareform(pdist(tmp2));
                tmp_min_dist = EventTextMinDistApart;
                tmp_shifty_by = 0.01;
                %tmp_shiftx_by = 0.01;
                tmp_y_change = (YLim(2)-YLim(1))*tmp_shifty_by;
                %tmp_x_change = (XLim(2)-XLim(1))*tmp_shiftx_by;
                for i = 1:size(eventtexthandpositions,1)
                    for j = i+1:size(eventtexthandpositions,1)
                        while tmp(i,j) < tmp_min_dist && eventtexthandpositions(j,2) > YLim(1) + tmp_y_change
                            %fprintf('%i is too close to %i (dist = %g)\n', j, i, tmp(i,j));
                            eventtexthandpositions(j,2) = eventtexthandpositions(j,2) - tmp_y_change;
                            %eventtexthandpositions(j,1) = eventtexthandpositions(j,1) + tmp_x_change;
                            tmp2 = eventtexthandpositions(:,1:2);
                            tmp2(:,1) = tmp2(:,1) ./ (XLim(2) - XLim(1));
                            tmp2(:,2) = tmp2(:,2) ./ (YLim(2) - YLim(1));
                            tmp = squareform(pdist(tmp2));
                        end
                        tmp_nattempts = 0;
                        while tmp(i,j) < tmp_min_dist && tmp_nattempts < 1000
                            eventtexthandpositions(j,2) = rand(1) * (YLim(2) - YLim(1)) + YLim(1);
                            tmp2 = eventtexthandpositions(:,1:2);
                            tmp2(:,1) = tmp2(:,1) ./ (XLim(2) - XLim(1));
                            tmp2(:,2) = tmp2(:,2) ./ (YLim(2) - YLim(1));
                            tmp = squareform(pdist(tmp2));
                            tmp_nattempts = tmp_nattempts + 1;
                        end
                    end
                end
                for i = 1:size(eventtexthandpositions,1)
                    set(eventtexthand(eventtexthandpositions(i,3)), 'Position', eventtexthandpositions(i,1:2));
                end
            end
            clear tmp tmp_*
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
            l = XLim(1) <= EventTimes(:,1) & EventTimes(:,1) <= XLim(2);
            set(eventtexthand(l), 'Visible', 'on');
            set(eventtexthand(~l), 'Visible', 'off');
        end
        
        update_cursorline();
        check_axes_fontangle_and_resize();
        update_infolabels();
        
    end


    function c = dynamicBusyYColor(s, smax)
        c = BusyYColor*s/smax;
    end

    function check_axes_fontangle_and_resize()
        % Newer versions of MATLAB automatically rotate the axes text
        % if the font is too large, so when we detect rotation, we
        % should reduce font size to fit
        ang = xtickangle(axehand);
        while ang ~= 0
            f_axesfont_dec([], []);
            drawnow;
            ang = xtickangle(axehand);
        end
    end


    % function f_chansel_auto(hObject, eventdata)
    %     try
    %         cids = GetChannelsWithSimilarStdev(Signal_postbutter);
    %         set(h_chansel_list, 'Value', cids);
    %         f_chansel_confirm([], []);
    %     catch
    %     end
    % end


    function f_infolabel_enable(hObject, eventdata)
        if InfoLabelEnable
            InfoLabelEnable = 0;
            set(h_infolabel_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
            update_infolabels();
            set(h_hintbar, 'String', '');
        else
            InfoLabelEnable = 1;
            update_infolabels();
            set(h_infolabel_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            tmp1 = 'Local deviations in the currently displayed time range are shown for each channel. A pure sine wave of amplitude 1 µV (2 µVpp) is 1.0 deviation';
            if ~isempty(infolabels_deviations)
                tab = [infolabels_deviations; selchan]';
                [~,ia] = setdiff(tab(:,2),chan2idx(ChanNames,nofiltchannames,1));
                tab = tab(ia,:);
                tab = flipud(sortrows(tab));
                thres = max( median(infolabels_deviations) + 2*1.4826*mad(infolabels_deviations,1), 1.1 * median(infolabels_deviations));
                tmp2 = ChanNames(tab(tab(:,1) > thres,2)');
                if ~isempty(tmp2)
                    tmp1 = [infolabels_text_highdeviations cell_to_string(tmp2, '|')];
                end
            end
            set(h_hintbar, 'String', tmp1);
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

            if isempty(SavedPointsTable) && evalin('base', ['exist(''signalviewer_SavedPointsTable_' signalHashStr ''',''var'')'])
                q = questdlg('There is already a list of saved minimum/maximum points for the exact same data. Load this list?', 'Load from base workspace?', 'Yes', 'No', 'Yes');
                if isequal(q, 'Yes')
                    SavedPointsTable = evalin('base', ['signalviewer_SavedPointsTable_' signalHashStr]);
                end
            end
            CursorEnable = 1;
            %set(h_cursor_state, 'String', text_on);
            set(h_cursor_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
            set(h_hintbar, 'String', 'Click on signal to get time/value. Hold Shift to get local max. Hold Ctrl to get local min. Single click on the min/max text label to save the point into a list.');
        end
    end

    % function f_linemarker_switch(hObject, eventdata)
    %     if LinemarkerEnable
    %         LinemarkerEnable = 0;
    %         set(h_linemarker_switch, 'Value', 0, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off);
    %         update_linemarker();
    %     else
    %         LinemarkerEnable = 1;
    %         set(h_linemarker_switch, 'Value', 1, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1);
    %         update_linemarker();
    %     end
    % end

    function f_chansel_list(hObject, eventdata)
        set(h_hintbar, 'String', ['Hold Ctrl to select/deselect multiple. Hold Shift to select contiguous range. Click "' get(h_chansel_confirm, 'String') '" to update. For advanced selection, use "' get(h_chansel_commandentry, 'String') '"']);
        if ~isequal(selchan, get(h_chansel_list, 'Value'))
            set(h_chansel_confirm, 'BackgroundColor', backgroundcolor_buttonalert);
            set(h_chansel_reset, 'BackgroundColor', backgroundcolor_buttonalert);
            set(h_chansel_warnconfirm, 'Visible', 'on');
        else
            set(h_chansel_confirm, 'BackgroundColor', backgroundcolor_buttondefault);
            set(h_chansel_reset, 'BackgroundColor', backgroundcolor_buttondefault);
            set(h_chansel_warnconfirm, 'Visible', 'off');
        end
    end


    function f_chansel_commandentry(hObject, eventdata)
        set(h_hintbar, 'String', 'Selection commands are parsed one by one in sequential order. Leave blank to skip a step.');
        tmp_prompt = {'(Optional) Step 1. Enter a variable name in the base workspace to select these channels (example: channames):', '(Optional) Step 2. Regexp pattern to SELECT channels:', '(Optional) Step 3. Regexp pattern to DESELECT channels:', '(Optional) Step 4. Regexp pattern to SELECT channels:'};
        tmp_fieldsize = [1 60; 1 60; 1 60; 1 60];
        tmp_defaultinput = {'', '', '', ''};
        if evalin('base', 'exist(''signalviewer_PlottedChanNames'',''var'')') && isfirsttime_commandentry
            tmp_defaultinput{1} = 'signalviewer_PlottedChanNames';
        end
        tmp_answer = inputdlg(tmp_prompt, 'Select/Deselect Channels to Plot', tmp_fieldsize, tmp_defaultinput);
        if length(tmp_answer) ~= 4
            set(h_hintbar, 'String', 'Canceled selection. No change was made.');
            return
        end
        tmp_VARset = tmp_answer{1};
        tmp_REinclude = tmp_answer{2};
        tmp_REexclude = tmp_answer{3};
        tmp_REincludeagain = tmp_answer{4};
        selchan = get(h_chansel_list, 'Value');
        if ~isempty(tmp_VARset)
            %tmp_importedchannames = cell(0);
            tmp_whosinbase = evalin('base','who');
            if ismember(tmp_VARset,tmp_whosinbase)
                try
                    tmp_importedchannames = evalin('base', tmp_VARset);
                    selchan = chan2idx(ChanNames, tmp_importedchannames, 1);
                end
            end
        end
        for i = 1:length(ChanNames)
            if ~isempty(regexp(ChanNames{i}, tmp_REinclude, 'match', 'once'))
                selchan = union(selchan, i);
            end
        end
        for i = 1:length(ChanNames)
            if ~isempty(regexp(ChanNames{i}, tmp_REexclude, 'match', 'once'))
                selchan = setdiff(selchan, i);
            end
        end
        for i = 1:length(ChanNames)
            if ~isempty(regexp(ChanNames{i}, tmp_REincludeagain, 'match', 'once'))
                selchan = union(selchan, i);
            end
        end
        
        set(h_chansel_list, 'Value', selchan);
        isfirsttime_commandentry = false;
        f_chansel_confirm(hObject, eventdata);

    end


    function f_eventlabels_showhide(hObject, eventdata)
        textvisible = get(eventtexthand, 'Visible');
        linevisible = get(eventplothand, 'Visible');
        if any([textvisible{:}] == 1)
            set(eventtexthand, 'Visible', 'off');
        elseif any([linevisible{:}] == 1)
            set(eventplothand, 'Visible', 'off');
        else
            redraw();
        end
    end


    function f_chansel_confirm(hObject, eventdata)
        set(h_hintbar, 'String', 'Plot is updating'); drawnow;
        set(h_chansel_confirm, 'BackgroundColor', backgroundcolor_buttondefault);
        set(h_chansel_reset, 'BackgroundColor', backgroundcolor_buttondefault);
        set(h_chansel_warnconfirm, 'Visible', 'off');
        selchan = get(h_chansel_list, 'Value');
        Nsch = length(selchan);
        set(axehand, 'YTick', [-chansep*Nsch:chansep:-chansep], 'YTickLabel', fliplr(ChanNames(selchan)));

        % Po240528: Added these updaters before redraw 
        % (they won't recalculate unnecessarily anymore)
        notch_update();
        bandpass_update();
        envelope_update();
        render_update(); % this already includes redraw()
        
        %Automatically zoom out if less-than-full channels fill the screen
        if Nsch < (YLim(2)-YLim(1))/chansep
            f_yzoomout(hObject, eventdata);
        end

        same_as_car_channels = isequal(find(PerChannelFilterStates(1,:)), setdiff(selchan,chan2idx(ChanNames, nofiltchannames, 1)));
        warntext = '';
        if RerefFilter.state
            if ~same_as_car_channels
            warntext = 'WARNING: Common average reference was NOT updated to match the selected channels. Click [CAR] again to update.';
            set(h_car_switch, 'ForegroundColor', fontcolor_on3);
            else
                set(h_car_switch, 'ForegroundColor', fontcolor_on1);
            end
        end
        set(h_hintbar, 'String', ['Plot is updated. ' warntext]); drawnow;
    end

    function f_chansel_reset(hObject, eventdata)
        set(h_chansel_list, 'Value', selchan);
        set(h_chansel_confirm, 'BackgroundColor', backgroundcolor_buttondefault);
        set(h_chansel_reset, 'BackgroundColor', backgroundcolor_buttondefault);
        set(h_chansel_warnconfirm, 'Visible', 'off');
    end

    function f_psd_plot(hObject, eventdata)
        if ~ishandle(viewhand_psd)
            viewhand_psd = figure;
            psd_dblim_allchans = [inf -inf];
            psd_held_lpxx = [];
            psd_held_fxx = [];
            psd_held_Time1 = NaN;
            psd_held_Time2 = NaN;
            psd_held_chancolor = [1 1 1];
            psd_held_channame = '';
            psd_held_filttext = '';
            viewhand_psd_axe = axes('parent',viewhand_psd);
            set(viewhand_psd, 'KeyPressFcn', @f_fig_keypress);
            Signal_psd_source = Signal_psd_source*0;
            update_psd(1);
            set(viewhand_psd, 'MenuBar', 'none', 'Toolbar', 'figure', 'HandleVisibility', 'callback');
        else
            psd_dblim_allchans = [inf -inf];
            figure(viewhand_psd);
            set(viewhand_psd_axe, 'YLimMode', 'auto');
            psd_ylim = get(viewhand_psd_axe, 'YLim');
            psd_ylim = [max([psd_ylim(1) psd_ylim_restrict(1)]) min([psd_ylim(2) psd_ylim_restrict(2)])];
            if ~isequal(psd_ylim, get(viewhand_psd_axe, 'YLim'))
                set(viewhand_psd_axe, 'YLim', psd_ylim);
            end
            update_psd(1);
        end
        set(h_hintbar, 'String', 'PSD updated. Arrow keys to scroll left/right/up/down. Press h to hold a PSD for comparison (press again to unhold). Press c to select another channel. Press p to toggle PSD on/off.');
    end


    % function update_linemarker()
    %     if LinemarkerEnable
    %         set(plothand, 'Marker', ZoomedInMarker);
    %     else
    %         set(plothand, 'Marker', ZoomedOutMarker);
    %     end
    % end


    function update_cursorline(sel_type)
        
        if ~exist('sel_type','var') || isempty(sel_type)
            sel_type = selected_cursortype;
        end
        
        set([cursorlinehand plottext_cursor_hand plottext_lmin_hand plottext_lmax_hand plotpip1hand plotpip2hand plotpip3hand], 'Visible', 'off');
        
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
                YDATA_fullres = Signal_postenvelope(t1:BLIM:t2,selchan(ch));
                if ZscoreFilter.state
                    YDATA_fullres = nanzscore(YDATA_fullres)*ZscoreFilter.multiplier - chansep*ch;
                else
                    YDATA_fullres = YDATA_fullres - nanmedian(YDATA_fullres) - chansep*ch;
                end
                
                ypoint(ch) = interp1(Time, Signal_postenvelope(:,selchan(ch)), selected_timepoint);
                
                % Locally search for min and max in each channel
                lookcenter = find(Time >= selected_timepoint, 1, 'first');
                loclookcenter = find(Time(t1:BLIM:t2) >= selected_timepoint, 1, 'first');
                
                if isempty(loclookcenter) || loclookcenter == 1
                    return
                end
                
                lookradius = ceil((t2 - t1) * lookradius_fraction);
                
                tmpunit = 'µV';
                if ZscoreFilter.state && EnvelopeFilter.state
                    tmpunit = 'µV² (without zscore)';
                elseif ZscoreFilter.state
                    tmpunit = 'µV (without zscore)';
                elseif EnvelopeFilter.state
                    tmpunit = 'µV²';
                end
                
                xr = diff(get(axehand, 'XLim'));
                
                if sel_type == 1
                    ypointdisplayed(ch) = YDATA_fullres(loclookcenter);
                    set(plotpip1hand(ch), 'XData', selected_timepoint, 'YData', ypointdisplayed(ch), 'Visible', 'on');
                    set(plottext_cursor_hand(ch), 'String', sprintf('%s (%i,%i)\n%.6gs\n%.6g%s',ChanNames{selchan(ch)}, lookcenter, ch, selected_timepoint, ypoint(ch),tmpunit), 'Position', [selected_timepoint + xr/100, ypointdisplayed(ch), 0], 'Visible', 'on');
                elseif sel_type == 2
                    if lookcenter - lookradius < 1 || lookcenter + lookradius > Ntp
                        return
                    end
                    tmp = Signal_postenvelope(lookcenter-lookradius:lookcenter+lookradius,selchan(ch));
                    [ylocmin(ch), ind] = min(tmp);
                    xindmin(ch) = ind - 1 + loclookcenter - lookradius;
                    xlocmin(ch) = Time(t1-1+xindmin(ch));
                    if xindmin(ch) < 1 || xindmin(ch) > length(YDATA_fullres)
                        continue
                    end
                    set(plotpip2hand(ch), 'XData', xlocmin(ch), 'YData', YDATA_fullres(xindmin(ch)), 'Visible', 'on');
                    %2024-12-17 Added more info to the local min display
                    %setappdata(plottext_lmin_hand(ch), 'originalString', sprintf('%.4gs', xlocmin(ch)));
                    setappdata(plottext_lmin_hand(ch), 'originalString', sprintf('%s @%.6gs=%.6g%s',ChanNames{selchan(ch)}, xlocmin(ch), ylocmin(ch),tmpunit));

                    if size(SavedPointsTable,2) >= 2 && size(SavedPointsTable,1) >= 1
                        tcpair = cell2mat(SavedPointsTable(:,1:2));
                        L = xlocmin(ch)==tcpair(:,1) & selchan(ch)==tcpair(:,2);
                    else
                        L = false;
                    end
                    if any(L)
                        set(plottext_lmin_hand(ch), 'String', sprintf('%s\n%s',getappdata(plottext_lmin_hand(ch),'originalString'),'(SAVED)'), 'Position', [xlocmin(ch) + xr/100, YDATA_fullres(xindmin(ch)), 0], 'Visible', 'on');
                    else
                        set(plottext_lmin_hand(ch), 'String', getappdata(plottext_lmin_hand(ch),'originalString'), 'Position', [xlocmin(ch) + xr/100, YDATA_fullres(xindmin(ch)), 0], 'Visible', 'on');
                    end
                    setappdata(plottext_lmin_hand(ch), 'chanind', selchan(ch));
                    setappdata(plottext_lmin_hand(ch), 'ChanName', ChanNames{selchan(ch)});
                    setappdata(plottext_lmin_hand(ch), 'tvalue', xlocmin(ch));
                    setappdata(plottext_lmin_hand(ch), 'yvalue', ylocmin(ch));

                elseif sel_type == 3
                    if lookcenter - lookradius < 1 || lookcenter + lookradius > Ntp
                        return
                    end
                    tmp = Signal_postenvelope(lookcenter-lookradius:lookcenter+lookradius,selchan(ch));
                    [ylocmax(ch), ind] = max(tmp);
                    xindmax(ch) = ind - 1 + loclookcenter - lookradius;
                    xlocmax(ch) = Time(t1-1+xindmax(ch));
                    if xindmax(ch) < 1 || xindmax(ch) > length(YDATA_fullres)
                        continue
                    end
                    set(plotpip3hand(ch), 'XData', xlocmax(ch), 'YData', YDATA_fullres(xindmax(ch)), 'Visible', 'on');

                    %2024-12-17 Added more info to the local max display
                    %setappdata(plottext_lmax_hand(ch), 'originalString', sprintf('%.4gs', xlocmax(ch)));
                    setappdata(plottext_lmax_hand(ch), 'originalString', sprintf('%s @%.6gs=%.6g%s',ChanNames{selchan(ch)}, xlocmax(ch), ylocmax(ch),tmpunit));
                    
                    if size(SavedPointsTable,2) >= 2 && size(SavedPointsTable,1) >= 1
                        tcpair = cell2mat(SavedPointsTable(:,1:2));
                        L = xlocmax(ch)==tcpair(:,1) & selchan(ch)==tcpair(:,2);
                    else
                        L = false;
                    end
                    if any(L)
                        set(plottext_lmax_hand(ch), 'String', sprintf('%s\n%s',getappdata(plottext_lmax_hand(ch),'originalString'),'(SAVED)'), 'Position', [xlocmax(ch) + xr/100, YDATA_fullres(xindmax(ch)), 0], 'Visible', 'on');
                    else
                        set(plottext_lmax_hand(ch), 'String', getappdata(plottext_lmax_hand(ch),'originalString'), 'Position', [xlocmax(ch) + xr/100, YDATA_fullres(xindmax(ch)), 0], 'Visible', 'on');
                    end
                    setappdata(plottext_lmax_hand(ch), 'chanind', selchan(ch));
                    setappdata(plottext_lmax_hand(ch), 'ChanName', ChanNames{selchan(ch)});
                    setappdata(plottext_lmax_hand(ch), 'tvalue', xlocmax(ch));
                    setappdata(plottext_lmax_hand(ch), 'yvalue', ylocmax(ch));
                end
                
            end
        end
        
        
        
    end


    function update_infolabels()
        set([plotinfolabel1hand plotinfolabel2hand plotinfolabel3hand], 'Visible', 'off');
        infolabels_deviations = zeros(1,Nsch);
        if InfoLabelEnable
            for ch = Nsch:-1:1
                %set(plotinfolabel1hand(selchan(ch)), 'String', sprintf('[min %+.3g, max %+.3g, sd %.3g]', min(Signal_postenvelope(t1:t2,selchan(ch))), max(Signal_postenvelope(t1:t2,selchan(ch))), std(Signal_postenvelope(t1:t2,selchan(ch)))), 'Position', [XLim(1), -chansep*ch, 0], 'Visible', 'on');
                infolabels_deviations(ch) = std(Signal_postenvelope(t1:t2,selchan(ch))) / (sqrt(2)/2);
                set(plotinfolabel1hand(selchan(ch)), 'String', sprintf('%.1f', infolabels_deviations(ch) ), 'Position', [XLim(1), -chansep*ch, 0], 'Visible', 'on');
            end
        end
    end


    function update_psd(SWwipe)
        if ~exist('SWwipe','var') || isempty(SWwipe)
            SWwipe = 0;
        end
        persistent psd_dblim psd_lastchanname
        if isempty(psd_dblim) || SWwipe
            psd_dblim = [inf -inf];
            psd_lastchanname = '';
        end
        if ishandle(viewhand_psd)
            if ~strcmpi(psd_xscale,get(viewhand_psd_axe,'XScale'))
                set(viewhand_psd_axe, 'XScale', psd_xscale);
            end

            if ~isempty(selected_plothand) && ishandle(selected_plothand)
                psd_now_channame = getappdata(selected_plothand, 'channame');
                chanind = getappdata(selected_plothand, 'chanind');
                psd_now_chancolor = get(selected_plothand, 'Color');

                if ~ismember(chanind, selchan)
                    %Po240712 The PSD channel is not a selected channel.
                    %Force it to be.
                    selected_plothand = plothand(selchan(1));
                    psd_now_channame = getappdata(selected_plothand, 'channame');
                    chanind = getappdata(selected_plothand, 'chanind');
                    psd_now_chancolor = get(selected_plothand, 'Color');
                end
                
                tmp = Signal_postenvelope(t1:t2,chanind);
                tmp = tmp(isfinite(tmp));
                if size(Signal_psd_source,1) == size(tmp,1) && size(Signal_psd_source,2) == size(tmp,2) && norm(Signal_psd_source - tmp) == 0
                    % These are the same signal.
                    return
                end
                tmp1 = fighand;
                set(0, 'CurrentFigure', viewhand_psd);
                
                set(viewhand_psd_axe, 'FontUnits', 'pixel', 'FontSize', 16);
                Signal_psd_source = tmp;
                [psd_now_pxx, psd_now_fxx] = pwelch(Signal_psd_source, [], [], pwelch_nfft, Fs);
                psd_now_lpxx = 10*log10(psd_now_pxx);

                %2024-10-09: Plot the psd held in memory
                plotting_a_held_psd = 0;
                if ~isempty(psd_held_lpxx) && numel(psd_held_lpxx) == numel(psd_held_fxx)
                    plotting_a_held_psd = 1;
                    psd_held_hand = plot(psd_held_fxx, psd_held_lpxx, 'Color', psd_held_chancolor, 'LineWidth', psd_held_linewidth);
                    if strcmpi(psd_xscale, 'log')
                        set(gca, 'XScale', psd_xscale);
                    end
                    hold on
                    psd_now_hand = plot(psd_now_fxx, psd_now_lpxx, 'Color', psd_now_chancolor, 'LineWidth', psd_now_linewidth);
                    hold off
                else
                    psd_now_hand = plot(psd_now_fxx, psd_now_lpxx, 'Color', psd_now_chancolor, 'LineWidth', psd_now_linewidth); %Po240712: PSD line color matches main window's line color
                    if strcmpi(psd_xscale, 'log')
                        set(gca, 'XScale', psd_xscale);
                    end
                end
                xlabel('Frequency (Hz)');
                ylabel('PSD (dB/Hz)');
                fc1 = max(0,ceil(BandPassFilter.cutoff(1)*.9));
                fc2 = min(ceil(BandPassFilter.cutoff(2)*1.1),Fs/2);
                if ~BandPassFilter.state
                    fc1 = 0;
                end
                if ~BandPassFilter.state
                    fc2 = Fs/2;
                end
                if EnvelopeFilter.state
                    fc1 = 0;
                    fc2 = min(ceil(EnvelopeFilter.cutoff*1.1),Fs/2);
                end
                

                % Applicable filters are: ICA, CAR, Notch, BPF, Env
                tmp_af = {'ICA', 'CAR', 'Notch', 'BPF', 'Env'};
                if use_ltmat
                    tmp_af{1} = 'Custom';
                end
                tmp_ftr = false(1,5);
                if use_ltmat
                    if isscalar(selica) && selica > 1
                        tmp_ftr(1) = 1;
                    end
                else
                    if length(selica) ~= size(ica_A,2)
                        tmp_ftr(1) = 1;
                    end
                end
                if RerefFilter.state
                    tmp_ftr(2) = 1;
                end
                if NotchFilter.state
                    tmp_ftr(3) = 1;
                end
                if BandPassFilter.state
                    tmp_ftr(4) = 1;
                end
                if EnvelopeFilter.state
                    tmp_ftr(4) = 1;
                end
                
                if any(tmp_ftr)
                    psd_now_filttext = ['(after ' cell_to_string(tmp_af(tmp_ftr), ', ') ')'];
                else
                    psd_now_filttext = ['(no filter applied)'];
                end

                try  %#ok<*TRYNC>
                    set(viewhand_psd_axe, 'XLim', [fc1 fc2]);

                    irange = find(psd_now_fxx>=fc1 & psd_now_fxx<=fc2);
                    [ymax,imax] = max(psd_now_lpxx(irange));
                    fmax = psd_now_fxx(irange(imax));
                    viewhand_psd_peak = text(fmax, ymax, sprintf('PEAK\n%.3g Hz\n%.3g dB/Hz', fmax, ymax), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8); %#ok<NASGU>

                end
                tmp2 = get(viewhand_psd_axe, 'YLim');
                if ~isequal(psd_now_channame, psd_lastchanname)
                    psd_lastchanname = psd_now_channame;
                    psd_dblim = tmp2;
                    if ~psd_ylim_auto
                        % Keep the widest YLim when navigating across channels
                        if psd_dblim(1) < psd_dblim_allchans(1)
                            psd_dblim_allchans(1) = psd_dblim(1);
                        else
                            psd_dblim(1) = psd_dblim_allchans(1);
                        end
                        if psd_dblim(2) > psd_dblim_allchans(2)
                            psd_dblim_allchans(2) = psd_dblim(2);
                        else
                            psd_dblim(2) = psd_dblim_allchans(2);
                        end
                    end
                else
                    % Keep the widest YLim when navigating across time
                    if tmp2(1) < psd_dblim(1)
                        psd_dblim(1) = tmp2(1);
                    end
                    if tmp2(2) > psd_dblim(2)
                        psd_dblim(2) = tmp2(2);
                    end
                end
                try  %#ok<*TRYNC>
                    set(viewhand_psd_axe, 'YLim', psd_dblim);
                end
                psd_ylim = get(viewhand_psd_axe, 'YLim');
                psd_ylim = [max([psd_ylim(1) psd_ylim_restrict(1)]) min([psd_ylim(2) psd_ylim_restrict(2)])];
                if ~isequal(psd_ylim, get(viewhand_psd_axe, 'YLim'))
                    set(viewhand_psd_axe, 'YLim', psd_ylim);
                end

                psd_now_Time1 = Time(t1);
                psd_now_Time2 = Time(t2);

                if ~ismember(chanind, selchan)
                    title(sprintf('Welch PSD in %s, %g - %g s\n%s', psd_now_channame, psd_now_Time1, psd_now_Time2, 'WARNING: PSD channel is not plotted in the main window.'));
                else
                    title(sprintf('Welch PSD in %s, %g - %g s\n%s', psd_now_channame, psd_now_Time1, psd_now_Time2, ' '));
                end
                if plotting_a_held_psd
                    legend(sprintf('%s in %g - %g s %s', psd_held_channame, psd_held_Time1, psd_held_Time2, psd_held_filttext), sprintf('%s %s', psd_now_channame, psd_now_filttext));
                else
                    legend(sprintf('%s %s', psd_now_channame, psd_now_filttext));
                end
                
                set(viewhand_psd, 'Name', ['PSD from ' num2str(psd_now_Time1) ' to ' num2str(psd_now_Time2) ' s']);
                set(0, 'CurrentFigure', tmp1);
                figure(viewhand_psd);
            else
                set(viewhand_psd, 'Name', 'Select a channel first by clicking on its signal.');
            end

        end
    end


    function StitchSignalCell()
        if any(isnan(Signal(:)))
            stitch_mult = 0;
        end
        if SignalIsStitched
            tmp = Signal;
            averageepochduration = mean(diff(EventTimePoints,[],2)+1) / Fs;
            if stitch_mult > 0
                Fcut_minfreq = 2*stitch_mult/averageepochduration;
                if Fcut_minfreq < Fs / 10000
                    Fcut_minfreq = Fs / 10000;
                end
            end
            imax = size(EventTimePoints,1); % Don't confuse this with EventTimeStamps. Time points are sample indices, not seconds.
            for i = 1:imax
                % Do a detrend
                try
                    tmp(round(EventTimePoints(i,1):EventTimePoints(i,2)),:) = detrend(Signal(round(EventTimePoints(i,1):EventTimePoints(i,2)),:), 'linear', 'omitmissing');
                catch
                    % If failed, tmp is just Signal in that segment
                end
                % Do a high-pass if stitch_mult > 0
                if stitch_mult > 0
                    reflect_len = - averageepochduration * Fs / 2;
                    Fcut_stitch = Fcut_minfreq;
                    FO_stitch = 2;
                    try
                        tmp(round(EventTimePoints(i,1):EventTimePoints(i,2)),:) = freqfilter(tmp(round(EventTimePoints(i,1):EventTimePoints(i,2)),:), Fs, [Fcut_stitch, FO_stitch], 'high', 'butter', reflect_len);
                    catch
                        try
                            tmp(round(EventTimePoints(i,1):EventTimePoints(i,2)),:) = tmp(round(EventTimePoints(i,1):EventTimePoints(i,2)),:);
                        catch
                            % If failed, tmp is just Signal in that segment
                        end
                    end
                end
                if nanaround_stitch_samples > 0
                    try
                        tmp(round(EventTimePoints(i,1)) - 1 + [1:nanaround_stitch_samples],:) = NaN;
                        tmp(round(EventTimePoints(i,2)) + 1 + [-nanaround_stitch_samples:-1],:) = NaN;
                    catch
                        % If failed, tmp is just Signal in that segment
                    end
                end

            end
            if nanaround_stitch_samples > 0
                tmp = fix_nans_for_filtering(tmp, bridgenans_method);
            end
            Signal = tmp;

            %2024-11-07 Also adjust splitnotch_factors
            %splitnotch_factors = zeros(round(EventTimePoints(end,2)),1);
            %2024-12-17 Workaround NaN in the EventTimePoints
            splitnotch_factors = zeros(round(max(max(EventTimePoints))),1);
            for i = 1:size(EventTimePoints,1)
                tmp_a = round(EventTimePoints(i,1));
                tmp_b = round(EventTimePoints(i,2));
                try
                    splitnotch_factors(tmp_a:tmp_b) = 1 - 1./(1+exp((tmp_a:tmp_b).' - (tmp_a+tmp_b)/2));
                catch
                    % If failed, tmp is just Signal in that segment
                end
            end

            clear tmp
        end
    end


    function f_icasel_list(hObject, eventdata)
        f_icasel_confirm(hObject, eventdata);
    end


    function f_icasel_confirm(hObject, eventdata)
        if FilterBusy
            return;
        end
        
        if isempty(ica_W)
            % Use default ICA arguments
            FilterBusy = 1;
            set(h_icasel_confirm, 'Enable', 'off', 'String', 'Wait');
            set(h_bigtext, 'Visible', 'on', 'String', sprintf('ICA: Calculating ICs... Progress can be monitored in the command window.\nNote that ICA uses all channels, including those not currently selected for plotting.')); drawnow;
            try
                fprintf('Running FastICA on Signal with NaN removed and with these parameters: \n  stabilization = %s\n  maxNumIterations = %g\n  approach = %s\n  g = %s\n', fastica_stabilization, fastica_maxNumIterations, fastica_approach, fastica_g);
                tmp_finites = Signal(~any(Signal_nonfinite,2),:).';
                [~, ica_A, ica_W] = fastica(tmp_finites, 'stabilization', fastica_stabilization, 'maxNumIterations', fastica_maxNumIterations, 'approach', fastica_approach, 'g', fastica_g, 'interactivePCA', fastica_interactivePCA);
                clear tmp_finites
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

            if use_ltmat
                icachans = [{'orig'} string_to_cell(sprintf('ltmat %i,', 1:size(ltmat,3)),',')];
                set(h_icasel_list, 'Max', 1); % Disable multi-select
                set(h_icasel_title, 'String', 'Custom Transform', 'FontSize', 0.3); % Change title
                set(h_icasel_list, 'String', icachans);
                set(h_icasel_list, 'Value', 1);
                selica = 1;
                set(h_icasel_view_export, 'String', 'Export TM');
                set(h_icasel_reset, 'Enable', 'off');
                set(h_icasel_reset, 'Visible', 'off');
                set(h_ica_switch, 'String', 'ltmat', 'FontSize', 0.8, 'Tooltip', 'Bottom right list. Select original or one of the linear transformations'); %Change button
            else
                icachans = string_to_cell(num2str(1:mica, 'i%i '), ' ');
                set(h_icasel_list, 'String', icachans);
                set(h_icasel_list, 'Value', 1:mica);
                set(h_icasel_reset, 'Enable', 'on');
            end

            
            set(h_icasel_confirm, 'String', 'Mix', 'Enable', 'on');
            set(h_icasel_view_sources, 'Enable', 'on');
            set(h_icasel_view_mixmat, 'Enable', 'on');
            set(h_icasel_view_sepmat, 'Enable', 'on');
            set(h_icasel_view_export, 'Enable', 'on');
            ICA_Initialized = 1;
            FilterBusy = 0;
            set(h_bigtext, 'Visible', 'off', 'String', ''); drawnow;
            if ~isequal(eventdata, 'startup')
                f_icasel_view_sources([], []);
            end
        else
            tmp = get(h_icasel_list, 'Value');
            if length(tmp) == length(selica) && all(tmp == selica)
                % No change
                return;
            else
                selica = tmp;
            end
            if length(selica) == size(ica_A,2) || (use_ltmat && isequal(selica,1))
                % All ICs selected (or in the case of ltmat, "orig" is selected)
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Mixing new ICs...'); drawnow;
                Signal_postica = Signal;
                set(h_ica_switch, 'ForegroundColor', fontcolor_off, 'FontWeight', fontweight_off, 'Value', 0);
                FilterBusy = 0;
                set(h_bigtext, 'Visible', 'off', 'String', '');
            else
                FilterBusy = 1;
                set(h_bigtext, 'Visible', 'on', 'String', 'Mixing new ICs...'); drawnow;

                if use_ltmat
                    if ltmat_do_not_detrend
                        tmp = Signal.';
                    else
                        tmp = detrend(Signal, 'constant', 'omitmissing').';
                        
                    end
                    Signal_postica = (ltmat(:,:,selica-1)*tmp).';
                else
                    tmp = ica_A;
                    tmp(:,setdiff(1:size(ica_A,2),selica)) = 0;
                    Signal_postica = (tmp*ica_sig).';
                end
                set(h_ica_switch, 'ForegroundColor', fontcolor_on1, 'FontWeight', fontweight_on1, 'Value', 1);
                FilterBusy = 0;
                set(h_bigtext, 'Visible', 'off', 'String', '');
            end
            reref_update();
            notch_update();
            bandpass_update();
            envelope_update();
            render_update();
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
            ic_names = string_to_cell(num2str(1:size(ica_sig,1),['i%0' num2str(npad2) 'i,']),',');
            svopts.disable_ica = true;
            viewhand_ica_sig = signalviewer(ica_sig.', Fs, ic_names, svopts, [], [], 'This figure window displays the ICA sources. Close this window to return to the original time signals. Exclude and remix ICs in the original time signals, not in here.');
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
            set(viewhand_ica_A_ax, 'Position', ICA_AxesPosition);
            colormap jet
            colorbar
            set(viewhand_ica_A_ax, 'FontUnits', 'normalized', 'FontSize', 0.02);
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
            set(viewhand_ica_W_ax, 'Position', ICA_AxesPosition);
            colormap jet
            colorbar
            set(viewhand_ica_W_ax, 'FontUnits', 'normalized', 'FontSize', 0.02);
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
        assignin(export_workspace_dest, 'signalviewer_ica_A', ica_A);
        assignin(export_workspace_dest, 'signalviewer_ica_W', ica_W);
        assignin(export_workspace_dest, 'signalviewer_ica_sig', ica_sig);
        assignin(export_workspace_dest, 'signalviewer_selica', selica);
        assignin(export_workspace_dest, 'signalviewer_ltmat', ltmat);
        assignin(export_workspace_dest, 'signalviewer_ica_export_timestamp', now);
        set(h_hintbar, 'String', ['ICA variables exported to ' export_workspace_dest ' workspace.']);
    end


    function f_signal_export(hObject, eventdata)
        assignin(export_workspace_dest, 'signalviewer_Signal', Signal);
        assignin(export_workspace_dest, 'signalviewer_Signal_nonfinite', Signal_nonfinite);
        assignin(export_workspace_dest, 'signalviewer_Signal_postica', Signal_postica);
        assignin(export_workspace_dest, 'signalviewer_Signal_postreref', Signal_postreref);
        assignin(export_workspace_dest, 'signalviewer_Signal_postnotch', Signal_postnotch);
        assignin(export_workspace_dest, 'signalviewer_Signal_postbutter', Signal_postbutter);
        assignin(export_workspace_dest, 'signalviewer_Signal_postenvelope', Signal_postenvelope);
        assignin(export_workspace_dest, 'signalviewer_ChanNames', ChanNames);
        assignin(export_workspace_dest, 'signalviewer_selchan', selchan);
        assignin(export_workspace_dest, 'signalviewer_PlottedChanNames', ChanNames(selchan));
        assignin(export_workspace_dest, 'signalviewer_RerefChanNames', ChanNames(RerefFilter.chanidx));
        assignin(export_workspace_dest, 'signalviewer_SavedPointsTable', SavedPointsTable);
        if exist('EventTimeStamps','var')
            assignin(export_workspace_dest, 'signalviewer_EventTimeStamps', EventTimeStamps);
        else
            assignin(export_workspace_dest, 'signalviewer_EventTimeStamps', {});
        end

        assignin(export_workspace_dest, 'signalviewer_psd_held_lpxx', psd_held_lpxx);
        assignin(export_workspace_dest, 'signalviewer_psd_held_fxx', psd_held_fxx);
        assignin(export_workspace_dest, 'signalviewer_psd_held_Time1', psd_held_Time1);
        assignin(export_workspace_dest, 'signalviewer_psd_held_Time2', psd_held_Time2);
        assignin(export_workspace_dest, 'signalviewer_psd_held_chancolor', psd_held_chancolor);
        assignin(export_workspace_dest, 'signalviewer_psd_held_channame', psd_held_channame);
        assignin(export_workspace_dest, 'signalviewer_psd_held_filttext', psd_held_filttext);

        assignin(export_workspace_dest, 'signalviewer_psd_lpxx', psd_now_lpxx);
        assignin(export_workspace_dest, 'signalviewer_psd_pxx', psd_now_pxx);
        assignin(export_workspace_dest, 'signalviewer_psd_fxx', psd_now_fxx);
        assignin(export_workspace_dest, 'signalviewer_psd_Time1', psd_now_Time1);
        assignin(export_workspace_dest, 'signalviewer_psd_Time2', psd_now_Time2);
        assignin(export_workspace_dest, 'signalviewer_psd_chancolor', psd_now_chancolor);
        assignin(export_workspace_dest, 'signalviewer_psd_channame', psd_now_channame);
        assignin(export_workspace_dest, 'signalviewer_psd_filttext', psd_now_filttext);
        
        assignin(export_workspace_dest, 'signalviewer_fighand_Number', get(fighand,'Number'));
        assignin(export_workspace_dest, 'signalviewer_signalHashStr', signalHashStr);
        assignin(export_workspace_dest, 'signalviewer_signal_export_timestamp', now);
        set(h_hintbar, 'String', ['Signal variables exported to ' export_workspace_dest ' workspace.']);
    end


    function f_axesfont_inc(hObject, eventdata)
        AxesFontSize = AxesFontSize + .001;
        set(axehand, 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
        set([plottext_cursor_hand plottext_lmin_hand plottext_lmax_hand], 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
        set([plotinfolabel1hand plotinfolabel2hand plotinfolabel3hand], 'FontUnits', 'normalized', 'FontSize', AxesFontSize*InfoLabelFontSizeScale);
        set(h_hintbar, 'String', sprintf('Axes font size set to %g units', AxesFontSize));
    end

    function f_axesfont_dec(hObject, eventdata)
        if AxesFontSize > .002
            AxesFontSize = AxesFontSize - .001;
        end
        set(axehand, 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
        set([plottext_cursor_hand plottext_lmin_hand plottext_lmax_hand], 'FontUnits', 'normalized', 'FontSize', AxesFontSize);
        set([plotinfolabel1hand plotinfolabel2hand plotinfolabel3hand], 'FontUnits', 'normalized', 'FontSize', AxesFontSize*InfoLabelFontSizeScale);
        set(h_hintbar, 'String', sprintf('Axes font size set to %g units', AxesFontSize));
    end

    function f_linewidth_inc(hObject, eventdata)
        PlotLineWidth = PlotLineWidth + .5;
        psd_held_linewidth = psd_held_linewidth + .5;
        psd_now_linewidth = psd_now_linewidth + .5;
        set(plothand, 'LineWidth', PlotLineWidth);
        set(h_hintbar, 'String', sprintf('Set line width to %g units', PlotLineWidth));
        if ishandle(psd_held_hand)
            set(psd_held_hand, 'LineWidth', psd_held_linewidth);
        end
        if ishandle(psd_now_hand)
            set(psd_now_hand, 'LineWidth', psd_now_linewidth);
        end
    end

    function f_linewidth_dec(hObject, eventdata)
        if PlotLineWidth > .5
            PlotLineWidth = PlotLineWidth - .5;
            psd_held_linewidth = psd_held_linewidth - .5;
            psd_now_linewidth = psd_now_linewidth - .5;
        end
        set(plothand, 'LineWidth', PlotLineWidth);
        set(h_hintbar, 'String', sprintf('Set line width to %g units', PlotLineWidth));
        if ishandle(psd_held_hand)
            set(psd_held_hand, 'LineWidth', psd_held_linewidth);
        end
        if ishandle(psd_now_hand)
            set(psd_now_hand, 'LineWidth', psd_now_linewidth);
        end
    end

    function f_eventfont_inc(hObject, eventdata)
        if EventEnable
            EventFontSize = EventFontSize + 1;
            set(eventtexthand, 'FontSize', EventFontSize);
            set(h_hintbar, 'String', sprintf('Event font size set to %g pixels', EventFontSize));
        end
    end

    function f_eventfont_dec(hObject, eventdata)
        if EventEnable
            if EventFontSize > 2
                EventFontSize = EventFontSize - 1;
            end
            set(eventtexthand, 'FontSize', EventFontSize);
            set(h_hintbar, 'String', sprintf('Event font size set to %g pixels', EventFontSize));
        end
    end

    function f_windowhsize_inc(hObject, eventdata)
        %AxesPosition = [0.0500    0.0600    0.86    0.93];
        AxesPosition = AxesPosition + [-0.01 0 0.01 0];
        set(axehand, 'Position', AxesPosition);
        set(h_hintbar, 'String', sprintf('Window horizontal size set to %g fraction', AxesPosition(3)));
    end

    function f_windowhsize_dec(hObject, eventdata)
        %AxesPosition = [0.0500    0.0600    0.86    0.93];
        if AxesPosition(3) > 0
            AxesPosition = AxesPosition - [-0.01 0 0.01 0];
        end
        set(axehand, 'Position', AxesPosition);
        set(h_hintbar, 'String', sprintf('Window horizontal size set to %g fraction', AxesPosition(3)));
    end

    function f_plothand_buttondown(hObject, eventdata)
        selected_plothand = hObject;
        update_psd();
        
        acp = get(axehand, 'CurrentPoint');
        selected_timepoint = acp(1);

        %Po240611: Snap to the closest actual time point.
        [~,ind] = min(abs(selected_timepoint - Time));
        selected_timepoint = Time(ind);

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

        %Po240621: Autofit to the current visible area, not the entire plot
        %Fs_int = round(Fs);
        %tmpr = 1+Fs_int:size(Signal_postenvelope,1)-Fs_int+1;
        tmpr = find(Time >= XLim(1),1):find(Time <= XLim(2),1,'last');

        nchanvisible = diff(YLim) / chansep;
        tmpy = - ((YLim + chansep/2) / chansep - 1);
        tmpcl = linspace(tmpy(2),tmpy(1),nchanvisible+1 );
        tmpcl = tmpcl(1:end-1);
        tmpcl = tmpcl(tmpcl>=1 & tmpcl<=Nsch);
        tmpc = selchan(tmpcl);

        if isempty(tmpr)
            tmpr = 1:size(Signal_postenvelope,1);
        end
        if isempty(tmpc)
            tmpc = selchan;
        end

        
        if ZscoreFilter.state
            tmp = quantile(  max(nanzscore(Signal_postenvelope(tmpr,tmpc))*ZscoreFilter.multiplier) - min(nanzscore(Signal_postenvelope(tmpr,tmpc))*ZscoreFilter.multiplier), 0.5 );
        else
            tmp = quantile(max(Signal_postenvelope(tmpr,tmpc)) - min(Signal_postenvelope(tmpr,tmpc)),0.5);
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
while containsValidString(s)
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



function y = fix_nans_for_filtering(y, bridgenans_method)
% 20180514 Remove NaNs before processing
ny = isnan(y);

if nnz(ny) > 0
    y = bridge_nans(y, bridgenans_method, '');
    
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
[h, c] = hist(s); %#ok<*HIST>
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

function list = stagger_odds_evens(list)
n = length(list);
list = list([1:2:n, 2:2:n]);
end

