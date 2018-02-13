% This function extracts the timing boundaries from an analog physical
% signal such as an electro-goniometer or a pressure foot switch.
% The filter settings assume these signals came from volitional movements.
% If multiple signals are provided, each column will be processed
% separately, and the output will be cells
%
% x = signals. Each column is a channel of time series signal.
% Fs_x = sampling frequency of x, in Hz.
%
% Optional parameters:
%
% SWinfo: Information level. 0 = silent. 1 = text. 2 = text and plot.
%
% SWskipmanualinstruction: Skip showing the popup that explains the manual
%                          threshold selection
%
% SWautoaccept: Automatically accepts threshold selection after loading
%
% JitTol: Jitter tolerance when counting the number of actions for plateau
%
% MinLen: remove an action if shorter than this duration (unit: seconds)
%
% MaxLen: remove an action if longer than this duration (unit: seconds)
%
% MergeLen: merge two adjacent actions if separated by less than this
%                 (time unit: seconds)
%
% MinHeight: remove an action if vertically shorter than this height
%
% MaxHeight: remove an action if vertically taller than this height
%
% Ymin: Minimum for threshold finding. Default is the global minimum
%
% Ymax: Maximum for threshold finding. Default is the global maximum
%
% ManualThreshold: 1 = manually draw thresholds. 
%                  2 = prompt even if there is a stored threshold profile.
%                  [] (empty) = Completely disallows ManualThreshold
%
% SWloadManualCoords: 1 = load 'ManualCoords' from parent workspace.
%                         setting this will also set ManualThreshold to 1
%                         when ManualCoords is successfully loaded
%                     0 = do not load (default)
%                     %-1 = also do not load persistent coordinates
%
% SWloadSignalOverlay: 0 or 1 = load 'SignalOverlay' from parent workspace
%                               if available
%                      -1 = never load
% 
% Fcutoff:         (N by 2 matrix, N = number of data sets)
%                  First column is low cutoff, second column is high cutoff
%                  NaN to use default value
%                  If Fcutoff(Dataset,:) = [0 0], filter is disabled
%
% Squared:         (Array of same size as the number of data sets, each
%                   element describing whether to square the signal on the
%                   corresponding data set)
%                   0 = Let the program decide whether to square a signal
%                  -1 = Will not square a signal
%                   1 = Will always square a signal
%
% Flipped:         (Array of same size as the number of data sets, each
%                   element describing whether to flip the signal on the
%                   corresponding data set)
%                   0 = Let the program decide whether to flip a signal
%                   upside down
%                  -1 = Will not flip signal
%                   1 = Will always flip signal
%
% ThresLevel       If a stable range is detected, from the bottom, thing.
%
% Comment          Comment in the title

% 20140325: Removed the -1 in Actions = [Bounds(:,1) Bounds(:,2)]./Fs_x;

function [Bounds_sec, StableWidth, thres_level, ManualCoords] = getbounds2 (x, Fs_x, varargin)
%%
persistent Persist

SWinfo = [];
JitTol = [];
MinLen = [];
MaxLen = [];
MergeLen = [];
MinHeight = [];
MaxHeight = [];
SWloadManualCoords = [];
SWskipmanualinstruction = 0;
SWautoaccept = 0;
SWedgetreatment = 1;
SWedgeaddleft = 1;
SWedgeaddright = 1;

if mod(length(varargin),2)
    error('Properties list must be in pairs, ie. property1 name, property1 value, ...');
end

if nargin > 2
    tmpgindisp = cell(1,floor((nargin-2)/2));
    for i = 1:2:nargin-2
        tmpgindisp{(i+1)/2} = [varargin{i} ' = [' num2str(varargin{i+1}) ']'];
        if strcmpi(varargin{i},'Comment')
            Comment = varargin{i+1};
        else
            eval([varargin{i} ' = [' num2str(varargin{i+1}) '];']);
        end
    end
end

% Verbosity level
if ~exist('SWinfo','var') || isempty(SWinfo)
    SWinfo = 0;
end

if SWinfo > 0 && nargin > 2 && exist('tmpgindisp','var')
    disp('getbounds2 command line arguments: ');
    for i = 1:length(tmpgindisp)
        disp(tmpgindisp{i});
    end
    disp(' ');
end

if ~exist('SWloadManualCoords','var') || isempty(SWloadManualCoords)
    SWloadManualCoords = 0;
end

loadedManualCoords = 0;
if SWloadManualCoords > 0
    try
        ManualCoords = evalin('caller', 'ManualCoords');
        if SWinfo > 0
            disp('loaded ManualCoords from caller workspace');
            loadedManualCoords = 1;
            ManualThreshold = max(ManualThreshold,1);
        end
    catch %#ok<CTCH>
        if SWinfo > 0
            disp('cannot load ManualCoords from caller workspace');
        end
    end
end

% Check input dimension
if length(x) ~= size(x,1)
    x = x.';
    if SWinfo >= 1
        fprintf('Transposed the signal matrix.\n');
    end
end

if ~exist('Fs_x','var')
    error('Sampling frequency of the signals must be specified.');
end

if ~exist('MinLen','var') || isempty(MinLen)
    MinLen = 0;
end

if ~exist('MaxLen','var') || isempty(MaxLen)
    MaxLen = inf;
end

if ~exist('MinHeight','var') || isempty(MinHeight)
    MinHeight = -inf;
end

if ~exist('MaxHeight','var') || isempty(MaxHeight)
    MaxHeight = inf;
end

if ~exist('MergeLen','var') || isempty(MergeLen)
    MergeLen = 0.050;
end

if islogical(x) || isinteger(x)
    isLogicalSignal = 1;
    if SWinfo >= 1
        fprintf('Signal is of logical or integer type.\n');
    end
    x = double(x);
else
    isLogicalSignal = 0;
end

Ndataset = size(x,2);
SignalIsFlipped = zeros(1,Ndataset);
SignalIsAbsed = zeros(1,Ndataset);

if isLogicalSignal
    % Logical signals cannot be flipped, filtered, or squared
    Flipped = -1 * ones(1,Ndataset);
    Squared = -1 * ones(1,Ndataset);
    Fcutoff = nan(Ndataset, 2);
end

if ~exist('Fcutoff','var') || isempty(Fcutoff) 
    Fcutoff = nan(Ndataset, 2);
end

if ~exist('Squared','var') || isempty(Squared) 
    Squared = zeros(1, Ndataset);
end

if ~exist('Flipped','var') || isempty(Flipped) 
    Flipped = zeros(1, Ndataset);
end

if ~exist('ThresLevel','var') || isempty(ThresLevel) %#ok<NODEF>
    ThresLevel = 0.5;
end

OverlayEnabled = 0;

if ~exist('SWloadSignalOverlay','var') || SWloadSignalOverlay >= 0
    try
        SignalOverlay = evalin('caller', 'SignalOverlay');
        if min(size(SignalOverlay) == size(x))
            OverlayEnabled = 1;
        end
    end
end




if length(Flipped) < Ndataset
    Flipped = Flipped(:).';
    Flipped(Ndataset-length(Flipped)+1:Ndataset) = 0;
end

if length(Squared) < Ndataset
    Squared = Squared(:).';
    Squared(Ndataset-length(Squared)+1:Ndataset) = 0;
end

if size(Fcutoff,1) < Ndataset
    Fcutoff(Ndataset-size(Fcutoff,1)+1:Ndataset,:) = NaN;
end

if size(Fcutoff,2) < 2
    Fcutoff(:,2) = NaN;
end

if Ndataset == 1
    Bounds_sec = [0 0];
    StableWidth = zeros(1,Ndataset);
    thres_level = zeros(1,Ndataset);
else
    Bounds_sec = cell(1,Ndataset);
    StableWidth = zeros(1,Ndataset);
    thres_level = zeros(1,Ndataset);
end


if ~exist('Comment','var') || isempty(Comment) %#ok<NODEF>
    Comment = '';
end


Fs_default = 256;

% Check if signal is audio (actual decision is inside loop)
if ~isLogicalSignal
    if Fs_x > Fs_default
        HFspec = signalpower(x, Fs_x, [0 128; 128 Fs_x/2]);
    else
        HFspec = zeros(2,Ndataset);
    end
else
    HFspec = zeros(2,Ndataset);
end

%20140328: Shifted start time to 0
t = (0:size(x,1)-1).'/Fs_x;
%x_filtered = resample(freqfilter(x, Fs_x, [0.035 35], 'pass', 'butter'), Fs_default, Fs_x);
%Fs_x = Fs_default;
%x_filtered = freqfilter(x, Fs_x, [0.035 35], 'pass', 'butter');
%x_filtered = freqfilter(x.^2, Fs_x, [0.01 35], 'pass', 'butter');

% 20130828: Feature disabled
% if ~isLogicalSignal
%     xsig = norm(x);
%     if isfield(Persist,'xsig') && Persist.xsig == xsig && SWloadManualCoords ~= -1
%         if ~loadedManualCoords && isfield(Persist,'ManualCoords')
%             ManualCoords = Persist.ManualCoords;
%         end
%     else
%         Persist.xsig = xsig;
%     end
% end

x2 = x;
x2_midxrange = ceil(size(x2,1)/100*[5 95]);
x2_orig = x2;
x_orig = x;


for Dataset = 1:Ndataset
    
    if SWinfo >= 1
        fprintf('Now processing time series %i ..\n', Dataset);
    end
    
    %x2 = x_filtered(:,Dataset);
    
    if ~isLogicalSignal
        
        % Is signal mostly high frequency?
        isHF = HFspec(2,Dataset) > HFspec(1,Dataset)*2;
        if isHF && SWinfo >= 1
            fprintf('Signal is mostly high frequency (probably audio signal).\n');
        end
        mustHF = Squared(Dataset) == 1;
        mustnotHF = Squared(Dataset) == -1;
        if isHF && mustnotHF && SWinfo >= 1
            fprintf('HF mode is disallowed.\n');
        end
        
        if mustHF || (~mustnotHF && isHF)
            if SWinfo >= 1
                fprintf('Switching to HF mode.\n');
            end
            if ~exist('JitTol','var') || isempty(JitTol)
                JitTol = 2;
            end
            if ~isfinite(Fcutoff(Dataset,1))
                Fcutoff(Dataset,1) = 0; %#ok<AGROW>
            end
            if ~isfinite(Fcutoff(Dataset,2))
                Fcutoff(Dataset,2) = 10; %#ok<AGROW>
            end
            
            [fcuts, ftype] = localfunc_determine_filtertype (Fcutoff, Dataset, Fs_x);
            if ~isempty(fcuts)
                x2(:,Dataset) = freqfilter(x(:,Dataset).^2, Fs_x, fcuts, ftype, 'butter');
                x2_orig = x2;
            else
                x2(:,Dataset) = x(:,Dataset);
                x2_orig = x2;
            end

        else
            if ~exist('JitTol','var') || isempty(JitTol)
                %JitTol = 15;
            end
            if ~isfinite(Fcutoff(Dataset,1))
                %Fcutoff(Dataset,1) = 0.035; %#ok<AGROW>
                Fcutoff(Dataset,1) = 0.01; %#ok<AGROW>
            end
            if ~isfinite(Fcutoff(Dataset,2))
                Fcutoff(Dataset,2) = 35; %#ok<AGROW>
            end
            
            
            [fcuts, ftype] = localfunc_determine_filtertype (Fcutoff, Dataset, Fs_x);
            if ~isempty(fcuts)
                x2(:,Dataset) = freqfilter(x(:,Dataset), Fs_x, fcuts, ftype, 'butter');
                x2_orig = x2;
            else
                x2(:,Dataset) = x(:,Dataset);
                x2_orig = x2;
            end
            
        end
    else
        % Signal is logical or integer, there is not much to do.
        x2(:,Dataset) = x(:,Dataset);
        x2_orig = x2;
        
        if ~exist('JitTol','var') || isempty(JitTol)
            %JitTol = 0;
        end
    end
    
    %t2 = (1:size(x2,1)).'/Fs_x;
    %20171011: Changed the start time to 0
    t2 = (0:size(x2,1)-1).'/Fs_x;
    
    me = median(x2(:,Dataset));
    
    isF = mean(x2(:,Dataset)) < me;
    if isF && SWinfo >= 1
        fprintf('Signal mean is less than median.\n');
    end
    mustF = Flipped(Dataset) == 1;
    mustnotF = Flipped(Dataset) == -1;
    if isF && mustnotF && SWinfo >= 1
        fprintf('Flipping signal is disallowed.\n');
    end
    
    if mustF || (~mustnotF && isF)
        if SWinfo >= 1
            fprintf('Flipped the signal upside-down.\n');
        end
        x2(:,Dataset) = -x2(:,Dataset);
        x(:,Dataset) = -x(:,Dataset);
        x2_orig = x2;
        x_orig = x;
        SignalIsFlipped(Dataset) = 1 - SignalIsFlipped(Dataset);
        me = -me;
    end
    
    %%
    
    
    if ~isLogicalSignal
        if exist('Ymin','var') && ~isempty(Ymin)
            x2min = Ymin;
        else
            x2min = min(x2(:,Dataset));
        end
        
        if exist('Ymax','var') && ~isempty(Ymax)
            x2max = Ymax;
        else
            x2max = max(x2(:,Dataset));
        end
        
        % You can increase the number of points to improve resolution at the
        % expense of time consumption
        xtest = linspace(x2min,x2max,201);
        
        TBounds = cell(1,length(xtest));
        parfor i = 1:length(xtest)
            [rise, fall] = localfunc_findedges (x2(:,Dataset) > xtest(i), SWedgetreatment, SWedgeaddleft, SWedgeaddright, size(x2,1));
            if ~isempty(rise) && length(rise) == length(fall)
                TBounds{i} = [rise fall];
            end
        end
        
        bn = zeros(1,length(TBounds));
        for i = 1:length(TBounds)
            bn(i) = length(TBounds{i});
        end
        
        bnd = [0 abs(diff(bn))];
        if isempty(JitTol)
            JitTol = mean(bnd);
        end
        if SWinfo >= 1
            fprintf('Using Jitter tolerance = %g\n', JitTol);
        end
        [rise, fall] = localfunc_findedges((bnd <= JitTol).');
        
    else
        %20140116: Added edge treatments
        [rise, fall] = localfunc_findedges(x2(:,Dataset)>0, SWedgetreatment, SWedgeaddleft, SWedgeaddright, size(x2,1));
    end
    
    %Has no first rising edge
    while ~isempty(rise) && ~isempty(fall) && rise(1) > fall(1)
        fall = fall(2:end);
    end
    
    %Has no last falling edge
    while ~isempty(fall) && ~isempty(rise) && fall(end) < rise(end)
        rise = rise(1:end-1);
    end
    
    if ~isLogicalSignal
        if ~isempty(rise) && length(rise) == length(fall)
            PBounds = [rise fall];
            [Y I] = max(diff(PBounds,[],2)); %#ok<ASGLU>
            xtestplateaubound = xtest(PBounds(I,:));
            StableWidth(Dataset) = diff(xtestplateaubound);
            
            % Use the threshold that is ThresLevel up from the bottom to determine boundary
            Ind = round(PBounds(I,1)+diff(PBounds(I,:))*ThresLevel);
            Bounds = TBounds{Ind};
            thres_level(Dataset) = xtest(Ind);
            
        else
            PBounds = [];
            Bounds = [];
            xtestplateaubound = [];
            StableWidth(Dataset) = 0;
        end
        if SWinfo >= 1
            if ~isempty(PBounds)
                fprintf('Plateau region found between [%.4g, %.4g]\n', xtest(PBounds(I,:)));
                fprintf('Will use %.4g as the threshold to detect actions.\n', xtest(Ind));
            else
                fprintf('No plateau detected.\n');
            end
        end
    else
        Bounds = [rise fall];
        xtestplateaubound = [];
        StableWidth(Dataset) = 0;
    end
    
    
    if SWinfo >= 1
        fprintf('Identified %i actions (an action is an instance of flexion, foot step, etc.).\n', size(Bounds,1));
    end
    
    if exist('ManualThreshold','var') && ~isempty(ManualThreshold) && ManualThreshold > 0
        if SWinfo >= 1
            fprintf('\n');
            fprintf('Manual threshold selection activated\n');
        end
        MTgood = 0;
        if exist('ManualCoords','var')
            if iscell(ManualCoords) && length(ManualCoords) == Ndataset && size(ManualCoords{Dataset},2) == 2
                MTgood = 1;
            elseif ~iscell(ManualCoords) && size(ManualCoords,2) == 2
                MTgood = 1;
            end
        end
        if MTgood
            if SWinfo >= 1
                fprintf('A previously saved manual threshold profile has been loaded.\n');
            end
        else
            ManualCoords = {};
        end
        if ManualThreshold > 1 || ~MTgood
            if SWinfo >= 1
                fprintf('A new manual threshold selection is in progress ..\n');
            end
            
            mcfh = figure(1370000000);
            set(mcfh, 'NumberTitle', 'off', 'Name', 'getbounds2');
            clf
            if ~SWautoaccept
                manx2ph = plot(t2, x2(:,Dataset), 'Color', [0 0 0]);
                if OverlayEnabled
                    hold on
                    overlayph = plot(t2, SignalOverlay(:,Dataset), '-', 'Color', [0 0.7 0], 'LineWidth', 2);
                    hold off
                end
                scrsz = get(0, 'ScreenSize');
                scrsz = scrsz(1,:);
                newpos = scrsz + [0 60 -3 -100];
                set(gcf, 'Position', newpos, 'MenuBar', 'none', 'Toolbar', 'none');
                %set(gcf, 'Position', [1 100 1600 800], 'MenuBar', 'none', 'Toolbar', 'none')
                set(gca, 'Position', [0.05 0.07 0.92 0.86])
                
                HELPMSG = [
                    'Please manually insert nodes to draw threshold lines. Signals above the threshold lines are logical high.' 10 10 ...
                    'UP/DOWN arrows: Zoom' 10 ...
                    'LEFT/RIGHT arrows: Pane' 10 ...
                    'Left mouse click: Insert new node at crosshair' 10 ...
                    'Right mouse click: Delete node closest to crosshair' 10 10 ...
                    'The pre-existing nodes are from automated threshold finding or from previously saved session (You can right click to remove them).' 10 10 ...
                    'Press ''k'' to clear all selections.' 10 ...
                    'Press ''f'' to flip the signal upside-down.' 10 ...
                    'Press ''b'' to take/untake absolute value of the signal.' 10 ...
                    'Press ''t'' to filter the signal.' 10 ...
                    'Press ''p'' to change parameters.' 10 ...
                    'Press ''e'' to confirm your current selections.' 10 ...
                    'Press ''c'' to print the coordinates at cursor.' 10 ...
                    'Press ''x'' to exit the program with an empty output.'
                    ];
                
                if ~SWskipmanualinstruction
                    localfunc_modalmsgbox(HELPMSG);
                end
                
                TITLE = '';
                titlehand = title(TITLE);
                XLABEL = 'Left/Right Click=Add/Delete node.  Arrow keys=Zoom/Pane.  ''e''=Done.  ''h''=Help message.';
                set(titlehand, 'FontSize', 12, 'FontWeight', 'bold');
                xlabelhand = xlabel(XLABEL);
                set(xlabelhand, 'FontSize', 12, 'FontWeight', 'bold');
            else
                text(0.5,0.5, sprintf('Processing channel %i of %i', Dataset, Ndataset), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
            ZinMat = ( [0.5 0.5; 0.5 0.5] + [0.25 -0.25; -0.25 0.25] );
            ZoutMat = ( [0.5 0.5; 0.5 0.5] + [1 -1; -1 1] );
            PrightMat = ( [-1 1; -1 1]/5 + [1 0; 0 1] );
            PleftMat = ( [1 -1; 1 -1]/5 + [1 0; 0 1] );
            InitialXLim = get(gca, 'XLim');
            
            
            if length(ManualCoords) < Dataset
                % Start off with detected ones.
                % ManualCoords{Dataset} = [t2(1) thres_level(Dataset); t2(end) thres_level(Dataset)];
                %20130824: Start off with blank.
                ManualCoords{Dataset} = zeros(0,2);
            end
            
            hold on
            if ~SWautoaccept
                manh = plot(InitialXLim - 1,[0 0],'s-r');
                %hold off
                set(manh, 'XData', ManualCoords{Dataset}(:,1), 'YData', ManualCoords{Dataset}(:,2));
                set(gca, 'XLim', InitialXLim);
                tbox1 = uicontrol(gcf, 'Style', 'text', 'Units', 'pixels', 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 11, 'BackgroundColor', [0.8 1 0.95], 'FontName', 'Consolas');
                %set(tbox1, 'Position', [0.838 0.835 0.162 0.165]);
                tbox1w = 248.184;
                tbox1h = 106.06;
                gcfpos = get(gcf, 'Position');
                set(tbox1, 'Position', [gcfpos(3)-tbox1w+1, gcfpos(4)-tbox1h+1, tbox1w, tbox1h]);
            end
            trying = 1;
            mb = -1;
            % Main loop
            while trying
                if ~SWautoaccept && ~ishandle(manh(1))
                    break
                end
                if mb == -1
                    mb = 117;
                elseif SWautoaccept
                    mx = 0;
                    my = 0;
                    mb = 101;
                else
                    [mx, my, mb] = ginput(1);
                end
                RunInstantUpdate = 0;
                RunRedraw = 0;
                if mb
                    switch mb
                        
                        case 1
                            % LMB
                            if ~any(ManualCoords{Dataset}(:,1) == mx)
                                ManualCoords{Dataset} = [ManualCoords{Dataset};[mx my]];
                                ManualCoords{Dataset} = unique(ManualCoords{Dataset},'rows');
                                ManualCoords{Dataset} = ManualCoords{Dataset}((ManualCoords{Dataset}(:,1)>=InitialXLim(1)) & (ManualCoords{Dataset}(:,1)<=InitialXLim(2)),:);
                                set(manh, 'XData', ManualCoords{Dataset}(:,1), 'YData', ManualCoords{Dataset}(:,2));
                                RunInstantUpdate = 1;
                                RunRedraw = 1;
                            end
                        case 3
                            % RMB
                            if ~isempty(ManualCoords{Dataset})
                                [Y I] = min((ManualCoords{Dataset}(:,1) - mx).^2 + (ManualCoords{Dataset}(:,2) - my).^2); %#ok<ASGLU>
                                ManualCoords{Dataset} = ManualCoords{Dataset}([1:I-1,I+1:end],:);
                                set(manh, 'XData', ManualCoords{Dataset}(:,1), 'YData', ManualCoords{Dataset}(:,2));
                                RunInstantUpdate = 1;
                                RunRedraw = 1;
                            end
                        case 30
                            % Up arrow key. Zoom in
                            XL = get(gca, 'XLim');
                            NewXrange = (XL(2) - XL(1)) * 0.50;
                            set(gca, 'XLim', mx + [-NewXrange/2, NewXrange/2]);
                            %set(gca, 'XLim', ZinMat * get(gca, 'XLim')');
                            XL = get(gca, 'XLim');
                            tmp_i = t2>=XL(1) & t2<=XL(2);
                            m1 = double(min(x2(tmp_i,Dataset)));
                            m2 = double(max(x2(tmp_i,Dataset)));
                            if OverlayEnabled
                                m1o = double(min(SignalOverlay(tmp_i,Dataset)));
                                m1 = min(m1,m1o);
                                m2o = double(max(SignalOverlay(tmp_i,Dataset)));
                                m2 = max(m2,m2o);
                            end
                            ymin = m1-(m2-m1)*0.01;
                            ymax = m2+(m2-m1)*0.01;
                            %keyboard
                            try
                                set(gca, 'YLim', [ymin ymax]);
                            catch
                                set(gca, 'YLimMode', 'auto');
                            end
                            RunRedraw = 1;
                            
                        case 31
                            % Down arrow key. Zoom out
                            set(gca, 'XLim', ZoutMat * get(gca, 'XLim')');
                            if diff(get(gca,'XLim')) > diff(InitialXLim)
                                set(gca, 'XLim', InitialXLim);
                            end
                            XL = get(gca, 'XLim');
                            tmp_i = t2>=XL(1) & t2<=XL(2);
                            m1 = double(min(x2(tmp_i,Dataset)));
                            m2 = double(max(x2(tmp_i,Dataset)));
                            if OverlayEnabled
                                m1o = double(min(SignalOverlay(tmp_i,Dataset)));
                                m1 = min(m1,m1o);
                                m2o = double(max(SignalOverlay(tmp_i,Dataset)));
                                m2 = max(m2,m2o);
                            end                            
                            ymin = m1-(m2-m1)*0.01;
                            ymax = m2+(m2-m1)*0.01;
                            try
                                set(gca, 'YLim', [ymin ymax]);
                            catch
                                set(gca, 'YLimMode', 'auto');
                            end
                            RunRedraw = 1;
                        case 28
                            % Left arrow key. Pane left
                            if [1 0]*PleftMat*get(gca, 'XLim')' >= InitialXLim(1)
                                set(gca, 'XLim', PleftMat * get(gca, 'XLim')');
                            end
                            XL = get(gca, 'XLim');
                            tmp_i = t2>=XL(1) & t2<=XL(2);
                            m1 = double(min(x2(tmp_i,Dataset)));
                            m2 = double(max(x2(tmp_i,Dataset)));
                            if OverlayEnabled
                                m1o = double(min(SignalOverlay(tmp_i,Dataset)));
                                m1 = min(m1,m1o);
                                m2o = double(max(SignalOverlay(tmp_i,Dataset)));
                                m2 = max(m2,m2o);
                            end                            
                            ymin = m1-(m2-m1)*0.01;
                            ymax = m2+(m2-m1)*0.01;
                            try
                                set(gca, 'YLim', [ymin ymax]);
                            catch
                                set(gca, 'YLimMode', 'auto');
                            end
                            RunRedraw = 1;
                        case 29
                            % Right arrow key. Pane right
                            if [0 1]*PrightMat*get(gca, 'XLim')' <= InitialXLim(2)
                                set(gca, 'XLim', PrightMat * get(gca, 'XLim')');
                            end
                            XL = get(gca, 'XLim');
                            tmp_i = t2>=XL(1) & t2<=XL(2);
                            m1 = double(min(x2(tmp_i,Dataset)));
                            m2 = double(max(x2(tmp_i,Dataset)));
                            if OverlayEnabled
                                m1o = double(min(SignalOverlay(tmp_i,Dataset)));
                                m1 = min(m1,m1o);
                                m2o = double(max(SignalOverlay(tmp_i,Dataset)));
                                m2 = max(m2,m2o);
                            end                            
                            ymin = m1-(m2-m1)*0.01;
                            ymax = m2+(m2-m1)*0.01;
                            try
                                set(gca, 'YLim', [ymin ymax]);
                            catch
                                set(gca, 'YLimMode', 'auto');
                            end
                            RunRedraw = 1;
                        case 101
                            % The 'e' key finalizes the manual selection.
                            trying = 0;
%                         case 105
%                             % The 'i' key resets the manual bounds
%                             ManualCoords{Dataset} = [t2(1) thres_level(Dataset); t2(end) thres_level(Dataset)];
%                             set(manh, 'XData', ManualCoords{Dataset}(:,1), 'YData', ManualCoords{Dataset}(:,2));
%                             RunInstantUpdate = 1;
%                             RunRedraw = 1;
                        case 107
                            % The 'k' key erases all manual bounds
                            ManualCoords{Dataset} = zeros(0,2);
                            set(manh, 'XData', ManualCoords{Dataset}(:,1), 'YData', ManualCoords{Dataset}(:,2));
                            RunInstantUpdate = 1;
                            RunRedraw = 1;
                        case 120
                            % The 'x' key erases all manual bounds and then
                            % exits
                            ManualCoords{Dataset} = zeros(0,2);
                            set(manh, 'XData', ManualCoords{Dataset}(:,1), 'YData', ManualCoords{Dataset}(:,2));
                            trying = 0;
                            RunInstantUpdate = 1;
                            RunRedraw = 1;
                        case 99
                            % The 'c' key prints the current XY coordinate
                            if exist('cursorprint','var') && ~isempty(cursorprint) && ishandle(cursorprint)
                                delete(cursorprint);
                                cursorprint = [];
                            end
                            if exist('cursormark','var') && ~isempty(cursormark) && ishandle(cursormark)
                                delete(cursormark);
                                cursormark = [];
                            end
                            cursormark = plot(mx, my, 'x', 'Color', [0.9 0.9 0], 'LineWidth', 1.5, 'MarkerSize', 16);
                            cursorprint = text(mx, my, ['CURSOR' 10 '(' num2str(mx, '%.4g') ', ' num2str(my, '%.4g') ')'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0 0.6 0]);
                            drawnow
                        case 102
                            % The 'f' key flips the signal
                            x2(:,Dataset) = -x2(:,Dataset);  x2_orig(:,Dataset) = -x2_orig(:,Dataset);
                            x(:,Dataset) = -x(:,Dataset);  x_orig(:,Dataset) = -x_orig(:,Dataset);
                            SignalIsFlipped(Dataset) = 1 - SignalIsFlipped(Dataset);
                            set(manx2ph, 'YData', x2(:,Dataset));
                            RunInstantUpdate = 1;
                            RunRedraw = 1;
                            %localfunc_modalmsgbox('Signal has been flipped upside-down. Press ''k'' to delete your selections, then click new selections.');
                        case 98
                            % The 'b' key abs the signal
                            if ~SignalIsAbsed(Dataset)
                                x2(:,Dataset) = abs(x2(:,Dataset));
                            else
                                x2 = x2_orig;
                            end
                            SignalIsAbsed(Dataset) = 1 - SignalIsAbsed(Dataset);
                            set(manx2ph, 'YData', x2(:,Dataset));
                            RunInstantUpdate = 1;
                            RunRedraw = 1;
                            %localfunc_modalmsgbox('Signal has been absed. Press ''k'' to delete your selections, then click new selections.');
                        case 116
                            % The 't' key filters the signal
                            DEF1 = num2str(Fcutoff(Dataset,1));
                            DEF2 = num2str(Fcutoff(Dataset,2));
                            dlgans = inputdlg({'High-pass filter cutoff (Hz, Enter 0 to disable)', 'Low-pass filter cutoff (Hz, Enter 0 to disable)'}, 'Name' , 1, {DEF1, DEF2});
                            if ~isempty(dlgans)
                                Fcutoff(Dataset,1) = str2double(dlgans{1});
                                Fcutoff(Dataset,2) = str2double(dlgans{2});
                                [fcuts, ftype] = localfunc_determine_filtertype (Fcutoff, Dataset, Fs_x);
                                if ~isempty(fcuts)
                                    try
                                        x2(:,Dataset) = freqfilter(x(:,Dataset), Fs_x, fcuts, ftype, 'butter');
                                        x2_orig = x2;
                                    catch
                                        Fcutoff(Dataset,:) = [0 0];
                                        x2(:,Dataset) = x(:,Dataset);
                                        x2_orig = x2;
                                        localfunc_modalmsgbox('Frequency filtering encountered an error. Filters are currently disabled to avoid the error. Please try again.');
                                    end
                                else
                                    x2(:,Dataset) = x(:,Dataset);
                                    x2_orig = x2;
                                end
                                set(manx2ph, 'YData', x2(:,Dataset));
                                RunInstantUpdate = 1;
                                RunRedraw = 1;
                                %localfunc_modalmsgbox('Signal has been filtered. Press ''k'' to delete your selections, then click new selections.');
                            end
                        case 112
                            % The 'p' key changes parameters
                            DEF_MERG = num2str(MergeLen);
                            DEF_MINL = num2str(MinLen);
                            DEF_MAXL = num2str(MaxLen);
                            DEF_MINH = num2str(MinHeight);
                            DEF_MAXH = num2str(MaxHeight);
                            DEF_JIT = num2str(JitTol);
                            dlgans = inputdlg({'JitTol=', 'MinLen=', 'MaxLen=', 'MinHeight=', 'MaxHeight=', 'MergeLen='}, 'Parameters' , 1, {DEF_JIT, DEF_MINL, DEF_MAXL, DEF_MINH, DEF_MAXH, DEF_MERG});
                            if ~isempty(dlgans)
                                ANS_JIT = str2double(dlgans{1});
                                ANS_MINL = str2double(dlgans{2});
                                ANS_MAXL = str2double(dlgans{3});
                                ANS_MINH = str2double(dlgans{4});
                                ANS_MAXH = str2double(dlgans{5});
                                ANS_MERG = str2double(dlgans{6});
                                if ~isnan(ANS_JIT), JitTol = ANS_JIT; end
                                if ~isnan(ANS_MINL), MinLen = ANS_MINL; end
                                if ~isnan(ANS_MAXL), MaxLen = ANS_MAXL; end
                                if ~isnan(ANS_MINH), MinHeight = ANS_MINH; end
                                if ~isnan(ANS_MAXH), MaxHeight = ANS_MAXH; end
                                if ~isnan(ANS_MERG), MergeLen = ANS_MERG; end
                                RunInstantUpdate = 1;
                                RunRedraw = 1;
                                %localfunc_modalmsgbox('Parameters have been changed.');
                            end
                        case 104
                            % The 'h' key displays help
                            localfunc_modalmsgbox(HELPMSG);
                        case 117
                            % The 'u' key refreshes
                            RunInstantUpdate = 1;
                            RunRedraw = 1;
                    end
                    
                    %% 20130506: Instant update
                    if RunInstantUpdate
                        if ~isempty(ManualCoords{Dataset})
                            ManualCoords{Dataset} = unique(ManualCoords{Dataset},'rows');
                            ManualCoords{Dataset} = ManualCoords{Dataset}((ManualCoords{Dataset}(:,1)>=InitialXLim(1)) & (ManualCoords{Dataset}(:,1)<=InitialXLim(2)),:);
                            tmp = [t2(1) ManualCoords{Dataset}(1,2); t2(end) ManualCoords{Dataset}(end,2)];
                            ManualCoords{Dataset} = [ManualCoords{Dataset}; tmp];
                            ManualCoords{Dataset} = unique(ManualCoords{Dataset},'rows');
                            [~, tmp_i] = unique(ManualCoords{Dataset}(:,1));
                            ManualCoords{Dataset} = ManualCoords{Dataset}(tmp_i,:);
                            MTthresline = interp1(ManualCoords{Dataset}(:,1),ManualCoords{Dataset}(:,2),t2,'linear');
                            [rise, fall] = localfunc_findedges (x2(:,Dataset) > MTthresline, SWedgetreatment, SWedgeaddleft, SWedgeaddright, size(x2,1));
                            if ~isempty(rise) && length(rise) == length(fall)
                                Bounds = [rise fall];
                            else
                                Bounds = [];
                            end
                        else
                            Bounds = [];
                        end
                        
                        %StableWidth(Dataset) = 0;
                        %thres_level(Dataset) = NaN;
                        %xtestplateaubound = [];
                        
                        
                        if ~isempty(Bounds)
                            %20171011: Changed to 0 start index
                            Actions = ([Bounds(:,1) Bounds(:,2)] - 1)./Fs_x;
                            % Process MinLen, MaxLen, MinHeight, MaxHeight
                            peak = zeros(size(Bounds,1),1);
                            peakloc = peak;
                            for k = 1:size(Bounds,1)
                                [tmp1, tmp2] = max(x2(Bounds(k,1):Bounds(k,2),Dataset));
                                if ~isempty(tmp1)
                                    peak(k) = tmp1;
                                    peakloc(k) = tmp2;
                                else
                                    peak(k) = NaN;
                                    peakloc(k) = NaN;
                                end
                                %20171011: Changed to 0 start index
                                peakloc(k) = (peakloc(k)-1) / Fs_x;
                            end
                            passed_minlen = diff(Actions,[],2) >= MinLen;
                            passed_maxlen = diff(Actions,[],2) <= MaxLen;
                            passed_minhei = peak >= MinHeight;
                            passed_maxhei = peak <= MaxHeight;
                            passed_all = passed_minlen & passed_maxlen & passed_minhei & passed_maxhei;
                            if SWinfo >= 1
                                fprintf('Deleted %i actions that are outside the duration and height limits.\n', nnz(~passed_all));
                            end
                            Actions = Actions(passed_all,:);
                            
                            % Process MergeLen
                            Downs = [Actions(1:end-1,2) Actions(2:end,1)];
                            to_merge = diff(Downs,[],2) < MergeLen;
                            Actions(to_merge,2) = NaN;
                            Actions(logical([0;to_merge(1:end)]),1) = NaN;
                            
                            Actions = [Actions(~isnan(Actions(:,1)),1) Actions(~isnan(Actions(:,2)),2)];
                            
                            if SWinfo >= 1
                                fprintf('Merged %i actions that are shorter than merge duration limit.\n', nnz(~to_merge));
                            end
                            
                        else
                            Actions = [];
                        end
                        if Ndataset == 1
                            Bounds_sec = Actions;
                        else
                            Bounds_sec{Dataset} = Actions;
                        end
                        
                        if ~isempty(Actions)
                            %20171011: Changed to 0 start index
                            %Bounds = round(Actions*Fs_x+ones(size(Actions,1),1)*[1 0]);
                            Bounds = round(Actions*Fs_x+1);
                        else
                            Bounds = [];
                        end
                        
                        
                        if exist('phand','var') && ~isempty(phand) && ishandle(phand(1))
                            delete(phand);
                        end
                        BoundColor = 0.9 - 0.15*jet(size(Bounds,1));
                        BoundColor = equalspaceorder(BoundColor);
                        
                        %                     for k = 1:size(Actions,1)
                        %                         max(x2(Actions(k,1):Actions(k,2),Dataset))
                        %                     end
                        
                        %m1 = double(min(x2(Fs_x*5:Fs_x*(t2(end)-5),Dataset)));
                        %m2 =
                        %double(max(x2(Fs_x*5:Fs_x*(t2(end)-5),Dataset)));p
                        
                        if ~SWautoaccept
                            m1 = double(min(x2(x2_midxrange(1):x2_midxrange(2),Dataset)));
                            m2 = double(max(x2(x2_midxrange(1):x2_midxrange(2),Dataset)));
                            ymin = m1-(m2-m1)*0.01;
                            ymax = m2+(m2-m1)*0.01;
                            phand = patch_bounds([], Actions, [ymin ymax], BoundColor, BoundColor*0.9);
                            tmp_c = get(gca, 'Children');
                            if OverlayEnabled
                                set(gca, 'Children', [overlayph; manx2ph; setdiff(tmp_c, [overlayph manx2ph])]);
                            else
                                set(gca, 'Children', [manx2ph; setdiff(tmp_c, manx2ph)]);
                            end
                            if SignalIsFlipped(Dataset), tmp_Flip = 'Yes'; else tmp_Flip = 'No'; end;
                            tboxcell = {
                                sprintf('Dataset# %05i  Flipped= %s', Dataset, tmp_Flip);
                                sprintf('  Fcutoff=[%.3g %.3g] Hz', Fcutoff(Dataset,1:2));
                                sprintf('   JitTol= %.3g', JitTol);
                                sprintf('LengthLim=[%.3g %.3g]', MinLen, MaxLen);
                                sprintf('HeightLim=[%.3g %.3g]', MinHeight, MaxHeight);
                                sprintf(' MergeLen= %.3g', MergeLen);
                                };
                            set(tbox1, 'String', tboxcell);
                            gcfpos = get(gcf, 'Position');
                            set(tbox1, 'Position', [gcfpos(3)-tbox1w+1, gcfpos(4)-tbox1h+1, tbox1w, tbox1h]);
                        end
                    end
                    
                    if ~SWautoaccept && RunRedraw
                        if ~exist('Actions','var')
                            Actions = [];
                        end
                        if exist('texthand0','var') && ~isempty(texthand0) && ishandle(texthand0(1))
                            delete(texthand0);
                            texthand0 = [];
                        end
                        
                        if exist('texthand1','var') && ~isempty(texthand1) && ishandle(texthand1(1))
                            delete(texthand1);
                            texthand1 = [];
                        end
                        if exist('texthand2','var') && ~isempty(texthand2) && ishandle(texthand2(1))
                            delete(texthand2);
                            texthand2 = [];
                        end
                        if exist('topedgehand','var') && ~isempty(topedgehand) && ishandle(topedgehand(1))
                            try delete(topedgehand); end;
                            topedgehand = [];
                        end
                        
                        info_Nact = size(Actions,1);
                        if info_Nact > 0
                            widths = diff(Actions,[],2);
                            info_Width = [mean(widths), std(widths), min(widths), median(widths), max(widths)];
                            peak = zeros(size(Bounds,1),1);
                            peakloc = peak;
                            for k = 1:size(Bounds,1)
                                [tmp1, tmp2] = max(x2(Bounds(k,1):Bounds(k,2),Dataset));
                                if ~isempty(tmp1)
                                    peak(k) = tmp1;
                                    peakloc(k) = tmp2;
                                else
                                    peak(k) = NaN;
                                    peakloc(k) = NaN;
                                end
                                %20171011: Changed to 0 start index
                                peakloc(k) = (peakloc(k)-1) / Fs_x;
                            end
                            info_Height = [nanmean(peak), nanstd(peak), nanmin(peak), nanmedian(peak), nanmax(peak)];
                            
                            ITI = diff(Actions(:,1),[],1);
                            info_ITI = [mean(ITI), std(ITI), min(ITI), median(ITI), max(ITI)];
                            
                            %zoomingy = get(gca, 'YLim');
                            %zoomingx = get(gca, 'XLim');
                            %texthand0(1) = text(mean(zoomingx), zoomingy(1), 'WIDTH OF EACH ACTION', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 36, 'FontWeight', 'bold', 'Color', [0.8 0.8 0.8]);
                            %texthand0(2) = text(mean(zoomingx), zoomingy(2), 'HEIGHT OF EACH ACTION', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 36, 'FontWeight', 'bold', 'Color', [0.8 0.8 0.8]);
                            
                            for k = size(Actions,1):-1:1
                                if ~isnan(peakloc(k))
                                    try
                                        topedgehand(k) = plot(Actions(k,1)+peakloc(k), peak(k), 'mo', 'LineWidth', 2, 'MarkerSize', 15);
                                    catch
                                        keyboard
                                    end
                                    %ctr = mean( [Actions(k,1), Actions(k,2)] );
                                    %texthand1(k) = text(ctr, zoomingy(1), ['#' num2str(k) 10 num2str(widths(k), '%.3g')], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold');
                                    %texthand2(k) = text(ctr, zoomingy(2), ['#' num2str(k) 10 num2str(peak(k), '%.3g')], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9, 'FontWeight', 'bold');
                                end
                            end
                            
                            
                        else
                            info_Width = [0 0 0 0 0];
                            info_Height = [0 0 0 0 0];
                            info_ITI = [0 0 0 0 0];
                        end
                        
                        TITLE = sprintf('%s: %i actions. (mean,sd,min,median,max) Width = (%.3g, %.3g, %.3g, %.3g, %.3g) s. \nHeight = (%.3g, %.3g, %.3g, %.3g, %.3g). ITI = (%.3g, %.3g, %.3g, %.3g, %.3g) s.', Comment, info_Nact, info_Width, info_Height, info_ITI);
                        set(titlehand, 'String', TITLE, 'Interpreter', 'none');
                    end
                    if ~SWautoaccept
                        drawnow
                    end
                    %%
                    
                end
                
            end
            
            if ~SWautoaccept
                close(mcfh);
                drawnow
            end
        else
            InitialXLim = t2([1 end]);
        end
        
        if ~isempty(ManualCoords{Dataset})
            ManualCoords{Dataset} = unique(ManualCoords{Dataset},'rows');
            ManualCoords{Dataset} = ManualCoords{Dataset}((ManualCoords{Dataset}(:,1)>=InitialXLim(1)) & (ManualCoords{Dataset}(:,1)<=InitialXLim(2)),:);
            tmp = [t2(1) ManualCoords{Dataset}(1,2); t2(end) ManualCoords{Dataset}(end,2)];
            ManualCoords{Dataset} = [ManualCoords{Dataset}; tmp];
            ManualCoords{Dataset} = unique(ManualCoords{Dataset},'rows');
            [~, tmp_i] = unique(ManualCoords{Dataset}(:,1));
            ManualCoords{Dataset} = ManualCoords{Dataset}(tmp_i,:);
            MTthresline = interp1(ManualCoords{Dataset}(:,1),ManualCoords{Dataset}(:,2),t2,'linear');
            [rise, fall] = localfunc_findedges (x2(:,Dataset) > MTthresline, SWedgetreatment, SWedgeaddleft, SWedgeaddright, size(x2,1));
            if ~isempty(rise) && length(rise) == length(fall)
                Bounds = [rise fall];
            else
                Bounds = [];
            end
        else
            Bounds = [];
        end
        
        StableWidth(Dataset) = 0;
        thres_level(Dataset) = NaN;
        xtestplateaubound = [];
        
        if SWinfo >= 1
            fprintf('Manual threshold selection completed. Voided automated selections in favor of manual selection.\n');
        end
        
    else
        ManualCoords{Dataset} = [];
    end
    
    if size(Bounds,2) == 2
        %20171011: Changed to 0 start index
        Actions = ([Bounds(:,1) Bounds(:,2)]-1)./Fs_x;
        % Process MinLen, MaxLen, MinHeight, MaxHeight
        peak = zeros(size(Bounds,1),1);
        peakloc = peak;
        for k = 1:size(Bounds,1)
            [tmp1, tmp2] = max(x2(Bounds(k,1):Bounds(k,2),Dataset));
            if ~isempty(tmp1)
                peak(k) = tmp1;
                peakloc(k) = tmp2;
            else
                peak(k) = NaN;
                peakloc(k) = NaN;
            end
            %20171011: Changed to 0 start index
            peakloc(k) = (peakloc(k)-1) / Fs_x;
        end
        passed_minlen = diff(Actions,[],2) >= MinLen;
        passed_maxlen = diff(Actions,[],2) <= MaxLen;
        passed_minhei = peak >= MinHeight;
        passed_maxhei = peak <= MaxHeight;
        passed_all = passed_minlen & passed_maxlen & passed_minhei & passed_maxhei;
        if SWinfo >= 1
            fprintf('Deleted %i actions that are outside the duration and height limits.\n', nnz(~passed_all));
        end
        Actions = Actions(passed_all,:);
        
        % Process MergeLen
        Downs = [Actions(1:end-1,2) Actions(2:end,1)];
        to_merge = diff(Downs,[],2) < MergeLen;
        Actions(to_merge,2) = NaN;
        Actions(logical([0;to_merge(1:end)]),1) = NaN;
        
        Actions = [Actions(~isnan(Actions(:,1)),1) Actions(~isnan(Actions(:,2)),2)];
        
        if SWinfo >= 1
            fprintf('Merged %i actions that are shorter than merge duration limit.\n', nnz(~to_merge));
        end

    else
        Actions = [];
    end
    if Ndataset == 1
        Bounds_sec = Actions;
    else
        Bounds_sec{Dataset} = Actions;
    end
    
    if ~isempty(Actions)
        %Bounds = Actions*Fs_x+ones(size(Actions,1),1)*[1 0];
        %20171011: Changed to 0 start index
        Bounds = Actions*Fs_x+1;
    else
        Bounds = [];
    end
    
    if SWinfo >= 1
        fprintf('\n\n');
        if ~isempty(Bounds)
            B.among.median = median(diff(Actions,[],2));
            B.among.sd = 1.4826*mad(diff(Actions,[],2),1);
            B.between.median = median(diff(Actions(:,1),[],1));
            B.between.sd = 1.4826*mad(diff(Actions(:,1),[],1),1);
            groupfindid = diff(Bounds(:,1),[],1) > (B.between.median + B.among.median);
            fprintf('Duration of each action = %.4g (%.4g) sec.\n', B.among.median/Fs_x, B.among.sd);
            fprintf('Time between the onsets of adjacent actions = %.4g (%.4g) sec.\n', B.between.median, B.between.sd);
        else
            groupfindid = 0;
            fprintf('No actions detected.\n');
        end
        fprintf('\n');
        if nnz(groupfindid) > 1
            B.group.median = median(diff(Actions(groupfindid,1),[],1));
            %[C LAGS] = autocorr( groupfindid );
            [C LAGS] = xcorr( double(groupfindid), double(groupfindid), 'coeff' );
            C = C(LAGS>=0);
            LAGS = LAGS(LAGS>=0);
            tmp = LAGS(C>mean(C));
            if length(tmp) > 1
                B.actions.per.group.median = tmp(2);
                B.number.groups = nnz(groupfindid);
                fprintf('Median time between the onsets of adjacent groups of actions = %g sec.\n', B.group.median);
                fprintf('Median number of actions in a group = %g\n', B.actions.per.group.median);
                fprintf('Total number of groups = %g\n', B.number.groups);
            else
                fprintf('No grouping detected.\n');
            end
        else
            fprintf('No grouping detected.\n');
        end
        fprintf('\n');
    end
    
    %%
    if SWinfo >= 2
        BoundColor = 0.90-0.20*rand(size(Bounds,1),3);
        
        figure(13701+(Dataset-1)*2);
        clf
        subplot(2,1,1);
        hold on
        if isLogicalSignal
            SigProcType = 'digital';
        else
            SigProcType = 'filtered';
        end
        
        title(['x(' num2str(Dataset) ',:) ' SigProcType ' signal (black) and threshold line (red)']);
        %patch_bounds([], t([1 end]).', xtestplateaubound, [0.4 1 0.4], [0.4 1 0.4]);
        if ~isempty(xtestplateaubound)
            patch_bounds([], Actions, xtestplateaubound, BoundColor, BoundColor);
        end
        if ~isnan(thres_level(Dataset))
            patch_bounds([], Actions, thres_level(Dataset)*[1 1], [1 0 0], [1 0 0]);
        end
        %plot(t, x(:,Dataset) - mean(x(:,Dataset)), 'Color', [0.8 0.8 0.8]);
        
        if exist('MTthresline','var')
            plot(t2, MTthresline, 'Color', [1 0 0]);
        end
        
        plot(t2, x2(:,Dataset), 'Color', [0 0 0]);
        plot(t([1 end]), double(me)*[1 1], 'm', 'Visible', 'off');
        %m1 = double(min(x2(Fs_x*5:Fs_x*(t2(end)-5),Dataset)));
        %m2 = double(max(x2(Fs_x*5:Fs_x*(t2(end)-5),Dataset)));
        m1 = double(min(x2(x2_midxrange(1):x2_midxrange(2),Dataset)));
        m2 = double(max(x2(x2_midxrange(1):x2_midxrange(2),Dataset)));
        ylim([m1-(m2-m1)*0.01 m2+(m2-m1)*0.01]);
        Xlim = [0 t(end)];
        xlim(Xlim);
        %figure(13702+(Dataset-1)*2);
        %clf
        subplot(2,1,2);
        hold on
        title(['x(' num2str(Dataset) ',:)' ' raw signal (black) and detected actions groups (random colors)']);
        
        patch_bounds([], Actions, [-1000 1000], BoundColor, BoundColor);
        plot(t, x(:,Dataset), 'Color', [0 0 0]);
        %m1 = double(min(x(Fs_x*5:Fs_x*(t2(end)-5),Dataset)));
        %m2 = double(max(x(Fs_x*5:Fs_x*(t2(end)-5),Dataset)));
        m1 = double(min(x2(x2_midxrange(1):x2_midxrange(2),Dataset)));
        m2 = double(max(x2(x2_midxrange(1):x2_midxrange(2),Dataset)));
        ylim([m1-(m2-m1)*0.01 m2+(m2-m1)*0.01]);
        xlim(Xlim);
        scrsz = get(0, 'ScreenSize');
        scrsz = scrsz(1,:);
        newpos = scrsz + [0 60 -3 -100];
        set(gcf, 'Position', newpos, 'MenuBar', 'none', 'Toolbar', 'none', 'NumberTitle', 'off', 'Name', 'getbounds2');
    end
    
end

if SWautoaccept
    close(mcfh);
end

% 20130828: Feature disabled
%Persist.ManualCoords = ManualCoords;


function localfunc_modalmsgbox (msg)
h = msgbox(msg, 'getbounds2 - manual mode', 'modal');
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 11 );
set(h,'Position', get(h,'Position') + [0 0 50 0]);
uiwait(h);



function [rise, fall] = localfunc_findedges (signal, SWedgetreatment, SWedgeaddleft, SWedgeaddright, Datasize)
% Simplified version of findedges.m
signal = (signal~=0);
for ch = 1:size(signal,2)
    sd = diff([0;signal(:,ch)]);
    rise = find(sd == 1);
    fall = find(sd == -1) - 1;
end

if nargin < 2
    SWedgetreatment = 0;
    SWedgeaddleft = 0;
    SWedgeaddright = 0;
end

if SWedgetreatment
    %Has no first rising edge
    while ~isempty(rise) && ~isempty(fall) && rise(1) > fall(1)
        if SWedgeaddleft
            rise = [1;rise];
        else
            fall = fall(2:end);
        end
    end
    
    %Has no last falling edge
    while ~isempty(fall) && ~isempty(rise) && fall(end) < rise(end)
        if SWedgeaddright
            fall = [fall;Datasize];
        else
            rise = rise(1:end-1);
        end
    end
end




function [fcuts, ftype] = localfunc_determine_filtertype (Fcutoff, Dataset, Fs_x)
if Fcutoff(Dataset,1) > 0 && Fcutoff(Dataset,2) < Fs_x/2 && Fcutoff(Dataset,2) > Fcutoff(Dataset,1)
    fcuts = Fcutoff(Dataset,1:2);
    ftype = 'pass';
elseif Fcutoff(Dataset,1) > 0
    fcuts = Fcutoff(Dataset,1);
    ftype = 'high';
elseif Fcutoff(Dataset,2) < Fs_x/2 && Fcutoff(Dataset,2) > Fcutoff(Dataset,1)
    fcuts = Fcutoff(Dataset,2);
    ftype = 'low';
else
    fcuts = [];
    ftype = 'pass';
end
