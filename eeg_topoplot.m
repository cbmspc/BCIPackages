function [xi, yi, zi, fhand] = eeg_topoplot (Filter, ChanNames, time, FullHeadElectrodes, OnlyShowElectrodes, SWnocontour, ManualClim, SWnocolorbar, FontSize)
% Plot topography (skull map).
% Filter is the coefficients in (chan,time) format
% ChanNames is a cell of strings
% time is the time points you wish to include in zi
% FullHeadElectrodes (optional) overrides the default 61-channel cap
% electrode names
% OnlyShowElectrodes (optional) shows the labels of only these electrodes
%    (shortcut) Enter 19 for 19-channel cap, 61 for 61-channel cap.
%
% EXAMPLE: 
% [xi, yi, zi, fhand] = eeg_topoplot (randn(19,1), {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4' ,'F8', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2'}, 1, {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4' ,'F8', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2'}, 19, 1, [], [], 24);
%

%Re-create the appearance of full head coverage by filling missing channels
%with zeros
if ~exist('FullHeadElectrodes','var') || isempty(FullHeadElectrodes)
    FullHeadElectrodes = {'AF3', 'AF4', 'AF7', 'AF8', 'AFz', ...
        'C1', 'C2', 'C3', 'C4', 'C5', 'C6', ...
        'CP1', 'CP2', 'CP3', 'CP4', 'CP5', 'CP6', 'CPz', 'Cz', ...
        'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', ...
        'FC1', 'FC2', 'FC3', 'FC4', 'FC5', 'FC6', 'FCz', 'FT7', 'FT8', ...
        'Fp1', 'Fp2', 'Fpz', 'Fz', 'O1', 'O2', 'Oz', ...
        'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', ...
        'PO3', 'PO4', 'PO7', 'PO8', 'POz', 'Pz', ...
        'T7', 'T8', 'TP7', 'TP8'};
end
if ~exist('OnlyShowElectrodes','var') || isempty(OnlyShowElectrodes)
    OnlyShowElectrodes = FullHeadElectrodes;
elseif ~iscell(OnlyShowElectrodes) && numel(OnlyShowElectrodes) == 1
    switch OnlyShowElectrodes
        case 19
            OnlyShowElectrodes = 19;
        case 20
            OnlyShowElectrodes = 19;
        case 61
            OnlyShowElectrodes = 61;
        case 62
            OnlyShowElectrodes = 61;
        case 63
            OnlyShowElectrodes = 61;
        case 64
            OnlyShowElectrodes = 61;
        otherwise
            OnlyShowElectrodes = 61;
    end
    if OnlyShowElectrodes == 19
        OnlyShowElectrodes = {'Fp1', 'Fp2', ...
            'F7', 'F3', 'Fz', 'F4' ,'F8', ...
            'T7', 'C3', 'Cz', 'C4', 'T8', ...
            'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2'};
    elseif OnlyShowElectrodes == 61
        OnlyShowElectrodes = FullHeadElectrodes;
    end
end
MissingElectrodes = FullHeadElectrodes(setdiff(1:length(FullHeadElectrodes),chan2idx(FullHeadElectrodes,ChanNames)));
Filter = [Filter; zeros(length(MissingElectrodes),size(Filter,2))];
ChanNames = [ChanNames(:); MissingElectrodes(:)];


if ~exist('SWnocontour','var') || isempty(SWnocontour)
    SWnocontour = 0;
end


SWrescale = 1;
Clim = [-1 1];
if ~exist('ManualClim','var') || isempty(ManualClim) || numel(ManualClim) ~= 2
    % If ManualClim is not specified or is empty: Full autoscale
    SWrescale = 1;
    Clim = [-1 1];
elseif numel(ManualClim) == 2 && ManualClim(2)-ManualClim(1) > 0
    % If ManualClim is specified as [Cmin, Cmax], scale to this
    SWrescale = 0;
    Clim = [ManualClim(1) ManualClim(2)];
elseif numel(ManualClim) == 2 && ManualClim(2)-ManualClim(1) == 0
    % If ManualClim is specified but is [0 0], no scaling at all
    SWrescale = 0;
    Clim = [NaN NaN];
end


if ~exist('SWnocolorbar','var') || isempty(SWnocolorbar)
    SWnocolorbar = 0;
end

if ~exist('FontSize','var') || isempty(FontSize)
    FontSize = 7;
end

Npts = 512;

if SWrescale
    Filter = Filter ./ max(abs(Filter(:)));
end

% datamax = max(Filter(:));
% datamin = min(Filter(:));
tidxmax = length(time);


% Do not plot if there are more than 1 time point entered.
SWdonotplot = 0;
if tidxmax > 1
    SWdonotplot = 1;
end

if ischar(ChanNames)
    ChanNames = string_to_cell(ChanNames,' ,-_');
end

% % Delete A1 A2 M1 M2
% [Filter,ChanNames] = eeg_deletechans(Filter,ChanNames,'A1 A2 M1 M2');


Nchan = length(ChanNames);
if Nchan ~= size(Filter,1)
    error('data length must be the same as number of channels.');
end

% Convert to projection view (delete channels if they cannot be interpreted)
CoordLookupTable = eeg_projtransform();
SqrGridCoord = eeg_electrode_to_gridcoord(ChanNames);
o = eeg_findcoord(CoordLookupTable.sqr,SqrGridCoord);
nanlist = find(isnan(o));
if ~isempty(nanlist)
    disp(['Cannot convert ' cell_to_string(ChanNames(nanlist),', ') ' into 10-10 international standard coordinate']);
    % Delete invalid channels
    DelChanNames = ChanNames(nanlist);
    [Filter,ChanNames] = eeg_deletechans(Filter,ChanNames,DelChanNames);
    SqrGridCoord = eeg_electrode_to_gridcoord(ChanNames);
    o = eeg_findcoord(CoordLookupTable.sqr,SqrGridCoord);
end
ProjGridCoord = CoordLookupTable.proj(o,:);


% Preload interpolation variables
ti = linspace(-1,1,Npts);
[xi,yi] = meshgrid(ti,ti);
zi = zeros(Npts,Npts,tidxmax);

% process each time point
for tidx = 1:tidxmax
    data = Filter(:,time(tidx));
    
    % Interpolate
    zi(:,:,tidx) = griddata(ProjGridCoord(:,1),ProjGridCoord(:,2),data,xi,yi,'cubic'); %#ok<FPARK>
    
    if SWdonotplot
    else
        % Plot
        fhand = figure;
        hold on
        set(gca,'DataAspectRatio',[1 1 1]);
        if ~isnan(Clim(1))
            caxis(Clim);
        end
        
        surf(xi,yi,zi,'EdgeColor','none');
        
        [temp ch] = contour(xi,yi,zi); %#ok<ASGLU>
        if SWnocontour
            set(ch, 'Visible', 'off');
        end
        cch = get(ch,'Children');
        if isempty(cch)
            cch = ch;
        end
        for i = 1:length(cch)
            xdata = get(cch(i),'xdata');
            ydata = get(cch(i),'ydata');
            MX = max(data);
            p3h = plot3(xdata,ydata,MX*ones(1,length(xdata))+1,'k');
            if SWnocontour
                set(p3h, 'Visible', 'off');
            end
        end
        
        if ~SWnocolorbar
            colorbar
        end
        axis equal
        axis square
        axis off
        
        % plotting nose
        plot3(0, 0.7, MX+2, 'k^', 'MarkerSize', 15);

        for i = 1:Nchan
            if ismember(upper(ChanNames{i}), upper(OnlyShowElectrodes))
                text(ProjGridCoord(i,1),ProjGridCoord(i,2),MX+2,ChanNames{i},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontWeight','bold','FontSize',FontSize,'FontName','VariableWidth');
            end
        end
        
        % plot head shape
%         Hspan = [min(ProjGridCoord(:,1)) max(ProjGridCoord(:,1))];
%         Vspan = [min(ProjGridCoord(:,2)) max(ProjGridCoord(:,2))];
%         Hspan = [min(ProjGridCoord(:,1)) max(ProjGridCoord(:,1))];
%         major = Vspan(2);
%         minor = Hspan(2);
%         t=linspace(0,2*pi,100);
%         x = major*cos(t-pi/2);
%         y = minor*sin(t-pi/2);
%         [th,r]=cart2pol(x,y);
%         [x,y]=pol2cart(th,r);
%         plot(x,y);
        
        hold off
    end
    
end
set(gca,'Xlim',[-0.7 0.7],'YLim',[-0.72 0.72]);


%disp(['Done!']);
