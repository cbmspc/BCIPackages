% New version of EEG topographic plot using the 10-10 or 10-5 international
% standard placement system
%

function [SurfHand, FlatXYCoords_active, ElecNames_active, SurfData, ContourHand] = eeg_topoplot2 (values, ChanNames, Opts)
SWnocontour = 1;
SWnosurf = 0;
SWnocolorbar = 0;
NumContourLevels = 10;
OnlyShowElectrodes = ChanNames;
FontSize = 11;
FontName = 'VariableWidth';
Npts = 512;
axes_hand = [];

if length(ChanNames) > 128
    FontSize = 6;
end

if exist('Opts', 'var') && isstruct(Opts)
    if isfield(Opts, 'contour')
        if Opts.contour
            SWnocontour = 0;
        else
            SWnocontour = 1;
        end
    end
    if isfield(Opts, 'nocontour')
        if Opts.nocontour
            SWnocontour = 1;
        else
            SWnocontour = 0;
        end
    end
    if isfield(Opts, 'surf')
        if Opts.surf
            SWnosurf = 0;
        else
            SWnosurf = 1;
        end
    end
    if isfield(Opts, 'nosurf')
        if Opts.nosurf
            SWnosurf = 1;
        else
            SWnosurf = 0;
        end
    end
    if isfield(Opts, 'nocolorbar')
        if Opts.nocolorbar
            SWnocolorbar = 1;
        else
            SWnocolorbar = 0;
        end
    end
    if isfield(Opts, 'numcontourlevels') && ~isempty(Opts.numcontourlevels)
        NumContourLevels = Opts.numcontourlevels;
    end
    if isfield(Opts, 'fontsize') && ~isempty(Opts.fontsize)
        FontSize = Opts.fontsize;
    end
    if isfield(Opts, 'fontname') && ~isempty(Opts.fontname)
        FontName = Opts.fontname;
    end
    if isfield(Opts, 'npts') && ~isempty(Opts.npts)
        Npts = Opts.npts;
    end
    if isfield(Opts, 'axes_hand') && ~isempty(Opts.axes_hand) && ishandle(Opts.axes_hand)
        axes_hand = Opts.axes_hand;
        axes(axes_hand);
    else
        axes_hand = [];
    end
end


[~, ElecNames, FlatXYCoords] = eeg_electrode_to_gridcoord_3d ('all');
ElecNames = ElecNames(1:336);
FlatXYCoords = FlatXYCoords(1:336,:);

ChanNames = regexprep(ChanNames, '^FCC([78])', 'FTT$1');
ChanNames = regexprep(ChanNames, '^CCP([78])', 'TTP$1');
ChanNames = regexprep(ChanNames, '^FC([78])', 'FT$1');
ChanNames = regexprep(ChanNames, '^CP([78])', 'TP$1');
if isempty(chan2idx(ChanNames,'TP9',1))
    ChanNames = regexprep(ChanNames, '^M1$', 'TP9');
end
if isempty(chan2idx(ChanNames,'TP10',1))
    ChanNames = regexprep(ChanNames, '^M2$', 'TP10');
end

cid = chan2idx(ElecNames, ChanNames);
if any(cid==0)
    % Delete any unrecognized channels and warn.
    ChanNames = ChanNames(cid>0);
    values = values(cid>0);
    cid = cid(cid>0);
    warning('EEG electrodes not recognized: %s', cell_to_string(ChanNames(cid==0), ' '));
end
data = nan(length(ElecNames),1);
data(cid,:) = values;
ProjGridCoord = FlatXYCoords;

ProjGridCoord = ProjGridCoord(~isnan(data),:);
data = data(~isnan(data),:);

txi = linspace(min(ProjGridCoord(:,1)),max(ProjGridCoord(:,1)),Npts);
tyi = linspace(min(ProjGridCoord(:,2)),max(ProjGridCoord(:,2)),Npts);
[xi,yi] = meshgrid(txi,tyi);

zi = griddata(ProjGridCoord(:,1),ProjGridCoord(:,2),data,xi,yi,'cubic');

% fhand = figure;
hold on
set(gca,'DataAspectRatio',[1 1 1]);
% if ~isnan(Clim(1))
%     caxis(Clim);
% end

radius = -min(min(FlatXYCoords));
drawcirc(gca, [0 0], radius, 2);

if ~isempty(axes_hand)
    SurfHand = surf(axes_hand, xi, yi, zi, 'EdgeColor', 'none');
    [~, ContourHand] = contour(axes_hand, xi, yi, zi, NumContourLevels, 'LineWidth', 2);
else
    SurfHand = surf(xi, yi, zi, 'EdgeColor', 'none');
    [~, ContourHand] = contour(xi, yi, zi, NumContourLevels, 'LineWidth', 2);
end

SurfData.xi = xi;
SurfData.yi = yi;
SurfData.zi = zi;

if SWnocontour
    set(ContourHand, 'Visible', 'off');
end
if SWnosurf
    set(SurfHand, 'Visible', 'off');
end

cch = get(ContourHand,'Children');
if isempty(cch)
    cch = ContourHand;
end
for i = 1:length(cch)
    xdata = get(cch(i),'xdata');
    ydata = get(cch(i),'ydata');
    MX = max(data);
    if ~SWnosurf
        % Draw the contour lines as black if surface plot is enabled
        p3h = plot3(xdata,ydata,MX*ones(1,length(xdata))+1,'k');
        if SWnocontour
            set(p3h, 'Visible', 'off');
        end
    end
end

if ~SWnocolorbar
    colorbar('FontSize', 16);
end
axis equal
axis square
axis off

% plotting nose
%plot3(0, radius+0.03, MX+2, 'k^', 'MarkerSize', 15);
nosesize = 0.3;
drawnose(gca, radius, nosesize, 2);

for i = 1:length(ElecNames)
    if ismember(upper(ElecNames{i}), upper(OnlyShowElectrodes))
        text(FlatXYCoords(i,1),FlatXYCoords(i,2),MX+2,ElecNames{i},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontWeight','bold','FontSize',FontSize,'FontName',FontName);
    end
end

set(gca,'Xlim',(radius+nosesize*sqrt(2)/2)*[-1 1],'YLim',(radius+nosesize*sqrt(2)/2)*[-1 1]);

if exist('Opts', 'var') && isstruct(Opts)
    if isfield(Opts, 'symmetricclim') && Opts.symmetricclim
        CLim = max(abs(data(:)))*[-1 1];
        set(gca, 'CLim', CLim);
    end
    if isfield(Opts, 'clim') && length(Opts.clim) == 2
        set(gca, 'CLim', Opts.clim);
    end
end

tmp = chan2lidx(upper(ElecNames),upper(OnlyShowElectrodes));
FlatXYCoords_active = FlatXYCoords(tmp,:);
ElecNames_active = ElecNames(tmp);


function drawcirc (ax, origin, radius, linewidth)
[x,y] = pol2cart(linspace(-pi,pi,3600), radius);
plot(ax, origin(1)+x, origin(2)+y, 'k', 'LineWidth', linewidth);

function drawnose(ax, radius, nosesize, linewidth)
plot(ax, [0, -nosesize/2], [radius+nosesize*sqrt(2)/2, radius], 'k', 'LineWidth', linewidth);
plot(ax, [0, nosesize/2], [radius+nosesize*sqrt(2)/2, radius], 'k', 'LineWidth', linewidth);
