% New version of EEG topographic plot using the 10-10 or 10-5 international
% standard placement system
%

function eeg_topoplot2 (values, ChanNames, Opts)
SWnocontour = 1;
SWnocolorbar = 0;
OnlyShowElectrodes = ChanNames;
FontSize = 12;
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

Npts = 512;
txi = linspace(min(ProjGridCoord(:,1)),max(ProjGridCoord(:,1)),Npts);
tyi = linspace(min(ProjGridCoord(:,2)),max(ProjGridCoord(:,2)),Npts);
[xi,yi] = meshgrid(txi,tyi);

zi = griddata(ProjGridCoord(:,1),ProjGridCoord(:,2),data,xi,yi,'cubic');

fhand = figure;
hold on
set(gca,'DataAspectRatio',[1 1 1]);
% if ~isnan(Clim(1))
%     caxis(Clim);
% end

radius = -min(min(FlatXYCoords));
drawcirc(gca, [0 0], radius);

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
%plot3(0, radius+0.03, MX+2, 'k^', 'MarkerSize', 15);
nosesize = 0.3;
drawnose(gca, radius, nosesize);

for i = 1:length(ElecNames)
    if ismember(upper(ElecNames{i}), upper(OnlyShowElectrodes))
        text(FlatXYCoords(i,1),FlatXYCoords(i,2),MX+2,ElecNames{i},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontWeight','bold','FontSize',FontSize,'FontName','VariableWidth');
    end
end

set(gca,'Xlim',(radius+nosesize*sqrt(2)/2)*[-1 1],'YLim',(radius+nosesize*sqrt(2)/2)*[-1 1]);

if exist('Opts', 'var') && isstruct(Opts)
    if isfield(Opts, 'symmetricclim') && Opts.symmetricclim
        CLim = max(abs(data(:)))*[-1 1];
        set(gca, 'CLim', CLim);
    end
end


function drawcirc (ax, origin, radius)
[x,y] = pol2cart(linspace(-pi,pi,3600), radius);
plot(ax, origin(1)+x, origin(2)+y, 'k');

function drawnose(ax, radius, nosesize)
plot(ax, [0, -nosesize/2], [radius+nosesize*sqrt(2)/2, radius], 'k');
plot(ax, [0, nosesize/2], [radius+nosesize*sqrt(2)/2, radius], 'k');
