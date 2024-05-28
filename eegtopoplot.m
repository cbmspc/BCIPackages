function eegtopoplot (data, ChanNames)

error('Outdated function. Should use eeg_topoplot2 instead');


Npts = 64;

data = data(:);
if isstr(ChanNames)
    ChanNames = string_to_cell(ChanNames,' ,-_');
end
Nchan = length(ChanNames);
if Nchan ~= length(data)
    error('data length must be the same as number of channels.');
end

DataGrid = zeros(13,11);

X = -5:5;
Y = -5:5;

for i = 1:Nchan
    [V,H] = eeg_electrode_to_gridcoord(ChanNames{i});
    GridCoord(i,:) = [V H];
    y = V+6;
    x = H+7;
    DataGrid(x,y) = data(i);
end

EEGGrid = DataGrid(2:end-1,:).';

ti = linspace(-5,5,Npts);
[xi,yi] = meshgrid(ti,ti);
zi = griddata(X,Y,EEGGrid,xi,yi,'cubic');

clf

for i = 1:Nchan
    text(GridCoord(i,2),GridCoord(i,1),ChanNames{i},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[1 1 1],'Color',[0 0 0],'FontWeight','bold','FontSize',12,'FontName','FixedWidth');
end

hold on

contour(xi,yi,zi);

hold off
