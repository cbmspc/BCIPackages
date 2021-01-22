% Cluster and identify electrodes from 2D powerpoint images
% Electrodes must be marked by filled red circles in powerpoint first

function [centroid, img] = clusterelectrodes (image_file, N, MainRGBList, RGBTolList)
img = imread(image_file);

if nargin < 5
    MainRGBList = [0 112 192; 68 114 196; 255 0 0; 192 0 0];
    RGBTolList = [5 5 5; 5 5 5; 3 3 3; 3 3 3];
end

for i = size(MainRGBList,1):-1:1
    idot(:,:,i) = img(:,:,1)>=MainRGBList(i,1)-RGBTolList(i,1) & img(:,:,1)<=MainRGBList(i,1)+RGBTolList(i,1) & ...
        img(:,:,2)>=MainRGBList(i,2)-RGBTolList(i,2) & img(:,:,2)<=MainRGBList(i,2)+RGBTolList(i,2) & ...
        img(:,:,3)>=MainRGBList(i,3)-RGBTolList(i,3) & img(:,:,3)<=MainRGBList(i,3)+RGBTolList(i,3);
end

idot = logical(sum(idot,3));

xylist = zeros(nnz(idot),2);
k = 0;
for i = 1:size(idot,1)
    for j = 1:size(idot,2)
        if idot(i,j)
            k = k + 1;
            xylist(k,:) = [i j];
        end
    end
end


% Random subsampling when list is too big
if size(xylist,1) > 5000
    s = randperm(size(xylist,1));
    xylist = xylist(s(1:5000),:);
end

close all
image(img);
try
    pause(0.00001);
    oldWarningState = warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');
    frame_h = get(handle(1),'JavaFrame');
    set(frame_h,'Maximized',1);
    warning(oldWarningState);
    pause(0.00001);
end

hold on
set(gca, 'DataAspectRatio', [1 1 1]);
pause(0.1);

%ids = dbscan(xylist,9);
ids = ecogfinger_clustremix (xylist, N);

classes = unique(ids);
% 
% img2 = img;
centroid = zeros(length(classes),2);
for c = 1:length(classes)
    subxy = xylist(ids == classes(c),:);
    centroid(c,:) = median(subxy);
end

centroid = sortrows(round(centroid),-1);



for i = 1:size(centroid,1)
    ph1(i) = plot(centroid(i,2), centroid(i,1), 'ok', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'b');
    th1(i) = text(centroid(i,2), centroid(i,1), lower(dec2base(i,36)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontName', 'Arial Narrow', 'FontSize', 11, 'Color', [1 1 0]);
end

set(gcf, 'Position', [1           1        1600        960]);


uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.010 0.960 0.055 0.022], 'BackgroundColor', [1 1 1], 'String', ['Subject ID = '], 'HorizontalAlignment', 'left', 'FontSize', 11);
h_subject = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.070 0.960 0.055 0.022], 'BackgroundColor', [1 1 .8], 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 12);

uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.130 0.960 0.035 0.022], 'BackgroundColor', [1 1 1], 'String', ['ImgFile = '], 'HorizontalAlignment', 'left', 'FontSize', 10);
h_infile = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.170 0.960 0.295 0.022], 'BackgroundColor', [1 1 .8], 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 10);

uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.130 0.935 0.035 0.022], 'BackgroundColor', [1 1 1], 'String', ['OutFile = '], 'HorizontalAlignment', 'left', 'FontSize', 10);
h_outfile = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.170 0.935 0.295 0.022], 'BackgroundColor', [1 1 .8], 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 10);


uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.010 0.800 0.200 0.09], 'BackgroundColor', [1 1 1], 'String', ['Enter the grid/strip name and the IDs (shown in image) in ascending order separated by space. One grid/strip at a time.' sprintf('\n') 'Example: LTG 7 6 5 4 b c d e' sprintf('\n') 'Then click Confirm'], 'HorizontalAlignment', 'left', 'FontSize', 11);
h_grid = uicontrol(gcf, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.010 0.830-0.180 0.200 0.066*2], 'BackgroundColor', [.8 1 .8], 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Max', 2, 'Min', 0);
h_confirm = uicontrol(gcf, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.010 0.610 0.050 0.03], 'BackgroundColor', [0.7 0.7 0.7], 'String', 'Confirm');


uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.010 0.560 0.200 0.022], 'BackgroundColor', [1 1 1], 'String', 'ElectrodeNameCoordinates = ', 'HorizontalAlignment', 'left', 'FontSize', 11);
h_content1 = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.010 0.05 0.090 0.5], 'BackgroundColor', [.7 .7 1], 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 8, 'FontName', 'Consolas', 'Max', 2, 'Min', 0);
h_content2 = uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.110 0.05 0.090 0.5], 'BackgroundColor', [.7 .7 1], 'String', '', 'HorizontalAlignment', 'left', 'FontSize', 8, 'FontName', 'Consolas', 'Max', 2, 'Min', 0);


set(h_confirm, 'Callback', @confirm);

fSubject = regexp(image_file, 'SJ\d{4,}', 'match', 'once');
if ~isempty(fSubject)
    set(h_subject, 'String', fSubject);
end

set(h_infile, 'String', image_file);
out_file = regexprep(image_file, '\....$', '.mat');
set(h_outfile, 'String', out_file);

ElectrodeNameCoordinates = {};

    function confirm (~, ~)
        Subject = get(h_subject, 'String');
        OutFile = get(h_outfile, 'String');
        Entry = get(h_grid, 'String');
        
        if isempty(Subject)
            uiwait(msgbox('Missing Subject ID field.', 'cluster', 'error', 'modal'));
            return
        end
        
        if isempty(OutFile)
            uiwait(msgbox('Missing OutFile field.', 'cluster', 'error', 'modal'));
            return
        end

        if isempty(Entry)
            uiwait(msgbox('Missing Entry field.', 'cluster', 'error', 'modal'));
            return
        end
        
        EntryC = string_to_cell(Entry, ', ');
        GridName = EntryC{1};
        
        IDs = base2dec(EntryC(2:end), 36);
        
        for i = 1:length(IDs)
            set(th1(IDs(i)), 'String', num2str(i), 'Color', [0 0 0]);
            set(ph1(IDs(i)), 'MarkerFaceColor', [0.3 0.3 0.3]);
            pause(0.1);
        end
        
        
        try
            SubjectElectrodeNameCoordinates = evalin('base', 'SubjectElectrodeNameCoordinates');
        catch
            SubjectElectrodeNameCoordinates = {};
        end
        
        for i = 1:length(IDs)
            fprintf('%s%02i\t%i\t%i\n', GridName, i, centroid(IDs(i),1), centroid(IDs(i),2));
            ElectrodeNameCoordinates = [
                ElectrodeNameCoordinates
                {sprintf('%s%02i', GridName, i), centroid(IDs(i),1), centroid(IDs(i),2)}
                ];
            SubjectElectrodeNameCoordinates = [
                SubjectElectrodeNameCoordinates
                {Subject, sprintf('%s%02i', GridName, i), centroid(IDs(i),1), centroid(IDs(i),2)}
                ];
        end
        ElectrodeNameCoordinates = uniquerows(ElectrodeNameCoordinates);
        SubjectElectrodeNameCoordinates = uniquerows(SubjectElectrodeNameCoordinates);
        
        save(OutFile, 'ElectrodeNameCoordinates', 'img');
        
        ENC_content1 = '';
        ENC_content2 = '';
        for i = 1:ceil(size(ElectrodeNameCoordinates,1)/2)
            ENC_content1 = [ENC_content1 sprintf('%s(%i,%i)\n', ElectrodeNameCoordinates{i,1}, ElectrodeNameCoordinates{i,2}, ElectrodeNameCoordinates{i,3})];
        end
        for i = ceil(size(ElectrodeNameCoordinates,1)/2)+1:size(ElectrodeNameCoordinates,1)
            ENC_content2 = [ENC_content2 sprintf('%s(%i,%i)\n', ElectrodeNameCoordinates{i,1}, ElectrodeNameCoordinates{i,2}, ElectrodeNameCoordinates{i,3})];
        end
        set(h_content1, 'String', ENC_content1);
        set(h_content2, 'String', ENC_content2);
        
        assignin('base', 'SubjectElectrodeNameCoordinates', SubjectElectrodeNameCoordinates);
        set(h_grid, 'String', '');
    end

end


