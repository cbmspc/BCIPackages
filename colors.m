function colors (Number_of_colors)
color_dir = strrep(cell_to_string({matlabroot 'toolbox' 'matlab' 'graphics' 'color'}, '|'), '|', filesep);
if exist(color_dir, 'dir')
    dirlist = dir([color_dir filesep '*.m']);
    flist = regexprep({dirlist.name}, '\.m$', '');
    ivc = false(1,length(flist));
    for f = 1:length(flist)
        if isequal(flist{f}, 'colororder')
            continue
        end
        try %#ok<TRYNC>
            clist = feval(flist{f},63);
            if size(uniquerows(clist),1) > 1
                ivc(f) = true;
            end
        end
    end
    flist = fliplr(flist(ivc));
else
    flist = {'jet' 'hsv' 'hot' 'pink' 'flag' 'gray' 'cool' 'copper' 'bone'};
end
F = length(flist);
fprintf('Creating color bars for %i different maps...\n', length(flist));
fighand = figure('Name', 'Colors!', 'ToolBar', 'none', 'MenuBar', 'none', 'NumberTitle', 'off', 'Color', [0 0 0]);

N = 63;
if exist('Number_of_colors','var') && ~isempty(Number_of_colors) && numel(Number_of_colors) == 1 && isnumeric(Number_of_colors)
    Number_of_colors = round(abs(Number_of_colors));
    if Number_of_colors > 0
        N = Number_of_colors;
    end
else
    fprintf('Hint: The default is to draw %i colors per map. You can change this by passing in a number.\n', N);

end

set(fighand,'NextPlot','add');
for f = 1:F
    if ~ishandle(fighand)
       break;
    end
    Kolor = feval(flist{f},N);
    axehand = axes(fighand, 'Position', [0 (f-1)/F 1 1/F], 'NextPlot', 'add', 'Color', [0 0 0], 'TickLength', [0 0], 'Box', 'on', 'XTickLabel', [], 'YTickLabel', [], 'YLim', [0 1], 'XLim', [0 1]);
    for i = 1:N
        patch(axehand, ([-1 0 0 -1]+i)/N,[0 0 1 1],Kolor(i,:),'EdgeColor','none');
    end
    text(axehand, 0.5, 0.5, [' ' flist{f} ' '], 'Color', [0 0 0], 'BackgroundColor', [0.7 0.7 0.7], 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Margin', 0.001, 'FontName', 'FixedWidth');
end
