function colors ()
flist = {'jet' 'hsv' 'hot' 'pink' 'flag' 'gray' 'cool' 'copper' 'bone'};
F = length(flist);
figure;
clf;
set(gcf, 'Name', 'Colors!');
set(gcf, 'ToolBar', 'none', 'MenuBar', 'none', 'NumberTitle', 'off');
set(gcf, 'Color', [0 0 0]);

for f = 1:F
    N = 255;
    Kolor = feval(flist{f},N);
    
    %Kolor = feval(flist{f},N*2);
    %Kolor = Kolor(N:-1:1,:);
    
    axes('Position', [0 (f-1)/F 1 1/F]);
    %axis equal
    set(gca, 'Color', [0 0 0]);
    set(gca, 'TickLength', [0 0]);
    set(gca, 'Box', 'on');
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    set(gca, 'YLim', [0 1], 'XLim', [0 1]);
    hold on
    for i = 1:N
        %text(floor((i-1)/10)+1,mod(i-1,10)+1, num2str(i,'%03i'),'Color',Kolor(i,:), 'HorizontalAlignment', 'center');
        patch(([-0.5 0.5 0.5 -0.5]+i)/N,[0 0 1 1],Kolor(i,:),'EdgeColor','none');
    end
    text(0.5, 0.5, [' ' flist{f} ' '], 'Color', [1 1 1], 'BackgroundColor', [0 0 0], 'HorizontalAlignment', 'center', 'FontSize', 12);
end
