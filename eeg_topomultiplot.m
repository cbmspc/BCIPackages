function eeg_topomultiplot (Filter, ChanNames, ~, TimePointsToPlot, FilePrefix, TimePointsLabels, FullHeadElectrodes, OnlyShowElectrodes, SWnocontour, SWtitle, ManualClim, SWnocolorbar, FontSize)
close all
if ~exist('SWtitle','var') || isempty(SWtitle)
    SWtitle = 1;
end

if ~exist('ManualClim','var')
    ManualClim = [];
end
if ~exist('SWnocolorbar','var') 
    SWnocolorbar = [];
end
if ~exist('FontSize','var')
    FontSize = [];
end
if ~exist('TimePointsLabels','var') || numel(TimePointsLabels) ~= numel(TimePointsToPlot)
    TimePointsLabels = TimePointsToPlot;
end
if ~iscell(TimePointsLabels)
    TimePointsLabels = string_to_cell(num2str(TimePointsLabels), ' ');
end
%if ~exist('SWfigfile','var') || isempty(SWfigfile)
SWfigfile = 1;
%end
for x = 1:length(TimePointsToPlot)
    xx = TimePointsToPlot(x);
    [tmp, tmp, tmp, fhand] = eeg_topoplot(Filter,ChanNames,xx,FullHeadElectrodes,OnlyShowElectrodes,SWnocontour,ManualClim,SWnocolorbar,FontSize); %#ok<ASGLU>
    if SWtitle
        %title(['Time: ' num2str(TimeOffset+xx/Fs*1000,'%.0f') ' msec'], 'Interpreter', 'none');
        title(['Time/freq: ' TimePointsLabels{x}], 'Interpreter', 'none');
    end
    scrsz = get(0,'ScreenSize');
    set(gcf,'Position',scrsz);
    pause(0.01);
    if ~isempty(FilePrefix)
        FileName = [FilePrefix ' point ' num2str(xx,'%03i')];
        if SWfigfile || ~ispc
            saveas(fhand,FileName,'fig');
        else
            print(fhand,'-dpng',FileName);
        end
        pause(0.01);
        close all
    end
end
