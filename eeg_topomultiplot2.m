function eeg_topomultiplot2 (values, ChanNames, SubplotNames, CLim, GridDim, headcolors, OtherOpts)
if exist('CLim', 'var') && length(CLim) == 2
    Opts.clim = CLim;
else
    Opts.clim = [-1 1]*max(abs(values(:)));
end

Opts.fontsize = 4;
Opts.normalizedfontsize = 0.03;
Opts.nocolorbar = 1;
Opts.nocontour = 1;
Opts.npts = 64;

if exist('GridDim', 'var') && length(GridDim) == 2 && size(values,2) == GridDim(1) * GridDim(2)
    [hands, ~] = subplotcompact(GridDim(1), GridDim(2));
else
    [hands, ~] = subplotcompact(floor(sqrt(size(values,2))), ceil(sqrt(size(values,2))));
end

if exist('OtherOpts', 'var') && isstruct(OtherOpts)
    fns = fieldnames(OtherOpts);
    for i = 1:length(fns)
        if ~isfield(Opts, fns{i}) && ~strcmp(fns{i}, 'title') && ~strcmp(fns{i}, 'footer')
            Opts.(fns{i}) = OtherOpts.(fns{i});
        elseif strcmp(fns{i}, 'normalizedfontsize') || strcmp(fns{i}, 'fontsize')
            Opts.(fns{i}) = OtherOpts.(fns{i});
        end
    end
end

for s = 1:size(values,2)
    Opts.axes_hand = hands(s);
    if exist('headcolors', 'var')
        Opts.headcolor = headcolors(mod(s-1,size(headcolors,1))+1,:,:);
    end
    if exist('OtherOpts', 'var') && isstruct(OtherOpts) && isfield(OtherOpts, 'bordercolorlocations') && size(OtherOpts.bordercolorlocations,1) == size(values,2) && size(OtherOpts.bordercolorlocations,2) == 7
        Opts.bordercolorlocation = OtherOpts.bordercolorlocations(mod(s-1,size(OtherOpts.bordercolorlocations,1))+1,:,:);
    end
    if exist('OtherOpts', 'var') && isstruct(OtherOpts) && isfield(OtherOpts, 'SuppressElectrodeLabelsPerSubplot') && numel(OtherOpts.SuppressElectrodeLabelsPerSubplot) == size(values,2) && iscell(OtherOpts.SuppressElectrodeLabelsPerSubplot{1})
        Opts.SuppressElectrodeLabels = OtherOpts.SuppressElectrodeLabelsPerSubplot{s};
    end

    % Only pass down channels that aren't infinity or nan
    thisvalue = values(:,s);
    thischannames = ChanNames;
    thisisfinite = isfinite(thisvalue);
    thisvalue = thisvalue(thisisfinite);
    thischannames = thischannames(thisisfinite);
    eeg_topoplot2(thisvalue, thischannames, Opts);
    %title(SubplotNames{s});
    title([SubplotNames{s} '       '], 'HorizontalAlignment','right','VerticalAlignment','cap','FontUnit','normalized','FontSize',0.06)
    pause(0.05);
end

try %#ok<TRYNC>
    uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.10 0.955 0.80 0.05], 'String', OtherOpts.title, 'FontUnits', 'normalized', 'FontSize', 0.7, 'Visible', 'on', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.05 0.015 0.90 0.035], 'String', OtherOpts.footer, 'FontUnits', 'normalized', 'FontSize', 0.7, 'Visible', 'on', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    %sgtitle(OtherOpts.title);
end

