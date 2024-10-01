function eeg_topomultiplot2 (values, ChanNames, SubplotNames, CLim, GridDim, headcolors, OtherOpts)
if exist('CLim', 'var') && length(CLim) == 2
    Opts.clim = CLim;
else
    Opts.clim = [-1 1]*max(abs(values(:)));
end

Opts.fontsize = 4;
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
        end
    end
end

for s = 1:size(values,2)
    Opts.axes_hand = hands(s);
    if exist('headcolors', 'var')
        Opts.headcolor = headcolors(mod(s-1,size(headcolors,1))+1,:,:);
    end

    % Only pass down channels that aren't infinity or nan
    thisvalue = values(:,s);
    thischannames = ChanNames;
    thisisfinite = isfinite(thisvalue);
    thisvalue = thisvalue(thisisfinite);
    thischannames = thischannames(thisisfinite);

    eeg_topoplot2(thisvalue, thischannames, Opts);
    %title(SubplotNames{s});
    title([SubplotNames{s} '       '], 'HorizontalAlignment','right','VerticalAlignment','cap')
    pause(0.05);
end

try %#ok<TRYNC>
    uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.10 0.96 0.80 0.04], 'String', OtherOpts.title, 'FontUnits', 'normalized', 'FontSize', 0.6, 'Visible', 'on', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    uicontrol(gcf, 'Style', 'text', 'Units', 'normalized', 'Position', [0.05 0.02 0.90 0.02], 'String', OtherOpts.footer, 'FontUnits', 'normalized', 'FontSize', 0.6, 'Visible', 'on', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    %sgtitle(OtherOpts.title);
end

