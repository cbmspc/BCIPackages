function eeg_topomultiplot2 (values, ChanNames, SubplotNames)

%values = abs(values);
%vflat = abs(values(:));
vflat = values(:);

%Opts.clim = [0 1]*quantile(vflat, .95);
Opts.clim = [-1 1]*max(abs(vflat));
%Opts.symmetricclim = 1;
Opts.fontsize = 4;
Opts.nocolorbar = 1;
Opts.nocontour = 1;
Opts.npts = 64;
[hands, underhand] = subplotcompact(floor(sqrt(size(values,2))), ceil(sqrt(size(values,2))));



for s = 1:size(values,2)
    Opts.axes_hand = hands(s);
    eeg_topoplot2(values(:,s), ChanNames, Opts);
    %title(SubplotNames{s});
    title([SubplotNames{s} '       '], 'HorizontalAlignment','right','VerticalAlignment','cap')
    pause(0.05);
end
