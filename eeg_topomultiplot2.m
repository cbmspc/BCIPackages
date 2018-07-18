function eeg_topomultiplot2 (values, ChanNames, SubplotNames)

values = abs(values);
vflat = abs(values(:));

Opts.clim = [0 1]*quantile(vflat, .95);
Opts.fontsize = 4;
Opts.nocolorbar = 1;
Opts.nocontour = 1;
Opts.npts = 64;
[hands, underhand] = subplotcompact(floor(sqrt(size(values,2))), ceil(sqrt(size(values,2))));



for s = 1:size(values,2)
    Opts.axes_hand = hands(s);
    eeg_topoplot2(values(:,s), ChanNames, Opts);
    pause(0.05);
end
