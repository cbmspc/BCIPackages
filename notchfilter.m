function Signal = notchfilter (Signal, Fs, NotchFreq, Order, QFactor, SplitNotch)
if ~exist('Order','var') || isempty(Order)
    Order = 4;
end
if ~exist('QFactor','var') || isempty(QFactor)
    QFactor = 10;
end
if ~exist('SplitNotch','var') || numel(SplitNotch) ~= 1
    SplitNotch = false;
end

persistent Hd_prev Funda_prev Order_prev QFactor_prev
NotchFreq = NotchFreq(:);
NotchFreq(~(isfinite(NotchFreq) & NotchFreq > 0 & mod(NotchFreq,Fs/2) > 0)) = NaN;
NotchFreq = unique(NotchFreq(~isnan(NotchFreq)));
PerceFreq = abs(NotchFreq - Fs .* round(NotchFreq ./ Fs));
if isempty(NotchFreq)
    return
end
if iscell(Hd_prev) && isequal(Funda_prev,NotchFreq) && Order_prev == Order && QFactor_prev == QFactor
    Hd = Hd_prev;
else
    if length(NotchFreq) == 1
        % Design notch filters for NotchFreq and all of its harmonics
        for h = floor(Fs/2/NotchFreq):-1:1
            % Create one for each harmonic
            d = fdesign.notch('N,F0,Q',Order,NotchFreq*h/(Fs/2),QFactor*h);
            Hd{h} = design(d);
        end
        Funda_prev = NotchFreq;
        Order_prev = Order;
        QFactor_prev = QFactor;
        Hd_prev = Hd;
    else
        % Design notch filters for the frequencies listed in NotchFreq
        for h = length(NotchFreq):-1:1
            % Create one for each frequency in the list
            d = fdesign.notch('N,F0,Q',Order,PerceFreq(h)/(Fs/2),QFactor/NotchFreq(1)*NotchFreq(h));
            Hd{h} = design(d);
        end
        Funda_prev = NotchFreq;
        Order_prev = Order;
        QFactor_prev = QFactor;
        Hd_prev = Hd;
    end
end

Signal = double(Signal);

if SplitNotch
    Signalfu = flipud(Signal);
    for h = 1:length(Hd)
        Signal = filter(Hd{h},Signal);
        Signalfu = filter(Hd{h},Signalfu);
    end
    Signalfu = flipud(Signalfu);
    Ntp = size(Signal,1);
    Nch = size(Signal,2);
    factors = 1 - 1./(1+exp((1:Ntp).'-Ntp/2));
    for ch = 1:Nch
        Signal(:,ch) = Signal(:,ch) .* factors + Signalfu(:,ch) .* (1-factors);
    end
else
    for h = 1:length(Hd)
        Signal = filtfilt(Hd{h}.sosMatrix,Hd{h}.ScaleValues,Signal);
    end
end

