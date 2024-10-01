function Signal = notchfilter (Signal, Fs, NotchFreq, Order, QFactor)
if ~exist('Order','var') || isempty(Order)
    Order = 4;
end
if ~exist('QFactor','var') || isempty(QFactor)
    QFactor = 10;
end

persistent Hd_prev Funda_prev Order_prev QFactor_prev

if iscell(Hd_prev) && isequal(Funda_prev,NotchFreq) && Order_prev == Order && QFactor_prev == QFactor
    Hd = Hd_prev;
else
    NotchFreq = NotchFreq(NotchFreq < Fs/2 & NotchFreq > 0);
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
            d = fdesign.notch('N,F0,Q',Order,NotchFreq(h)/(Fs/2),QFactor*h);
            Hd{h} = design(d);
        end
    end
end

for h = 1:length(Hd)
    Signal = filtfilt(Hd{h}.sosMatrix,Hd{h}.ScaleValues,double(Signal));
end

