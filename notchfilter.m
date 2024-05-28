function Signal = notchfilter (Signal, Fs, PowerLineFrequency, Order, QFactor)
if ~exist('Order','var') || isempty(Order)
    Order = 4;
end
if ~exist('QFactor','var') || isempty(QFactor)
    QFactor = 10;
end

Funda = PowerLineFrequency;
persistent Hd_prev Funda_prev Order_prev QFactor_prev

if iscell(Hd_prev) && Funda_prev == Funda && Order_prev == Order && QFactor_prev == QFactor
    Hd = Hd_prev;
else
    %fprintf('Designing notch filters.\n');
    for h = floor(Fs/2/Funda):-1:1
        % Create one for each harmonic
        d = fdesign.notch('N,F0,Q',Order,Funda*h/(Fs/2),QFactor*h);
        Hd{h} = design(d);
    end
    Funda_prev = Funda;
    Order_prev = Order;
    QFactor_prev = QFactor;
    Hd_prev = Hd;
end

for h = 1:length(Hd)
    %fprintf('Filtering at harmonic# %i of %i.\n', h, length(Hd));
    Signal = filtfilt(Hd{h}.sosMatrix,Hd{h}.ScaleValues,double(Signal));
end

