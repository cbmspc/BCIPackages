function Signal = notchfilter (Signal, Fs, PowerLineFrequency, Order, QFactor)
if ~exist('Order','var') || isempty(Order)
    Order = 4;
end
if ~exist('QFactor','var') || isempty(QFactor)
    QFactor = 10;
end

Funda = PowerLineFrequency;

fprintf('Designing notch filters.\n');
for h = 1:floor(Fs/2/Funda)
    % Create one for each harmonic
    d = fdesign.notch('N,F0,Q',Order,Funda*h/(Fs/2),QFactor*h);
    Hd{h} = design(d);
end


for h = 1:length(Hd)
    fprintf('Filtering at harmonic# %i of %i.\n', h, length(Hd));
    Signal = filtfilt(Hd{h}.sosMatrix,Hd{h}.ScaleValues,Signal);
end

