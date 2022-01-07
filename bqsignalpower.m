% Calculates the spectral power of time segments in a real time series
% signal at particular frequency ranges using causal biquad bandpass filter
% bank

function [signal_power, opts] = bqsignalpower (signal, Fs, FreqRange, TimeRange, opts)
if size(signal,1) < size(signal,2)
    signal = signal.';
end

if ~exist('opts','var') || ~isstruct(opts) || isempty(opts)
    opts = struct;
end

if ~isfield(opts, 'FilterIter')
    opts.FilterIter = 2;
end

if ~isfield(opts, 'ReflectSignal')
    opts.ReflectSignal = 0;
end


SigLen = size(signal,1);
MinFreq = min(FreqRange(:));
opts.ReflectSigLen = round(5/MinFreq*Fs);
if opts.ReflectSigLen > SigLen
    opts.ReflectSigLen = SigLen;
end
Nchan = size(signal,2);
Ntimebin = size(TimeRange,1);
Nfreqbin = size(FreqRange,1);
TimeIndexRange = TimeRange*Fs+1;
TimeIndexRange(:,2) = TimeIndexRange(:,2) - 1;
signal_power = nan(Nchan, Nfreqbin, Ntimebin);
if opts.ReflectSignal
    signal_reflected = [flipud(signal(1:opts.ReflectSigLen,:)); signal];
end

%silent = 1;
a0 = 1;

%wh = waitbar(0, ['bqsignalpower Band 0/' num2str(Nfreqbin)]);



for b = 1:Nfreqbin
    Fc1 = FreqRange(b, 1);
    Fc2 = FreqRange(b, 2);
    [a1, a2, scale] = Calculate4thOrderBandpassFilterCoefficients(Fs, Fc1, Fc2);
    opts.filtcoefs(b,:) = [a1, a2, scale];
    b0 = 1*scale;
    b1 = 0*scale;
    b2 = -1*scale;
    if opts.ReflectSignal
        signal_reflected_filtered = filter([b0 b1 b2], [a0 a1 a2], signal_reflected);
        %signal_reflected_filtered = CustomFilter(a1, a2, scale, signal_reflected);
        signal_filtered = signal_reflected_filtered(opts.ReflectSigLen+1:end,:);
        for i = 2:opts.FilterIter
            signal_filtered = filter([b0 b1 b2], [a0 a1 a2], signal_filtered);
            %signal_filtered = CustomFilter(a1, a2, scale, signal_filtered);
        end
    else
        signal_filtered = filter([b0 b1 b2], [a0 a1 a2], signal);
        %signal_filtered = CustomFilter(a1, a2, scale, signal);
        for i = 2:opts.FilterIter
            signal_filtered = filter([b0 b1 b2], [a0 a1 a2], signal_filtered);
            %signal_filtered = CustomFilter(a1, a2, scale, signal_filtered);
        end
    end
    signal_filtered_squared = signal_filtered.^2;
    
    for k = 1:Ntimebin
        t1 = TimeIndexRange(k,1);
        t2 = TimeIndexRange(k,2);
        signal_power(:,b,k) = mean(signal_filtered_squared(t1:t2,:),1);
    end
    %waitbar(b/Nfreqbin, wh, ['bqsignalpower Band ' num2str(b) '/' num2str(Nfreqbin)]);
    
end
%delete(wh);


function [a1, a2, scale] = Calculate4thOrderBandpassFilterCoefficients(Fs, Fc1, Fc2)
r1 = pi*Fc1/Fs;
r2 = pi*Fc2/Fs;
zp0 = (1-r1) / (1+r1);
zp2 = (1-r2) / (1+r2);
rc = pi * (Fc1 + Fc2) / 2 / Fs;
zc = (1 + 1i*rc) / (1 - 1i*rc);
de = (zc - zp0) * (zc - zp2);
nu = (zc - 1) * (zc + 1);
ka = de / nu;
a1 = -(zp0 + zp2);
a2 = zp0 * zp2;
scale = abs(ka);


function [y, z1, z2] = Step(b0, b1, b2, a1, a2, s0, x, z1, z2)
y = b0 * s0 * x + z1;
z1 = b1 * s0 * x - a1 * y + z2;
z2 = b2 * s0 * x - a2 * y;

function y = CustomFilter(a1, a2, scale, x)
y = x;
for ch = 1:size(x,2)
    z1 = 0;
    z2 = 0;
    for t = 1:size(x,1)
        [y(t,ch), z1, z2] = Step(1, 0, -1, a1, a2, scale, x(t,ch), z1, z2);
    end
end
