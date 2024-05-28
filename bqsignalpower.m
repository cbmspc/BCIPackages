% Calculates the spectral power of time segments in a real time series
% signal at particular frequency ranges using causal biquad bandpass filter
% bank

function [signal_power, opts, signal_filtered_in_each_freqbin] = bqsignalpower (signal, Fs, FreqRange, TimeRange, opts)
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

if nargout >= 3
    % User wants signal_filtered_in_each_bin
    signal_filtered_in_each_freqbin = nan(SigLen,Nchan,Nfreqbin);
end

%silent = 1;

%wh = waitbar(0, ['bqsignalpower Band 0/' num2str(Nfreqbin)]);



for b = 1:Nfreqbin
    Fc1 = FreqRange(b, 1);
    Fc2 = FreqRange(b, 2);
    [a1, a2, scale] = Calculate4thOrderBandpassFilterCoefficients(Fs, Fc1, Fc2);
    a0 = a2*0+1;
    opts.filtcoefs{b} = [a1, a2, scale];
    b0 = 1*scale;
    b1 = 0*scale;
    b2 = -1*scale;
    if opts.ReflectSignal
        %signal_reflected_filtered = filter([b0(1) b1(1) b2(1)], [a0(1) a1(1) a2(1)], signal_reflected);
        signal_reflected_filtered = CustomFilter(a1(1), a2(1), scale(1), signal_reflected);
        signal_filtered = signal_reflected_filtered(opts.ReflectSigLen+1:end,:);
        for i = 2:opts.FilterIter
            %signal_filtered = filter([b0(i) b1(i) b2(i)], [a0(i) a1(i) a2(i)], signal_filtered);
            signal_filtered = CustomFilter(a1(i), a2(i), scale(i), signal_filtered);
        end
    else
        %signal_filtered = filter([b0(1) b1(1) b2(1)], [a0(1) a1(1) a2(1)], signal);
        signal_filtered = CustomFilter(a1(1), a2(1), scale(1), signal);
        for i = 2:opts.FilterIter
            %signal_filtered = filter([b0(i) b1(i) b2(i)], [a0(i) a1(i) a2(i)], signal_filtered);
            signal_filtered = CustomFilter(a1(i), a2(i), scale(i), signal_filtered);
        end
    end

    if exist('signal_filtered_in_each_bin','var')
        signal_filtered_in_each_freqbin(:,:,b) = signal_filtered;
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
