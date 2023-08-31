% Calculates the power of a time series signal in each frequency range, or
% all range if frange is not specified
% signal must be time x chan, or a tall vector if there is only 1 channel
% Fs is the sampling frequency in Hz
% frange is optional
% NFFT is optional
% method can be 'fft' or 'pwelch'
% pwelch is the default, as it requires less nfft


function P = signalpower (signal, Fs, frange, NFFT, method)
if ~exist('frange','var') || isempty(frange) || ~isnumeric(frange) || size(frange,2) ~= 2
    %frange = [-inf inf];
    frange = [0 Fs/2];
end
if ~exist('method','var') || isempty(method)
    method = '';
end
if ~exist('NFFT','var') || ~isnumeric(NFFT)
    NFFT = [];
end

if size(signal,1) <= 2 && size(signal,2) > 10
    % probably need to transpose
    signal = signal.';
    warning('Signal width > height. Transposed signal.');
end

% De-mean the signal for better estimation
mu = mean(signal,1);
signal = detrend(signal);
L = size(signal,1);

switch method
    case 'fft'
        P = via_fft_method(signal, Fs, frange, NFFT, L);
    case 'pwelch'
        P = via_pwelch_twosided_method(signal, Fs, frange, NFFT, L);
    case 'pwelchonesided'
        P = via_pwelch_onesided_method(signal, Fs, frange, NFFT, L);
    case 'pwelchtwosided'
        P = via_pwelch_twosided_method(signal, Fs, frange, NFFT, L);
    case 'pwelchcentered'
        P = via_pwelch_centered_method(signal, Fs, frange, NFFT, L);
    otherwise
        P = via_pwelch_twosided_method(signal, Fs, frange, NFFT, L);
end

% add DC power back
for i = 1:size(frange,1)
    if frange(i,1) == 0
        P(i,:) = P(i,:) + mu.^2;
    end
end

return


function P = via_fft_method (signal, Fs, frange, NFFT, L)
X = fft(signal,NFFT,1);
NFFT = size(X,1);
d = Fs/NFFT;
f = 0:d:Fs/2;
P = zeros(size(frange,1), size(signal,2));
for i = 1:size(frange,1)
    Y = X(f >= frange(i,1) & f < frange(i,2),:);
    P(i,:) = 2 * sum(abs(Y).^2,1) / NFFT / L;
end
return


function P = via_pwelch_twosided_method (signal, Fs, frange, NFFT, ~)
[X, f] = pwelch(signal, [], [], NFFT, Fs, 'twosided', 'psd');
NFFT = length(f);
P = zeros(size(frange,1), size(signal,2));
for i = 1:size(frange,1)
    Y = X(f >= frange(i,1) & f < frange(i,2),:);
    P(i,:) = 2 * sum(Y,1) / NFFT * Fs;
end
return


function P = via_pwelch_onesided_method (signal, Fs, frange, NFFT, ~)
[X, f] = pwelch(signal, [], [], NFFT, Fs, 'onesided', 'psd');
NFFT = (length(f)-1)*2;
P = zeros(size(frange,1), size(signal,2));
for i = 1:size(frange,1)
    Y = X(f >= frange(i,1) & f < frange(i,2),:);
    P(i,:) = sum(Y,1) / NFFT * Fs;
end
return


function P = via_pwelch_centered_method (signal, Fs, frange, NFFT, ~)
[X, f] = pwelch(signal, [], [], NFFT, Fs, 'centered', 'psd');
NFFT = length(f);
P = zeros(size(frange,1), size(signal,2));
for i = 1:size(frange,1)
    Y = X(f >= frange(i,1) & f < frange(i,2),:);
    P(i,:) = 2 * sum(Y,1) / NFFT * Fs;
end
return

