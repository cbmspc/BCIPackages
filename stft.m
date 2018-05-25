% Anchor = 
%       middle
%       causal

function [pxx, fxx] = stft (x, Fs, ws, NFFT, Causal, ZeroPad)
% stft
if ~exist('Causal','var') || isempty(Causal)
    Causal = 0;
end
if ~exist('ZeroPad','var') || isempty(ZeroPad)
    ZeroPad = 0;
end

if ZeroPad
    x2 = [flipud(x(1:ws,:)); x; flipud(x(end-ws+1:end,:))];
else
    x2 = x;
end

n = size(x2,1);

pxx = nan(NFFT/2+1,n);
tmp = nan(ws,n);

if Causal
    for i = 1+ws-1:n
        tmp(:,i) = detrend(x2(i-ws+1:i));
    end
else
    for i = 1+ws/2:n-ws/2+1
        tmp(:,i) = detrend(x2(i-ws/2:i+ws/2-1));
    end
end

tmp = fft(tmp, NFFT, 1);
pxx(:,1:size(tmp,2)) = abs(tmp(1:NFFT/2+1,:));


if ZeroPad
    pxx = pxx(:,ws+1:end-ws);
end

fxx = linspace(0, Fs/2, NFFT/2+1).';
