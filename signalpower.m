% Calculates the power of a real time series signal at a particular
% frequency range (or all range if not specified)

function P = signalpower (signal, Fs, frange, NFFT)
if ~exist('frange','var') || isempty(frange)
    frange = [-inf inf];
end
Min_NFFT = Fs/min(diff(frange,[],2));
L = size(signal,1);
if ~exist('NFFT','var') || isempty(NFFT) || ~isnumeric(NFFT)
    %NFFT = 2^nextpow2(max(L,Min_NFFT));
    NFFT = 8*2^nextpow2(Min_NFFT);
end
X = fft(signal,NFFT,1);
d = Fs/NFFT;
%f = [0:d:Fs/2, -Fs/2+d:d:-d];
f = 0:d:Fs/2;
P = zeros(size(frange,1), size(signal,2));
for i = 1:size(frange,1)
    %Y = X((f >= frange(i,1) & f <= frange(i,2)) | (f < -frange(i,1) & f > -frange(i,2)),:);
    Y = X(f >= frange(i,1) & f <= frange(i,2),:);
    P(i,:) = 2 * sum(abs(Y).^2,1) / NFFT / L;
end
return

% X = fft(signal,[],1);
% if ~exist('frange','var') || size(frange,2) ~= 2
%     Y = X;
%     P = sum(abs(Y).^2) / size(signal,1).^2;
% else
%     P = zeros(size(frange,1), size(signal,2));
%     for i = 1:size(frange,1)
%         Y = isolatefreq(X,Fs,frange(i,:),'pass');
%         P(i,:) = sum(abs(Y).^2,1) / size(signal,1).^2;
%     end
% end


% function Y = isolatefreq (Y, Fs, Fc, Type)
% 
% Ratio = (size(Y,1)) / Fs;
% Fs = (size(Y,1));
% Fc = Fc * Ratio;
% 
% % from DC (0) to Nyquist frequency (Fs/2)
% Fall = 1:floor(Fs/2+1);
% 
% RS = Y(Fall,:);
% 
% LS = flipud(Y(setdiff(1:Fs,Fall),:));
% 
% RSFall = Fall;
% LSFall = fliplr(-(setdiff(1:Fs,Fall)-Fs-1));
% 
% % Sanity check
% for i = 1:length(Fc)
%     if Fc(i) >= Fall(end)
%         Fc(i) = Fall(end)-1;
%     elseif Fc(i) < 0
%         Fc(i) = 0;
%     end
% end
% 
% if strcmp(Type,'pass')
%     Type = 'band';
% end
% 
% % Filtering
% switch Type
%     case 'low'
%         RSFstop = (RSFall > Fc+1);
%         LSFstop = (LSFall > Fc);
%         RS(RSFstop,:) = 0;
%         LS(LSFstop,:) = 0;
%     case 'high'
%         RSFstop = (RSFall < Fc+1);
%         LSFstop = (LSFall < Fc);
%         RS(RSFstop,:) = 0;
%         LS(LSFstop,:) = 0;
%     case 'band'
%         RSFstop = union(find(RSFall < Fc(1)+1),find(RSFall > Fc(2)+1));
%         LSFstop = union(find(LSFall < Fc(1)),find(LSFall > Fc(2)));
%         RS(RSFstop,:) = 0;
%         LS(LSFstop,:) = 0;
%     case 'stop'
%         RSFstop = intersect(find(RSFall >= Fc(1)+1),find(RSFall <= Fc(2)+1));
%         LSFstop = intersect(find(LSFall >= Fc(1)),find(LSFall <= Fc(2)));
%         RS(RSFstop,:) = 0;
%         LS(LSFstop,:) = 0;
%     otherwise
%         error('Unknown filter type (4th parameter). Must be "low", "high", "band", or "stop"');
% end
% 
% Y = cat(1,RS,flipud(LS));
% 
