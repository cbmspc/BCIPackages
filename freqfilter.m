% Frequency filter based on discrete fourier transform
%
% yout = freqfilter (y, Fs, Fc, Type, Method)
% 
% y: Time domain signal (time x channel)
%
% Fs: Sampling frequency (even number)
%
% Fc: Cutoff frequency or frequencies [lowbound highbound].
%     Additional options available when Method = 'butter':
%      For low/high pass:
%       Fc = [freq_Pass freq_Stop dB_pass dB_stop]
%              example: [40 50 1 60]
%         or [freq_Cutoff Order] (Default 3 dB roll off)
%         or [freq_Cutoff] (Default Order = 1)
%      For band pass/stop:
%       Fc = [freq_PassLeft freq_PassRight freq_StopLeft freq_StopRight ...
%             dB_pass dB_stop]
%              example: [0.02 40 0.001 55 1 60]
%         or [freq_CutoffLeft freq_CutoffRight Order] (Default 3 dB)
%         or [freq_CutoffLeft freq_CutoffRight] (Default Order = 1)
%
% Type: 
%        'low' = low pass, passes f <= Fc
%       'high' = high pass, passes Fc <= f
%       'band' = band pass, passes Fc(1) <= f <= Fc(2)
%       'stop' = band stop, stops Fc(1) <= f <= Fc(2)
%
% Method:
%      default = fft
%        'fft' = converts signals to freq domain via fft, zeros stop bands,
%                and converts back to time domain. This is a non-linear
%                method.
%     'butter' = uses Butterworth filter and filtfilt
%     'buttercausal' = uses Butterworth filter and filter
%

function [yout, FilterInfo] = freqfilter (y, Fs, Fc, Type, Method)

if isempty(who('y'))
    error('Time signal (1st parameter) not specified.');
end

if isempty(who('Fs'))
    error('Sampling frequency (2nd parameter) not specified.');
end

if isempty(who('Fc'))
    error('Cutoff frequency (3rd parameter) not specified.');
end

if isempty(who('Type'))
    error('Filter type (4th parameter) not specified.');
end

if ~exist('Method', 'var')
    Method = 'fft';
end

switch Method
    case 'butter'
        [yout, FilterInfo] = butter_filter(y, Fs, Fc, Type, [], @filtfilt);
    case 'buttercausal'
        [yout, FilterInfo] = butter_filter(y, Fs, Fc, Type, [], @filter);
    otherwise
        [yout, FilterInfo] = fft_filter(y, Fs, Fc, Type);
end


function [yout, FilterInfo] = butter_filter(y, Fs, Fc, Type, FilterInfo, FilterCommand)
% Takes four parameters per cutoff: Pass band, Stop band, dB pass, dB stop
if ~exist('FilterInfo','var') || isempty(FilterInfo)
    FilterInfo = [];
    FilterInfo.FilterMethod = 'butter';
    FilterInfo.FilterCommand = FilterCommand;
    FilterInfo.ButterOrder = [0 0];
    FilterInfo.ButterUnstable = [0 0];
    FilterInfo.FilterA = {};
    FilterInfo.FilterB = {};
end

Type = lower(Type);
if strcmp(Type,'pass')
    Type = 'band';
end

if ~isfield(FilterInfo,'FilterType')
    FilterInfo.FilterType = Type;
end

if strcmp(Type,'low') || strcmp(Type,'high')
    if length(Fc) == 4
        % Specified passband, stopband, dB of pass, dB of stop
        F.pass = Fc(1);
        F.stop = Fc(2);
        dB.pass = Fc(3);
        dB.stop = Fc(4);
        [N Wn] = buttord(F.pass/(Fs/2), F.stop/(Fs/2), dB.pass, dB.stop);
    elseif length(Fc) == 2
        % Specified only cutoff and order
        Wn = Fc(1)/(Fs/2);
        N = Fc(2);
    elseif length(Fc) == 1
        % Specified only cutoff. default to lowest order
        Wn = Fc(1)/(Fs/2);
        N = 1;
    else
        % Unknown specs. Defaults to lowest order
        Wn = Fc(1)/(Fs/2);
        N = 1;
        warning('FREQFILTER:WRONGFCPARM1','This is a high/low pass filter. Fc needs to have either 4, 2, or 1 elements. Defaulting to order N=1');
    end
elseif strcmp(Type,'band')
    if length(Fc) == 6
        % Specified passbands, stopbands, dB of pass, dB of stop
        [yout FilterInfo] = butter_filter(y, Fs, Fc([1 3 5 6]), 'high', FilterInfo, FilterCommand);
        [yout FilterInfo] = butter_filter(yout, Fs, Fc([2 4 5 6]), 'low', FilterInfo, FilterCommand);
    elseif length(Fc) == 3
        % Specified only cutoffs and order
        [yout FilterInfo] = butter_filter(y, Fs, Fc([1 3]), 'high', FilterInfo, FilterCommand);
        [yout FilterInfo] = butter_filter(yout, Fs, Fc([2 3]), 'low', FilterInfo, FilterCommand);
    elseif length(Fc) == 2
        % Specified only cutoffs. default to lowest order
        [yout FilterInfo] = butter_filter(y, Fs, Fc(1), 'high', FilterInfo, FilterCommand);
        [yout FilterInfo] = butter_filter(yout, Fs, Fc(2), 'low', FilterInfo, FilterCommand);
    else
        % Unknown specs. Default to lowest order
        [yout FilterInfo] = butter_filter(y, Fs, Fc(1), 'high', FilterInfo, FilterCommand);
        [yout FilterInfo] = butter_filter(yout, Fs, Fc(2), 'low', FilterInfo, FilterCommand);
        warning('FREQFILTER:WRONGFCPARM2','This is a band pass filter. Fc needs to have either 6, 3, or 2 elements. Defaulting to order 2N=2');
    end
elseif strcmp(Type,'stop')
    if length(Fc) == 6
        % Specified passbands, stopbands, dB of pass, dB of stop
        [yout1 FilterInfo] = butter_filter(y, Fs, Fc([1 3 5 6]), 'low', FilterInfo, FilterCommand);
        [yout2 FilterInfo] = butter_filter(y, Fs, Fc([2 4 5 6]), 'high', FilterInfo, FilterCommand);
        yout = yout1 + yout2;
    elseif length(Fc) == 3
        % Specified only cutoffs and order
        [yout1 FilterInfo] = butter_filter(y, Fs, Fc([1 3]), 'low', FilterInfo, FilterCommand);
        [yout2 FilterInfo] = butter_filter(y, Fs, Fc([2 3]), 'high', FilterInfo, FilterCommand);
        yout = yout1 + yout2;
    elseif length(Fc) == 2
        % Specified only cutoffs. default to lowest order
        [yout1 FilterInfo] = butter_filter(y, Fs, Fc(1), 'low', FilterInfo, FilterCommand);
        [yout2 FilterInfo] = butter_filter(y, Fs, Fc(2), 'high', FilterInfo, FilterCommand);
        yout = yout1 + yout2;
    else
        % Unknown specs. Default to lowest order
        [yout1 FilterInfo] = butter_filter(y, Fs, Fc(1), 'low', FilterInfo, FilterCommand);
        [yout2 FilterInfo] = butter_filter(y, Fs, Fc(2), 'high', FilterInfo, FilterCommand);
        yout = yout1 + yout2;
        warning('FREQFILTER:WRONGFCPARM3','This is a band stop filter. Fc needs to have either 6, 3, or 2 elements. Defaulting to order 2N=2');
    end
% elseif strcmp(Type,'band') || strcmp(Type,'stop')
%     if length(Fc) == 6
%         % Specified passbands, stopbands, dB of pass, dB of stop
%         F.pass = [Fc(1) Fc(2)];
%         F.stop = [Fc(3) Fc(4)];
%         dB.pass = Fc(5);
%         dB.stop = Fc(6);
%         [N Wn] = buttord(F.pass/(Fs/2), F.stop/(Fs/2), dB.pass, dB.stop);
%     elseif length(Fc) == 3
%         % Specified only cutoffs and order
%         Wn = Fc(1:2)/(Fs/2);
%         N = Fc(3);
%     elseif length(Fc) == 2
%         % Specified only cutoffs. default to lowest order
%         Wn = Fc/(Fs/2);
%         N = 1;
%     else
%         warning('This is a band pass/stop filter. Fc needs to have either 6, 3, or 2 elements.');
%     end
end

switch Type
    case 'low'
        [B A] = butter(N, Wn, 'low');
        FilterInfo.ButterOrder(1) = N;
        h1=dfilt.df2(B,A);
        if ~isstable(h1)
            warning('FREQFILTER:UNSTABLE1',['This Butterworth low pass filter of order ' num2str(N) ' is not stable.']);
            FilterInfo.ButterUnstable(1) = 1;
        end
        yout = FilterCommand(B, A, y);
    case 'high'
        [B A] = butter(N, Wn, 'high');
        FilterInfo.ButterOrder(2) = N;
        h1=dfilt.df2(B,A);
        if ~isstable(h1)
            warning('FREQFILTER:UNSTABLE2',['This Butterworth high pass filter of order ' num2str(N) ' is not stable.']);
            FilterInfo.ButterUnstable(2) = 1;
        end
        yout = FilterCommand(B, A, y);
    %case 'band'
    %    [B A] = butter(N, Wn);
    %case 'stop'
    %    [B A] = butter(N, Wn, 'stop');
end
if exist('A','var') && exist('B','var')
    FilterInfo.FilterA{length(FilterInfo.FilterA)+1} = A;
    FilterInfo.FilterB{length(FilterInfo.FilterB)+1} = B;
end

function [yout, FilterInfo] = fft_filter (y, Fs, Fc, Type, FilterInfo)

if ~exist('FilterInfo','var')
    FilterInfo = [];
    FilterInfo.FilterMethod = 'fft';
    FilterInfo.FilterCommand = '';
end

Ratio = (size(y,1)) / Fs;
Fs = (size(y,1));
Fc = Fc * Ratio;

Y = fft(y,[],1);

% from DC (0) to Nyquist frequency (Fs/2)
Fall = 1:floor(Fs/2+1);

RS = Y(Fall,:);

LS = flipud(Y(setdiff(1:Fs,Fall),:));

RSFall = Fall;
LSFall = fliplr(-(setdiff(1:Fs,Fall)-Fs-1));

% Sanity check
for i = 1:length(Fc)
    if Fc(i) >= Fall(end)
        Fc(i) = Fall(end)-1;
    elseif Fc(i) < 0
        Fc(i) = 0;
    end
end

Type = lower(Type);
if strcmp(Type,'pass')
    Type = 'band';
end

if ~isfield(FilterInfo,'FilterType')
    FilterInfo.FilterType = Type;
end

% Filtering
switch Type
    case 'low'
        %RSFstop = find(RSFall > Fc+1);
        %LSFstop = find(LSFall > Fc);
        %RS(RSFstop,:) = 0;
        %LS(LSFstop,:) = 0;
        RS(RSFall > Fc+1,:) = 0;
        LS(LSFall > Fc,:) = 0;
    case 'high'
        %RSFstop = find(RSFall < Fc+1);
        %LSFstop = find(LSFall < Fc);
        %RS(RSFstop,:) = 0;
        %LS(LSFstop,:) = 0;
        RS(RSFall < Fc+1,:) = 0;
        LS(LSFall < Fc,:) = 0;
    case 'band'
        RSFstop = union(find(RSFall < Fc(1)+1),find(RSFall > Fc(2)+1));
        LSFstop = union(find(LSFall < Fc(1)),find(LSFall > Fc(2)));
        RS(RSFstop,:) = 0;
        LS(LSFstop,:) = 0;
    case 'stop'
        RSFstop = intersect(find(RSFall >= Fc(1)+1),find(RSFall <= Fc(2)+1));
        LSFstop = intersect(find(LSFall >= Fc(1)),find(LSFall <= Fc(2)));
        RS(RSFstop,:) = 0;
        LS(LSFstop,:) = 0;
    otherwise
        error('Unknown filter type (4th parameter). Must be "low", "high", "band", or "stop"');
end

Y = cat(1,RS,flipud(LS));

yout = ifft(Y,[],1);

