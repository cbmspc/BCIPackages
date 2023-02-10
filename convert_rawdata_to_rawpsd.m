function [rawpxx, fxx] = convert_rawdata_to_rawpsd (rawdataOReegdata, Fs, Nfft)
% Converts rawdata (chan x time x trial) or eegdata (cell of time x chan)
% into pxx and fxx pair (chan x freq x trial)

if isempty(rawdataOReegdata)
    error('rawdata is empty');
end

if ~iscell(rawdataOReegdata) && size(rawdataOReegdata,3) > 10
    rawdata = rawdataOReegdata;
elseif iscell(rawdataOReegdata) && length(rawdataOReegdata) > 10
    eegdata = rawdataOReegdata;
else
    %warning('Cannot determine whether the input is rawdata (3d matrix) or eegdata (1d cell). Assuming it is rawdata');
    rawdata = rawdataOReegdata;
end

if ~exist('Nfft','var') || isempty(Nfft)
    Nfft = 1024;
end

if exist('rawdata','var')
    [Nch, ~, Ntr] = size(rawdata);
elseif exist('eegdata','var')
    Nch = size(eegdata{1},2);
    Ntr = length(eegdata);
end

Nfp = Nfft/2+1;
if floor(Nfft/2) ~= Nfft/2
    Nfp = (Nfft+1)/2;
end
rawpxx = nan(Nch, Nfp, Ntr);
if Ntr > 100
    fprintf(' Converting time domain signals to PSD ..');
end

% With 4 threads, parfor has 60% speed increase.
if exist('rawdata','var')
    [~, fxx] = pwelch(rawdata(1,:,1).', [], [], Nfft, Fs);
    parfor tr = 1:Ntr
        for ch = 1:Nch
            [pxx, ~] = pwelch(rawdata(ch,:,tr).', [], [], Nfft, Fs);
            rawpxx(ch,:,tr) = pxx.';
        end
        
    end
elseif exist('eegdata','var')
    [~, fxx] = pwelch(eegdata{1}(:,1), [], [], Nfft, Fs);
    parfor tr = 1:Ntr
        for ch = 1:Nch
            [pxx, ~] = pwelch(eegdata{tr}(:,ch), [], [], Nfft, Fs);
            rawpxx(ch,:,tr) = pxx.';
        end
    end
end
if Ntr > 100
    fprintf(' done \n');
end
