function [rawsigpow, frange] = convert_rawdata_to_rawsigpow (rawdataOReegdata, Fs, frange, NFFT, method)
% Converts rawdata (chan x time x trial) or eegdata (cell of time x chan)
% into rawsigpow (chan x freq x trial)

if isempty(rawdataOReegdata)
    error('rawdata is empty');
end

if ~exist('NFFT','var') || isempty(NFFT)
    NFFT = [];
end

if ~exist('method','var') || isempty(method)
    method = '';
end

if ~iscell(rawdataOReegdata) && size(rawdataOReegdata,3) > 10
    rawdata = rawdataOReegdata;
elseif iscell(rawdataOReegdata) && length(rawdataOReegdata) > 10
    eegdata = rawdataOReegdata;
else
    %warning('Cannot determine whether the input is rawdata (3d matrix) or eegdata (1d cell). Assuming it is rawdata');
    rawdata = rawdataOReegdata;
end

if ~exist('frange','var') || isempty(frange)
    frange = [[0:2:38];[2:2:40]]';
end

if exist('rawdata','var')
    [Nch tmp Ntr] = size(rawdata);
elseif exist('eegdata','var')
    Nch = size(eegdata{1},2);
    Ntr = length(eegdata);
end

Ntp = size(frange,1);

rawsigpow = nan(Nch, Ntp, Ntr);
%fprintf(' Converting time domain signals to spectral powers .. %3i%% done', 0);
% if Ntr > 100
%     fprintf(' Converting time domain signals to spectral powers ..');
% end

if exist('rawdata','var')
    for tr = 1:Ntr
        rawsigpow(:,:,tr) = signalpower(rawdata(:,:,tr).', Fs, frange, NFFT).';
        %fprintf('\b\b\b\b\b\b\b\b\b%3i%% left', round(tr/Ntr*100));
    end
elseif exist('eegdata','var')
    for tr = 1:Ntr
        rawsigpow(:,:,tr) = signalpower(eegdata{tr}, Fs, frange, NFFT).';
        %fprintf('\b\b\b\b\b\b\b\b\b%3i%% left', round(tr/Ntr*100));
    end
end
% if Ntr > 100
%     fprintf(' done \n');
% end
