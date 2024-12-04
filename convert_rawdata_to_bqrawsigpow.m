function [rawsigpow, frange, opts] = convert_rawdata_to_bqrawsigpow (rawdataOReegdata, Fs, frange, opts)
% Converts rawdata (chan x time x trial) or eegdata (cell of time x chan)
% into rawsigpow (chan x freq x trial)

if isempty(rawdataOReegdata)
    error('rawdata is empty');
end

if ~exist('opts', 'var') || ~isstruct(opts)
    opts = struct;
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
    [Nchan, Ntimepertrial, Ntrial] = size(rawdata);
    Ntime = Ntimepertrial * Ntrial;
elseif exist('eegdata','var')
    Nchan = size(eegdata{1},2);
    cout = cellfun(@size,eegdata,'UniformOutput',false);
    sizes = cat(1,cout{:});
    Ntime = sizes(:,1);
    Ntrial = length(eegdata);
end

Nfreq = size(frange,1);

%rawsigpow = nan(Nchan, Nfreq, Ntrial);
%fprintf(' Converting time domain signals to spectral powers .. %3i%% done', 0);
if Ntrial > 100
    fprintf(' Converting time domain signals to spectral powers ..');
end

% Concatenate to continuous time series
signal = nan(sum(Ntime),Nchan); 
TimeRange = nan(Ntrial,2);

if exist('rawdata','var')
    for tr = 1:Ntrial
        signal((tr-1)*Ntime+(1:Ntime),:) = rawdata(:,:,tr).';
        TimeRange(tr,:) = [tr-1, tr].*Ntime/Fs;
    end
elseif exist('eegdata','var')
    for tr = 1:Ntrial
        a = sum(Ntime(1:tr-1)) + 1;
        b = sum(Ntime(1:tr-1)) + Ntime(tr);
        signal(a:b,:) = eegdata{tr};
        TimeRange(tr,:) = [a b]/Fs;
    end
end
% Run bqsignalpower
[rawsigpow, opts] = bqsignalpower(signal, Fs, frange, TimeRange, opts);


% % With 4 threads, parfor has 60% speed increase.
% if exist('rawdata','var')
%     parfor tr = 1:Ntrial
%         rawsigpow(:,:,tr) = signalpower(rawdata(:,:,tr).', Fs, frange).';
%     end
% elseif exist('eegdata','var')
%     parfor tr = 1:Ntrial
%         rawsigpow(:,:,tr) = signalpower(eegdata{tr}, Fs, frange).';
%     end
% end
if Ntrial > 100
    fprintf(' done \n');
end