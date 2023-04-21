% 1) Run a band-pass filter on a continuous uncut time series data (eeg)
% 2) Segment the filtered data into segments within the specified time
%    ranges (specified in BoundsSec)
%    BoundsSec is (Number of segments, 2), where the first column is the
%    start time of each segment, and the second column is the end time of 
%    each segment.
% (1 and 2 simulate a filter bank that does not reset)
% 3) Calculate the frequency band power within each segment. There is one
%    power value per channel per segment. 
%    frange is (Number of frequency bands, 2)
% 4) Output format: chan x band x segment

function rawsigpow = convert_eeg_to_rawsigpow (eeg, Fs, frange, BoundsSec, opts)

filtorder = 4;
filtname = 'buttercausal';

if exist('opts','var') && isstruct(opts)
    if isfield(opts, 'filtorder') && ~isempty(opts.filtorder)
        filtorder = opts.filtorder;
    end
    if isfield(opts, 'filtname') && ~isempty(opts.filtname)
        filtname = opts.filtname;
    end
end

eeg = detrend(eeg, 'constant');
nseg = size(BoundsSec,1);
nchan = size(eeg,2);
nband = size(frange,1);
rawsigpow = nan(nchan,nband,nseg);

for band = 1:nband
    
    % Filter then segment
    
    % Filter
    fb = [1 1];
    if frange(band,1) <= 0
        fb(1) = 0;
    end    
    if frange(band,2) >= Fs/2
        fb(2) = 0;
    end    
    if fb(1) && fb(2)
        filtparm = [frange(band,1:2) filtorder];
        filttype = 'pass';
        bandeeg = freqfilter(eeg, Fs, filtparm, filttype, filtname);
    elseif fb(1) && ~fb(2)
        filtparm = [frange(band,1) filtorder];
        filttype = 'high';
        bandeeg = freqfilter(eeg, Fs, filtparm, filttype, filtname);
    elseif ~fb(1) && fb(2)
        filtparm = [frange(band,2) filtorder];
        filttype = 'low';
        bandeeg = freqfilter(eeg, Fs, filtparm, filttype, filtname);
    else
        bandeeg = eeg;
    end
    
    % Segment
    
    
    for seg = 1:nseg
        ta = BoundsSec(seg,1);
        tb = BoundsSec(seg,2);
        ka = round(ta*Fs+1);
        kb = round(tb*Fs);
        segbandeeg = bandeeg(ka:kb,:);
        rawsigpow(:,band,seg) = mean(segbandeeg.^2,1) / (frange(band,2)-frange(band,1));
    end
    
    
end