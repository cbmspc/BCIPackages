% 1) Runs a band-pass filter on a continuous uncut time series data
% 2) Segment the filtered data into equal sized segments within the
%    specified time range   
% (1 and 2 simulate a filter bank that does not reset)
% 3) Calculate the power within each segment. There is one power value per
%    channel per segment
% 4) Output format: chan x band x segment
% Note: All segments ideally have equal size. Possible gap between last and tend.
%   It is possible to have unequal segment sizes by setting tsegsize to
%   larger than trange because a minimum of one segment will always be
%   created for each trange.
%
%
% Example 1: To convert the first epoch of rawdata (chan x time x epoch)
% that is 10 seconds long per epoch into six segments:
%   * excluded: first 0.5 seconds in the epoch
%   * segment #1: 0.5 - 2.0 seconds into the one epoch
%   * segment #2: 2.0 - 3.5 seconds into the one epoch
%   * segment #3: 3.5 - 5.0 seconds into the one epoch
%   * segment #4: 5.0 - 6.5 seconds into the one epoch
%   * segment #5: 6.5 - 8.0 seconds into the one epoch
%   * segment #6: 8.0 - 9.5 seconds into the one epoch
%   * excluded: last 0.5 seconds in the epoch
%
%  Fs = 256
%  frange = [(8:2:38).' (10:2:40).']
%  trange = [0.5 9.5]
%  tsegsize = 1.5
%  opts = struct  (all opts below are optional)
%  opts.filtorder = 4
%  opts.filtname = 'buttercausal'
%  opts.detrend_n = 'constant'
%       %For backward compatibility, buttercausal and constant are the default
%  segpow = convert_timeseriesdata_to_segmentedpowervalues(rawdata(:,:,1).', Fs, frange, trange, tsegsize, opts);
%    %Transposed because rawdata(:,:,1) is chan x time x 1 but this function expects time x chan.
%    %Don't forget to take the logarithm if appropriate
%
%
%  Example 2: Convert a time series data "signal" into segment powers
%   Reminder: signal is (time x chan)
%   trange = [ 0 20
%             30 35
%             35 40
%             50 60]; % four epochs in the signal. They don't 
%                     % have to have the same length
%   frange = [ 4  8
%            [11 25 ]; % In Hz, the frequency bands
%   tsegsize = 2; %seconds
%   Example syntax: 
%      opts.filtname = 'butter'
%      [segpow, segtimesamplesrange, segepochindex] = convert_timeseriesdata_to_segmentedpowervalues(signal, Fs, frange, trange, tsegsize, opts);
%      This will create a (chan x band x seg) output. The original epoch
%      index (rows in trange) where each segment belonged is kept in
%      segepochindex. The original time sample indices where each segment
%      is cut from is kept in segtimesamplesrange.
%      If an epoch is shorter than tsegsize, the segment for that epoch
%      will be just that epoch.
%
% See the wrapper function convert_rawdata_to_segmentedpowervalues

function [segpow, segtimesamplesrange, segepochindex] = convert_timeseriesdata_to_segmentedpowervalues (data, Fs, frange, trange, tsegsize, opts)

if numel(Fs) > 1
    Fs = Fs(1);
end

if exist('opts', 'var') && isscalar(opts) && isnumeric(opts)
    % Old style input for compatibility
    filtorder = opts;
    opts = struct;
    opts.filtorder = filtorder;
    clear filtorder
end

if ~exist('opts', 'var') || isempty(opts) || ~isstruct(opts)
    opts = struct;
end


if ~isfield(opts,'filtorder') || numel(opts.filtorder) ~= 1 || ~isnumeric(opts.filtorder)
    opts.filtorder = 4;
end

if ~isfield(opts,'filtname') || isempty(opts.filtname) || ~ischar(opts.filtname)
    opts.filtname = 'buttercausal';
end

if ~isfield(opts,'detrend_n') || isempty(opts.detrend_n)
    opts.detrend_n = 'constant';
end

output_segtimesamplesrange_only = false;
if isfield(opts, 'output_segtimesamplesrange_only') && isscalar(opts.output_segtimesamplesrange_only) && opts.output_segtimesamplesrange_only
    output_segtimesamplesrange_only = true;
end


data = detrend(data, opts.detrend_n);
nepoch = size(trange,1);
trange(trange(:,2) > size(data,1)/Fs,2) = size(data,1)/Fs;
nseg = floor(diff(trange(:,1:2),[],2)/tsegsize);
nseg(nseg==0) = 1;
segpow = nan(size(data,2),size(frange,1),sum(nseg));
segtimesamplesrange = nan(sum(nseg),2);
segepochindex = nan(sum(nseg),1);


frange(frange(:,1) < 0,1) = 0;
frange(frange(:,2) > Fs/2,2) = Fs/2;


for band = 1:size(frange,1)
    % Filter then segment
    % Filter
    fb = [1 1];
    if frange(band,1) <= 0
        fb(1) = 0;
    end    
    if frange(band,2) >= Fs/2
        fb(2) = 0;
    end    
    if output_segtimesamplesrange_only
        banddata = data;
    elseif fb(1) && fb(2)
        filtparm = [frange(band,1:2) opts.filtorder];
        filttype = 'pass';
        banddata = freqfilter(data, Fs, filtparm, filttype, opts.filtname);
    elseif fb(1) && ~fb(2)
        filtparm = [frange(band,1) opts.filtorder];
        filttype = 'high';
        banddata = freqfilter(data, Fs, filtparm, filttype, opts.filtname);
    elseif ~fb(1) && fb(2)
        filtparm = [frange(band,2) opts.filtorder];
        filttype = 'low';
        banddata = freqfilter(data, Fs, filtparm, filttype, opts.filtname);
    else
        banddata = data;
    end

    kseg = 0;
    % Epoch
    for epoch = 1:nepoch
        tstart = trange(epoch,1);
        tend = trange(epoch,2);

        % Segment
        for seg = 1:nseg(epoch)
            kseg = kseg + 1;
            ta = tstart + (seg-1)*tsegsize;
            if isnan(ta)
                ta = tstart;
            end
            tb = ta + tsegsize;
            if tb > tend
                tb = tend;
            end
            ka = round(ta*Fs+1);
            kb = round(tb*Fs);
            segtimesamplesrange(kseg,:) = [ka kb];
            segepochindex(kseg) = epoch;
            if ~output_segtimesamplesrange_only
                segbanddata = banddata(ka:kb,:);
                %2024-12-16: output is now calibrated to a sine wave
                %segpow(:,band,kseg) = mean(segbanddata.^2,1) / (frange(band,2)-frange(band,1));
                segpow(:,band,kseg) = mean(segbanddata.^2,1);
            end
        end

    end

    if output_segtimesamplesrange_only
        break
    end
    
    
end

