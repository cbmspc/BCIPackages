% 1) Runs a band-pass filter on a continuous uncut time series data
% 2) Segment the filtered data into equal sized segments within the
%    specified time range   
% (1 and 2 simulate a filter bank that does not reset)
% 3) Calculate the power within each segment. There is one power value per
%    channel per segment
% 4) Output format: chan x band x segment
% Note: All segments have equal size. Possible gap between last and tend.
%
%
% Example: To convert the first epoch of rawdata (chan x time x epoch) that
% is 10 seconds long per epoch into six segments:
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


function segpow = convert_timeseriesdata_to_segmentedpowervalues (data, Fs, frange, trange, tsegsize, opts)

if numel(Fs) > 1
    Fs = Fs(1);
end

if exist('opts', 'var') && numel(opts) == 1 && isnumeric(opts)
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

tstart = trange(1);
tend = trange(2);

if tend > size(data,1)/Fs
    tend = size(data,1)/Fs;
end

data = detrend(data, opts.detrend_n);

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
    if fb(1) && fb(2)
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
    
    % Segment
    nseg = floor((tend-tstart)/tsegsize);
    for seg = nseg:-1:1
        ta = tstart + (seg-1)*tsegsize;
        tb = ta + tsegsize;
        ka = round(ta*Fs+1);
        kb = round(tb*Fs);
        segbanddata = banddata(ka:kb,:);
        segpow(:,band,seg) = mean(segbanddata.^2,1) / (frange(band,2)-frange(band,1));
    end
    
    
end


