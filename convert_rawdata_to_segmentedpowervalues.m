% [segpowdata, seglabels] = convert_rawdata_to_segmentedpowervalues (rawdata, Fs, frange, trange, tsegsize, opts)
% This is a wrapper function that calls
% convert_timeseriesdata_to_segmentedpowervalues to convert each epoch in
% rawdata (chan x time x epoch) into segments and get their band powers
% 
% Example use case:
%  "rawdata" has 60 channels, 10 seconds data per epoch sampled at 256 Hz
%  (2560 time points per epoch), and 100 epochs.
%  The corresponding "labels" numeric variable labels the epochs.
%  Want to get 2-Hz band powers from 8 to 40 Hz from each segment after
%  splitting each epoch into six equal-sized 1.5-s segments, skipping the
%  first and last 0.5 seconds.
%  Fs = 256
%  frange = [(8:2:38).' (10:2:40).']
%  trange = [0.5 9.5]
%  tsegsize = 1.5
%  opts = struct  (all opts below are optional)
%  opts.filtorder = 4
%  opts.filtname = 'butter'
%  opts.detrend_n = 'linear'
% [segpowdata, seglabels] = convert_rawdata_to_segmentedpowervalues(rawdata, labels, Fs, frange, trange, tsegsize, opts);
%  %Don't forget to take logarithm if appropriate
%

function [segpowdata, seglabels] = convert_rawdata_to_segmentedpowervalues (rawdata, labels, Fs, frange, trange, tsegsize, opts)

if ~exist('labels', 'var') || size(rawdata,3) ~= length(labels)
    error('Missing required input argument: labels must have the same length as the number of epochs in rawdata');
end

labels = labels(:);
if size(rawdata,3) ~= length(labels)
    error('Missing required input argument: labels must have the same length as the number of epochs in rawdata');
end

if ~exist('Fs', 'var') || isempty(Fs) || ~isnumeric(Fs)
    error('Missing required input argument: Fs must be the sampling rate in Hz');
end

if ~exist('frange', 'var') || isempty(frange) || size(frange,2) ~= 2 || ~isnumeric(frange)
    error('Missing required input argument: frange must be numeric matrix of size N x 2');
end

if ~exist('trange', 'var') || numel(trange) ~= 2 || ~isnumeric(trange)
    error('Missing required input argument: trange must be numeric array of length 2');
end

if ~exist('tsegsize', 'var') || numel(tsegsize) ~= 1 || ~isnumeric(tsegsize) || ~isfinite(tsegsize)
    error('Missing required input argument: tsegsize must be a numeric finite scalar');
end

valid_time_max = size(rawdata,2) / Fs;
valid_time_min = 0;

if trange(1) < valid_time_min || trange(2) > valid_time_max || trange(2) <= trange(1)
    error('Specified trange [%g %g] is out of valid range [%g %g]', trange(1), trange(2), valid_time_min, valid_time_max);
end

if tsegsize > trange(2) - trange(1)
    error('Specified tsize [%g] is too long', tsegsize);
end

if tsegsize < 1/Fs
    error('Specified tsize [%g] is too short', tsegsize);
end

segpow = cell(1,size(rawdata,3));
seglabels = cell(1,length(labels));
for e = 1:size(rawdata,3)
    segpow{e} = convert_timeseriesdata_to_segmentedpowervalues(rawdata(:,:,e).', Fs, frange, trange, tsegsize, opts);
    seglabels{e} = repmat(labels(e),[1,size(segpow{e},3)]);
end
segpowdata = cat(3,segpow{:});
seglabels = cat(2,seglabels{:});





