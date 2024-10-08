% timeseries = convert_rawdata_to_timeseries (rawdata)
% Rawdata format is (chan, time, epoch)
% Time series is (time, chan)
% This function stacks the epochs in time in the same order as they appear
% in rawdata, and then transposes so that the output ends up in the 
% (time, chan) orientation.

function rawdata = convert_timeseries_to_rawdata (timeseries, NumEpoch)
NumChan = size(timeseries,2);
NumTotalTime = size(timeseries,1);
NumEpochTime = NumTotalTime / NumEpoch;

if rem(NumEpoch,1) ~= 0 || NumEpoch < 0
    error('Number of epochs must be a positive integer.');
end
if rem(NumEpochTime,1) ~= 0 || NumEpochTime < 0
    error('The total number of time samples must be an integer multiple of NumEpoch.');
end

rawdata = nan(NumChan, NumEpochTime, NumEpoch);

for ep = 1:NumEpoch
    rawdata(:,:,ep) = timeseries((ep-1)*NumEpochTime + (1:NumEpochTime), :).';
end

