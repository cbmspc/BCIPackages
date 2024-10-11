% rawdata = convert_timeseries_to_rawdata (timeseries, NumEpoch)
% Time series is (time, chan)
% Rawdata format is (chan, time, epoch)
% This function reverses the action done by convert_rawdata_to_timeseries
%

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

