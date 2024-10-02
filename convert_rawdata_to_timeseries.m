% timeseries = convert_rawdata_to_timeseries (rawdata)
% Rawdata format is (chan, time, epoch)
% Time series is (time, chan)
% This function stacks the epochs in time in the same order as they appear
% in rawdata, and then transposes so that the output ends up in the 
% (time, chan) orientation.

function timeseries = convert_rawdata_to_timeseries (rawdata)
NumChan = size(rawdata,1);
NumEpochTime = size(rawdata,2);
NumEpoch = size(rawdata,3);

timeseries = nan(NumEpochTime*NumEpoch, NumChan);

for ep = 1:NumEpoch
    timeseries((ep-1)*NumEpochTime + (1:NumEpochTime), :) = rawdata(:,:,ep).';
end

