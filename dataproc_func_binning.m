% binning.m
% Syntax: BinData = binning(Data, binsize)
% Data format is (chan#, sample points, trial#)
% binsize is in sample points unit
% BinData is the binned data with all channels column-concatenated
% binning([1 2 3; 4 5 6], 1) ==> [1 4 2 5 3 6]'
%
% Data = windowed data or any raw data, as long as data is in the form of 
% (chan#, sample#, trial#).
%
% The unit of binsize is number of samples. So if Fs = 200 Hz (T = 5 ms),
% a 15 ms bin is equal to binsize = 3.
%
% This is a brand new binning code, not compatible with the old one

function BinData = binning (Data, binsize)

Ntrials = size(Data, 3);
Nsamples = size(Data, 2);
Nchan = size(Data, 1);

BinData = zeros(floor(Nsamples/binsize)*Nchan,Ntrials);

for j = 1:Ntrials
    i = 0;
    for a = 1:binsize:Nsamples
        b = a + binsize - 1;
        if b > Nsamples
            break
        end
        i = i + 1;
        TimeBin(:,i) = mean(Data(:, a:b, j), 2);
    end
    BinData(:,j) = reshape(TimeBin, 1, []);
end

