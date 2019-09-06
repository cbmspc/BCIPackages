% 1) Runs a band-pass filter on a continuous uncut time series data
% 2) Segment the filtered data into equal sized segments within the
%    specified time range   
% (1 and 2 simulate a filter bank that does not reset)
% 3) Calculate the power within each segment. There is one power value per
%    channel per segment
% 4) Output format: chan x band x segment
% Note: All segments have equal size. Possible gap between last and tend.

function segpow = convert_timeseriesdata_to_segmentedpowervalues (data, Fs, frange, trange, tsegsize, filtorder)

if ~exist('filtorder','var') || isempty(filtorder)
    filtorder = 4;
end

filtname = 'buttercausal';

tstart = trange(1);
tend = trange(2);

if tend > size(data,1)/Fs
    tend = size(data,1)/Fs;
end

data = detrend(data, 'constant');

for band = 1:size(frange,1)
    


    
%     % Segment then filter
%     nseg = floor((tend-tstart)/tsegsize);
%     for seg = 1:nseg
%         ta = tstart + (seg-1)*tsegsize;
%         tb = ta + tsegsize;
%         ka = round(ta*Fs+1);
%         kb = round(tb*Fs);
%         segdata = data(ka:kb,:);
%         
%         fb = [1 1];
%         if frange(band,1) <= 0
%             fb(1) = 0;
%         end
%         if frange(band,2) >= Fs/2
%             fb(2) = 0;
%         end
%         if fb(1) && fb(2)
%             filtparm = [frange(band,1:2) filtorder];
%             filttype = 'pass';
%             segbanddata = freqfilter(segdata, Fs, filtparm, filttype, filtname);
%         elseif fb(1) && ~fb(2)
%             filtparm = [frange(band,1) filtorder];
%             filttype = 'high';
%             segbanddata = freqfilter(segdata, Fs, filtparm, filttype, filtname);
%         elseif ~fb(1) && fb(2)
%             filtparm = [frange(band,2) filtorder];
%             filttype = 'low';
%             segbanddata = freqfilter(segdata, Fs, filtparm, filttype, filtname);
%         else
%             segbanddata = segdata;
%         end
%         segpow(:,band,seg) = mean(segbanddata.^2,1) / (frange(band,2)-frange(band,1));
%     end

    
    
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
        banddata = freqfilter(data, Fs, filtparm, filttype, filtname);
    elseif fb(1) && ~fb(2)
        filtparm = [frange(band,1) filtorder];
        filttype = 'high';
        banddata = freqfilter(data, Fs, filtparm, filttype, filtname);
    elseif ~fb(1) && fb(2)
        filtparm = [frange(band,2) filtorder];
        filttype = 'low';
        banddata = freqfilter(data, Fs, filtparm, filttype, filtname);
    else
        banddata = data;
    end
    
    % Segment
    nseg = floor((tend-tstart)/tsegsize);
    for seg = 1:nseg
        ta = tstart + (seg-1)*tsegsize;
        tb = ta + tsegsize;
        ka = round(ta*Fs+1);
        kb = round(tb*Fs);
        segbanddata = banddata(ka:kb,:);
        segpow(:,band,seg) = mean(segbanddata.^2,1) / (frange(band,2)-frange(band,1));
    end
    
    
end