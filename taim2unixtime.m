% TAIM (TAI time in minutes) to unixtime (seconds)
% Reminder: TAI, TAIM, and Unixtime all do not care about timezones or
% daylight saving time. 
% TAI (atomic clock time) ignores leap seconds, whereas Unixtime repeats
% them. Unixtime and UTC are exactly correlated since 1972-01-01.
%
function u = taim2unixtime (t)
if iscell(t)
    u = cellfun(@taim2unixtime,t);
    return
end
if ischar(t)
    % in the form tttttttt:fffff
    t = str2double(split(t, ':'));
    if numel(t) < 2
        t(2) = 0;
    end
    t = t(1) + t(2)/60000;
end
u = tai2unixtime(t * 60);