% TAIM (TAI time in minutes) to MATLAB datenum (days)
% The use of datenum is not recommended for high-precision comparisons due
% to the precision limit of the "double" data type
% 
function n = taim2datenum (t, nonutc)
if ischar(t)
    % in the form tttttttt:fffff
    t = str2double(split(t, ':'));
    if numel(t) < 2
        t(2) = 0;
    end
    t = t(1) + t(2)/60000;
end
u = tai2unixtime(t * 60);
n = unixtime2datenum(u);
if exist('nonutc', 'var') && nonutc
    n = n - tzoffsethour()/24;
end
