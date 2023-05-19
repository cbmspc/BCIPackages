% TAIM (TAI time in minutes) to MATLAB datenum (days)
function n = taim2datenum (t, nonutc)
u = tai2unixtime(t * 60);
n = unixtime2datenum(u);
if exist('nonutc', 'var') && nonutc
    n = n - tzoffsethour()/24;
end
