% TAIM (TAI time in minutes) to MATLAB datenum (days)
function n = taim2datenum (t)
u = tai2unixtime(t * 60);
n = unixtime2datenum(u);
