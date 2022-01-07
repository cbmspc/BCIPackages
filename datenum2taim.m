% MATLAB datenum (days) to TAIM (TAI time in minutes)
function t = datenum2taim (n)
u = datenum2unixtime(n);
t = unixtime2taim(u);
