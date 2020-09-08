% Convert from MATLAB's datenum (local time) to Unixtime (always in UTC)
function u = datenum2unixtime (n)
u = (n - now)*86400 + unixtime;
