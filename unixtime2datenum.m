% Convert from Unixtime (always in UTC) to MATLAB's datenum (local time)
function n = unixtime2datenum (u)
n = (u - unixtime)/86400 + now;
