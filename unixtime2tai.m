% Convert from Unixtime (seconds) to TAI time (seconds)
% TAI clock predates unixtime by 12 years (378691200 seconds).
function t = unixtime2tai (u)
t = u + 378691200;