% Convert from TAI time (seconds) to Unixtime (seconds)
% TAI clock predates unixtime by 12 years (378691200 seconds).
function u = tai2unixtime (t)
u = t - 378691200;