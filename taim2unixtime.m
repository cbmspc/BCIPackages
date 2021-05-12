% TAIM (TAI time in minutes) to unixtime (seconds)
function u = taim2unixtime (t)
u = tai2unixtime(t * 60);