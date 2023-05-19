% Time zone offset from UTC, in hours (accounts for daylight saving).
function h = tzoffsethour()
h = round((unixtime2datenum(0) - 719529)*1440)/60;
