% Convert from MATLAB's datenum (local time) to Unixtime (always in UTC)
function u = datenum2unixtime (n)
dt = datetime(n, 'ConvertFrom', 'datenum');
dt.TimeZone = 'local';
u = datetime2unixtime(dt);
