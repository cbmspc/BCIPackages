% Be aware that the output DateTime object has TimeZone set to UTC
function D = unixtime2datetime(U)
D = datetime(U, 'ConvertFrom', 'epochtime', 'TimeZone', 'Etc/UTC');


