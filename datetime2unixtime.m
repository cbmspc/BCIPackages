% First make sure the TimeZone property in the datetime object is set properly
% For example, D = datetime(sprintf('%s-%s-%s %s:%s:%s', b{1}, b{2}, b{3}, b{4}, b{5}, b{6}), 'TimeZone', 'Etc/UTC');

function U = datetime2unixtime(D)
D2 = datetime(0, 'ConvertFrom', 'epochtime', 'TimeZone', 'Etc/UTC');
if isempty(D.TimeZone)
    D.TimeZone = 'local';
end
Dur = D - D2;
U = seconds(Dur);


