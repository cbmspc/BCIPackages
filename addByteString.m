% Automatically convert bytes to the appropriate binary multiple,
% optionally with sprintf format specifier
function str = addByteString(bytes, formatSpecifier)
multipliers = [
2^0
2^10
2^20
2^30
2^40
2^50
2^60
2^70
2^80
];

multiplier_suffix = {
'B'
'KiB'
'MiB'
'GiB'
'TiB'
'PiB'
'EiB'
'ZiB'
'YiB'
};

if ~exist('formatSpecifier','var') || isempty(formatSpecifier) || ~ischar(formatSpecifier) || ~contains(formatSpecifier,'%')
    formatSpecifier = '%.3g';
end

suffix = '';
value = bytes;

for i = length(multipliers):-1:1
    if bytes >= multipliers(i)
        suffix = multiplier_suffix{i};
        value = bytes / multipliers(i);
        break
    end
end

str = sprintf([formatSpecifier ' %s'], value, suffix);