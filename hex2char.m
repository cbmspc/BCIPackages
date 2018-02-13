function str = hex2char (hex)
% hex must be in two-digit per byte format

strlength = length(hex)/2;
if strlength < 1
    str = '';
    return
end
for i = 1:strlength
    str(i) = char(hex2dec(hex((i-1)*2+[1:2])));
end

