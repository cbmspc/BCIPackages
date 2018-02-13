% Converts bit array of single-precision floating-point format (IEEE 754)
% to Matlab variable. Based on
% https://en.wikipedia.org/wiki/Single-precision_floating-point_format .
% Bit array is an array of 1 or 0 of length 32, read like a human would
% example, 3.1416 would be 01000000010010010000111111010000
% and b = '01000000010010010000111111010000' - '0';
%
% If IsByteLittleEndian = 1, processes the byte array as little endian,
% i.e. 3.1416 would be '11010000000011110100100101000000' - '0' as stored
% in the memory of a modern PC
% 

function d = bitarray2float (b, IsByteLittleEndian)
if ~exist('IsByteLittleEndian','var')
    IsByteLittleEndian = 0;
end
d = nan(1,floor(length(b)/32));
for i = 0:floor(length(b)/32)-1
    d(i+1) = bitarray2float_once(b(i*32+(1:32)), IsByteLittleEndian);
end


function d = bitarray2float_once (b, IsByteLittleEndian)
if IsByteLittleEndian
    tmp = b;
    tmp(1:8) = b(25:32);
    tmp(9:16) = b(17:24);
    tmp(17:24) = b(9:16);
    tmp(25:32) = b(1:8);
    b = tmp;
end
d = (-1)^b(1)*(1+b(10:32)*2.^-(1:23)')*2^(b(2:9)*2.^(7:-1:0)'-127);

