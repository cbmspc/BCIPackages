% Converts bit array of single-precision floating-point format (IEEE 754)
% to Matlab variable. Based on
% https://en.wikipedia.org/wiki/Single-precision_floating-point_format .
% Bit array is an array of 1 or 0 of length 32, read like a human would
% example, 3.1416 would be 01000000010010010000111111010000
% and b = '01000000010010010000111111010000' - '0';

function d = bitarray2float (b)
d = (-1)^b(1)*(1+b(10:32)*2.^-(1:23)')*2^(b(2:9)*2.^(7:-1:0)'-127);

