% Converts little endian byte array to float. 4 bytes make a float.
function F = bytearray2float (B)
F = nan(1,floor(length(B)/4));
for i = 0:floor(length(B)/4)-1
    F(i+1) = bitarray2float(reshape(dec2bin(B(i*4+(1:4)),8).',1,[]) - '0', 1);
end