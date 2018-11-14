% DCT-based signal decompression, unpacked from uint8 array per channel

function signal = dctdecompress (packed)
if size(packed,1) == 1
    packed = packed.';
end

for ch = size(packed,2):-1:1
    N = typecast(packed(1:4,ch), 'uint32');
    N8 = ceil(N/8)*8;
    FMIN = typecast(packed(5:8,ch), 'single');
    FDYN = typecast(packed(9:12,ch), 'single');
    BMAP8 = packed(13:12+N8/8,ch).';
    payload = packed(13+N8/8:end,ch).';
    signal(:,ch) = dctdecompress_1(payload, BMAP8, FMIN, FDYN, N);
end



function signal = dctdecompress_1 (payload, BMAP8, FMIN, FDYN, N)
dc2b = double(payload) / 255 * FDYN + FMIN;
BMAP = logical(uint8(dec2bin(BMAP8(:), 8)) - '0');
BMAP = reshape(BMAP.',[],1);
dc2 = zeros(length(BMAP),1);
dc2(BMAP) = dc2b;
signal = fastidct(dc2, N);
return

