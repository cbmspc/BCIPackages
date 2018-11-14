% DCT-based signal compression, packed into uint8 array per channel
% signal = (time, chan)
% KF = compression factor. Keeps N/KF coefficients.

function packed = dctcompress (signal, KF)
if size(signal,1) == 1
    signal = signal.';
end
[dc1, N] = fastdct(signal);
N8 = ceil(N/8)*8;
N = uint32(N);

NK = floor(N8 / KF);

for ch = size(signal,2):-1:1
    [FMIN, FDYN, BMAP8, payload] = dctcompress_1(dc1(:,ch), NK);
    
    packed(:,ch) = [
        typecast(N, 'uint8').';
        typecast(FMIN, 'uint8').';
        typecast(FDYN, 'uint8').';
        BMAP8;
        payload
    ];

end



function [FMIN_low, FDYN_low, BMAP8, payload] = dctcompress_1 (dc1, NK)

[~, ind] = sort(abs(dc1), 'descend');
NK = floor(NK);

dc2 = dc1*0;
dc2(ind(1:NK)) = dc1(ind(1:NK));

BMAP = false(size(dc1,1),1);
BMAP(ind(1:NK)) = 1;
BMAP8 = uint8(reshape(BMAP,8,[]).'*[128 64 32 16 8 4 2 1].');

NSTEP = 256;

FMAX = max(dc2);
FMIN = min(dc2);
FDYN = FMAX - FMIN;

FDYN_low = single(FDYN);
FMIN_low = single(FMIN);

dc2b = round((dc2 - FMIN) / FDYN * (NSTEP-1));
payload = uint8(dc2b(BMAP));

return
