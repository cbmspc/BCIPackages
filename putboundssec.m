% Opposite of getbounds, in seconds
function Signal = putboundssec (BoundsSec, Fs, SigLen)
Bounds = round(BoundsSec*Fs+1);
Bounds(:,2) = Bounds(:,2) - 1;
Signal = putbounds(Bounds, SigLen);